package fr.inserm.umr1087.jvarkit.tools.cnv20120907;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;

import org.apache.commons.math.analysis.interpolation.LoessInterpolator;

import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

import fr.inserm.umr1087.jvarkit.segments.TidStartEnd;
import fr.inserm.umr1087.jvarkit.segments.bdb.TidStartEndBinding;
import fr.inserm.umr1087.jvarkit.segments.bdb.TidStartEndSorter;
import fr.inserm.umr1087.jvarkit.tools.genomegraph.Graphics2DGenomeGraphDrawer;
import fr.inserm.umr1087.jvarkit.util.intervalparser.IntervalParser;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class CNV20120927
	{
	private int BUFFER_SIZE=1000000;//1E6
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
	private int minQual=30;
	private TidStartEndBinding tidStartEndBinding=new TidStartEndBinding();	
	
	private enum Step {
		CALC_LOESS,
		APPLY_LOESS_AND_GET_MIN,
		SUBSTRACT_MIN_AND_GET_MEDIAN,
		DIV_MEDIAN,
		NORMALIZE_MEDIAN,
		RUN_MED
		}
	
	
	private ReferenceSequenceFile reference=null;
	private int windowSize=100;
	private int windowStep=50;
	private List<BamBuffer> bams=new ArrayList<BamBuffer>();
	private Interval targetInterval=null;
	private File dbHome=null;
	private Environment environment=null;
	private Database region2depths;
	private Transaction txn=null;
	
	private class BamBuffer
		{
		File file;
		BamBuffer(File file)
			{
			//super(new SAMFileReader(file));
			
			this.file=file;
			}
		}
	
	private static class Spline
		{
		double x_val[];
		double y_val[];
		Spline(int n)
			{
			x_val=new double[n];
			y_val=new double[n];
			}
		private int lowerBound(double L)
			{
			int first=0;
		    int len = x_val.length;
		    while (len > 0)
		            {
		            int half = len / 2;
		            int middle = first + half;
		
		            if ( x_val[middle]< L  )
		                    {
		                    first = middle + 1;
		                    len -= half + 1;
		                    }
		            else
		                    {
		                    len = half;
		                    }
		            }
		    return first;
			}

		double value(double x)
			{
			int n=lowerBound(x);
			if(n>=y_val.length) return y_val[y_val.length-1];
			return y_val[n];
			}	
		}


	private CNV20120927()
		{
		
		}

	
	private void open()
		{
		String dbName1="region2depths";
		EnvironmentConfig envConfig= new EnvironmentConfig();
		envConfig.setAllowCreate(true);
		envConfig.setReadOnly(false);
		envConfig.setConfigParam(EnvironmentConfig.LOG_FILE_MAX,"250000000");
		envConfig.setTransactional(true);
		this.environment= new Environment(dbHome, envConfig);
		this.txn=this.environment.beginTransaction(null, null);
		for(String db:this.environment.getDatabaseNames())
			{
			if(db.equals(dbName1)) this.environment.removeDatabase(txn, db);
			}
		DatabaseConfig cfg= new DatabaseConfig();
		cfg.setAllowCreate(true);
		cfg.setReadOnly(false);
		cfg.setTransactional(true);
		cfg.setBtreeComparator(TidStartEndSorter.class);
		this.region2depths= this.environment.openDatabase(txn,dbName1,cfg);
		}
	
	
	
	private void close()
		{
		if(this.region2depths!=null)
			{
			try{region2depths.close();} catch(Exception err) {}
			this.region2depths=null;
			}
		if(this.txn!=null)
			{
			try{txn.commit();} catch(Exception err) {}
			this.txn=null;
			}
		if(this.environment!=null)
			{
			try{environment.close();} catch(Exception err) {}
			this.environment=null;
			}
		}
	
	
	private void dump(File file) throws IOException
		{
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		Cursor c=this.region2depths.openCursor(txn, null);
		PrintWriter pw=new PrintWriter(file);
		pw.print("CHROM\tSTART\tEND\tGC");
		for(int bx=0;
		bx< this.bams.size();
		++bx)
			{
			pw.print("\t"+this.bams.get(bx).file.getName().replace(".bam", ""));
			}
		pw.println();
		while(c.getNext(key, data, null)==OperationStatus.SUCCESS)
			{
			TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
			CNVRow row=CNVRow.BINDING.entryToObject(data);
			
			pw.print(
					this.reference.getSequenceDictionary().getSequence(seg.getChromId()).getSequenceName()+
					"\t"+seg.getStart()+
					"\t"+seg.getEnd()+
					"\t"+row.getGcPercent()
					);
			for(int i=0;i< this.bams.size();++i)
				{
				pw.print("\t"+row.getDepth(i));
				}
			pw.println();
				
			}
		c.close();
		pw.flush();
		pw.close();
		}
	
	private void dumpGnuplot(File file) throws IOException
		{
		
		ZipOutputStream zout=new ZipOutputStream(new FileOutputStream(file));
		
		for(int bx=0;
			bx< this.bams.size();
			++bx)
			{
			long shift=0L;
			ZipEntry zentry=new ZipEntry(this.bams.get(bx).file.getName()+".tsv");
			zout.putNextEntry(zentry);
			PrintStream pout=new PrintStream(zout);
			for(SAMSequenceRecord chrom: this.reference.getSequenceDictionary().getSequences())
				{
				if(targetInterval!=null && !chrom.getSequenceName().equals(targetInterval.getSequence())) continue;
				long chromMax=0L;
				DatabaseEntry key=new DatabaseEntry();
				DatabaseEntry data =new DatabaseEntry();
				Cursor c=this.region2depths.openCursor(txn, null);
				while(c.getNext(key, data, null)==OperationStatus.SUCCESS)
					{
					
					TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
					if(seg.getChromId()!=chrom.getSequenceIndex())
						{
						continue;
						}
					int chromStart=seg.getStart();				
					CNVRow row=CNVRow.BINDING.entryToObject(data);
					if(row.containsNanDepth()) continue;
					chromMax=Math.max(chromMax, seg.getEnd());
				    double y=row.getDepth(bx);
					pout.println(""+(shift+chromStart)+"\t"+y);
					}
				c.close();
				shift+=chromMax;
				}
			pout.flush();
			zout.closeEntry();
			}
		
		
		ZipEntry zentry=new ZipEntry("jeter.gnuplot");
		zout.putNextEntry(zentry);
		PrintStream gnuplotout=new PrintStream(zout);
		gnuplotout.println("set term postscript");
		gnuplotout.println("set output jeter.ps");
		
		for(int bx=0;
		bx< this.bams.size();
		++bx)
			{
			gnuplotout.println("set title \"SAMPLE-"+(1+bx)+" "+this.bams.get(bx).file.getName()+"\"");
			gnuplotout.println("set xlabel \"Position\"");
			gnuplotout.println("set ylabel \"Depth\"");
			gnuplotout.println("set yrange [-3:3]");
			gnuplotout.println("plot \""+ this.bams.get(bx).file.getName()+".tsv\" using 1:2 notitle\n");
			}
		gnuplotout.flush();
		zout.closeEntry();
		
		zout.finish();
		zout.close();
		}
	
	private void run() throws Exception
		{
		open();
		LoessInterpolator loessInterpolator=new LoessInterpolator();
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		SAMSequenceDictionary 	dict=this.reference.getSequenceDictionary();
		/** loop over each chromosome */
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			if(targetInterval!=null && !targetInterval.getSequence().equals(chrom.getSequenceName()))
				{
				continue;
				}
			/** loop over the bases of the chromosomes */
			byte bases[]=this.reference.getSequence(chrom.getSequenceName()).getBases();
			int  start0=(targetInterval!=null?targetInterval.getStart():0);
			while(start0+this.windowSize<
					(targetInterval!=null?targetInterval.getEnd():bases.length))
				{
				int gc=0;
				int N=0;
				for(int i=start0;N==0 && i<start0+this.windowSize;++i)
					{
					byte base=bases[i];
					switch(base)
						{
						case 'n': case 'N': ++N;break;
						case 's': case 'S':
						case 'c': case 'C':
						case 'g': case 'G': ++gc; break;
						default:break;
						}
					}
				if(N>0)
					{
					/* if there is a 'N' in the window, skip to next base */
					++start0;
					continue;
					}
				TidStartEnd tidStartEnd=new TidStartEnd(
						chrom.getSequenceIndex(),
						start0,
						start0+windowSize
						);
				CNVRow row=new CNVRow(this.bams.size());
				row.setGcPercent(gc/(double)windowSize);
				
				
				
				
				tidStartEndBinding.objectToEntry(tidStartEnd, key);
				CNVRow.BINDING.objectToEntry(row, data);
				
				if(this.region2depths.put(txn, key, data)!=OperationStatus.SUCCESS)
					{
					close();
					throw new RuntimeException("Cannot insert");
					}

				start0+=this.windowStep;
				}
			bases=null;
			
			/* get the depth for all bams */
			short depth0[]=new short[chrom.getSequenceLength()];
			for(int bx=0;
				bx< this.bams.size();
				++bx)
					{
					LOG.info("Getting depth for "+chrom.getSequenceName()+" "+this.bams.get(bx).file);
					Arrays.fill(depth0, (short)0);
					SAMFileReader sf=new SAMFileReader(this.bams.get(bx).file);
					sf.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
					
					/* create a SamLocusIterator to get the DEPTH at each position */
					Interval interval;
					if(targetInterval==null)
						{
						interval=new Interval(
			            		chrom.getSequenceName(),
			            		1,
			            		depth0.length
			            		);//Coordinates are 1-based closed ended.
						}
					else
						{
						interval=new Interval(
			            		chrom.getSequenceName(),
			            		Math.max(1, this.targetInterval.getStart()),
			            		Math.min(this.targetInterval.getEnd(),depth0.length)
			            		);//Coordinates are 1-based closed ended.
						}
					IntervalList iL=new IntervalList(sf.getFileHeader());
		            iL.add(interval);
		            SamLocusIterator sli=new SamLocusIterator(sf,iL,true);           
		            sli.setMappingQualityScoreCutoff(this.minQual);
		            for(Iterator<SamLocusIterator.LocusInfo>  iter=sli.iterator();
		                iter.hasNext();
		                )
		                {
		                SamLocusIterator.LocusInfo locusInfo=iter.next();
		                int heredepth=locusInfo.getRecordAndPositions().size();
		                if(heredepth>Short.MAX_VALUE)
		                	{
		                	LOG.info("WARNING depth >"+Short.MAX_VALUE+" at "+locusInfo.getPosition());
		                	}
		                depth0[locusInfo.getPosition()-1]=(short)Math.min(heredepth, Short.MAX_VALUE);
		                }
		            
		            sf.close();
		            
		            /* calc the mean depth at each window */
					key=new DatabaseEntry();
					data=new DatabaseEntry();
					Cursor c=this.region2depths.openCursor(txn, null);
					boolean first=false;
					
					for(;;)
						{
						OperationStatus status;
						if(first) 
							{
							first=false;
							TidStartEnd init=new TidStartEnd(chrom.getSequenceIndex(),
										0,
										0
										);
							
							this.tidStartEndBinding.objectToEntry(init, key);
							status=c.getSearchKeyRange(key, data,null);
							LOG.info("starting with "+this.tidStartEndBinding.entryToObject(key));
							}
						else
							{
							status=c.getNext(key, data,null);
							}
						if(status!=OperationStatus.SUCCESS) break;
						TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
						if(seg.getChromId()<chrom.getSequenceIndex()) continue;
						if(seg.getChromId()>chrom.getSequenceIndex()) break;
						
						int count=0;
						double sum=0.0;
						for(int i=seg.getStart();i<seg.getEnd() &&
								i < depth0.length;
								++i)
							{
							sum+=depth0[i];
							++count;
							}
						
						CNVRow row=CNVRow.BINDING.entryToObject(data);
						row.setDepth(bx, (count==0?Double.NaN:sum/count));
						CNVRow.BINDING.objectToEntry(row,data);
						if(c.putCurrent(data)!=OperationStatus.SUCCESS)
							{
							c.close();
							close();
							throw new RuntimeException("Cannot update "+row);
							}
						}
					
					c.close();
					}
			depth0=null;
			}
		
		dump(new File("/commun/data/users/lindenb/_ignore.backup/jeter.tsv"));
		
		Cursor c=this.region2depths.openCursor(txn, null);
		List<Point2D.Double> points=null;
		OperationStatus status;
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			if(targetInterval!=null && !chrom.getSequenceName().equals(targetInterval.getSequence())) continue;
			
			for(int bx=0;
					bx< this.bams.size();
					++bx)
				{
			
				Spline splineFun=null; 
				double minDepth= Integer.MAX_VALUE;
				List<Double> vData=null;
				double medianValue=Double.NaN;
				
				/* do each step for the BAM */
				for(Step step: Step.values())
					{
					LOG.info("step "+step+" "+this.bams.get(bx).file.getName());
					key=new DatabaseEntry();
					data=new DatabaseEntry();
					int rowIndex=-1;
					boolean first=true;
					points=null;
					
					/* initialize this step */
					switch(step)
						{
						case APPLY_LOESS_AND_GET_MIN:
							{
							minDepth= Double.MAX_VALUE;
							break;
							}
						case CALC_LOESS:
							{
							points=new ArrayList<Point2D.Double>();
							break;
							}
						case SUBSTRACT_MIN_AND_GET_MEDIAN:
						case NORMALIZE_MEDIAN:
							{
							vData=new ArrayList<Double>();
							break;
							}
						}
					
					/* loop over the key for this chromosome */
					for(;;)
						{
					
						if(first) 
							{
							first=false;
							TidStartEnd init=new TidStartEnd(chrom.getSequenceIndex(), 0, 0);
							this.tidStartEndBinding.objectToEntry(init, key);
							status=c.getSearchKeyRange(key, data,null);
							}
						else
							{
							status=c.getNext(key, data,null);
							}
						if(status!=OperationStatus.SUCCESS) break;
						TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
						if(seg.getChromId()<chrom.getSequenceIndex()) continue;
						if(seg.getChromId()>chrom.getSequenceIndex()) break;
						++rowIndex;//initialized to '-1'
						
						CNVRow row=CNVRow.BINDING.entryToObject(data);
						
						
						
						switch(step)
							{
							case CALC_LOESS:
								{
								points.add(new Point2D.Double(
									row.getDepth(bx),
									row.getGcPercent()
									));
								break;
								}
							case APPLY_LOESS_AND_GET_MIN:
							case SUBSTRACT_MIN_AND_GET_MEDIAN:
							case DIV_MEDIAN:
							case RUN_MED:
								{
								double newDepth=0;
								switch(step)
									{
									case APPLY_LOESS_AND_GET_MIN:
										{
										newDepth=splineFun.value(row.getDepth(bx));
										if(minDepth> newDepth) minDepth=newDepth;
										break;
										}
									case SUBSTRACT_MIN_AND_GET_MEDIAN:
										{	
										newDepth= row.getDepth(bx)- minDepth;
										vData.add(newDepth);
										break;
										}
									case DIV_MEDIAN:
										{
										newDepth= row.getDepth(bx) / medianValue;
										break;
										}
									case RUN_MED:
										{
										if(vData.get(rowIndex)!=row.getDepth(bx))
											{
											throw new RuntimeException("Err "+vData.get(rowIndex)+"/"+row.getDepth(bx));
											}
										List<Double> runmed=new ArrayList<Double>(11);
										for(int i=Math.max(0,rowIndex-5);
												i<= rowIndex+5 && i< vData.size();
												++i)
											{
											runmed.add(vData.get(i));
											}
										Collections.sort(runmed);
										newDepth=Math.log(runmed.get(runmed.size()/2));
										if(Double.isNaN(newDepth)) newDepth=0;
										break;
										}
									default: throw new RuntimeException();
									}
								row.setDepth(bx,newDepth);
								CNVRow.BINDING.objectToEntry(row,data);
								if(c.putCurrent(data)!=OperationStatus.SUCCESS)
									{
									c.close();
									close();
									throw new RuntimeException("Cannot update "+row);
									}
								break;
								}
							case NORMALIZE_MEDIAN:
								{
								if(bx==0)
									{
									double v[]=new double[this.bams.size()-1];
									for(int i=1;i<this.bams.size();++i)
										{
										v[i-1]=row.getDepth(i);
										}
									Arrays.sort(v);
									double mediane=v[v.length/2];
									if(mediane>0.2)
										{
										for(int i=0;i< this.bams.size();++i)
											{
											row.setDepth(i, row.getDepth(i)/mediane);
											}
										CNVRow.BINDING.objectToEntry(row,data);
										if(c.putCurrent(data)!=OperationStatus.SUCCESS)
											{
											c.close();
											close();
											throw new RuntimeException("Cannot update "+row);
											}
										}
									else
										{
										c.delete();
										break;//don't add value to vData below
										}
									}
								vData.add(row.getDepth(bx));
								break;
								}
							default: throw new RuntimeException(step.toString());
							}
						}
					
					
					if(step==Step.CALC_LOESS && points.isEmpty()) break;
					
					/* post step processing */
					switch(step)
						{
						case CALC_LOESS:
							{
							
							Collections.sort(points,new Comparator<Point2D.Double>()
									{
									@Override
									public int compare(Point2D.Double o1, Point2D.Double o2)
										{
										if(o1.getX()==o2.getX())
											{
											if(o1.getY()==o2.getY()) return 0;
											return o1.getY() < o2.getY() ? -1  : 1 ;
											}
										return o1.getX() < o2.getX() ? -1  : 1 ;
										}
									});
							splineFun= new Spline(points.size());
							
							for(int i=0;i< points.size();++i)
								{
								splineFun.x_val[i]=points.get(i).x;
								splineFun.y_val[i]=points.get(i).y;
								
								if(i>0 && splineFun.x_val[i]<=splineFun.x_val[i-1])//loess plante si xi==x(i-1)
									{
									splineFun.x_val[i]=splineFun.x_val[i-1]+1E-9;
									}
								}
							points=null;
							splineFun.y_val=loessInterpolator.smooth(splineFun.x_val, splineFun.y_val);
							break;
							}
						case APPLY_LOESS_AND_GET_MIN:
							{
							splineFun=null;
							break;
							}
						case SUBSTRACT_MIN_AND_GET_MEDIAN:
							{
							Collections.sort(vData);
							medianValue = vData.get(vData.size()/2);
							vData=null;
							break;
							}
						case RUN_MED:
							{
							vData=null;
							break;
							}
						}	//end-switch
					if(bx+1==bams.size())
						{
						dump(new File("/commun/data/users/lindenb/_ignore.backup/jeter."+step+".tsv"));
						}
					}
				
				
				}
			}
		
	
		
		
		for(int bx=0;
			bx< this.bams.size();
			++bx)
			{
			Graphics2DGenomeGraphDrawer drawer=new Graphics2DGenomeGraphDrawer(this.reference.getSequenceDictionary(),this.targetInterval);
			drawer.setMinY(-3.0);
			drawer.setMaxY(3.0);
			BufferedImage img=drawer.createImage();
			Graphics2D g=(Graphics2D)img.getGraphics();
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			drawer.init(g);
			
			File imgFile=new File("/commun/data/users/lindenb/_ignore.backup/"+this.bams.get(bx).file.getName()+".jpg");
			key=new DatabaseEntry();
			status=c.getFirst(key, data, null);
			while(status==OperationStatus.SUCCESS)
				{
				
				TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
				int tid=seg.getChromId();
				int chromStart=seg.getStart();
				int chromEnd=seg.getEnd();
				
				CNVRow row=CNVRow.BINDING.entryToObject(data);
				if(!row.containsNanDepth())
					{
				    double y=row.getDepth(bx);
				    
				    double pix_y= drawer.convertValueToPixelY(y);
				    double pix_x= drawer.convertPositionToPixelX(tid, (chromStart+chromEnd)/2);
				    g.setColor(Color.BLACK);
				    g.fillOval((int)pix_x-2,(int)pix_y-2,5,5);
					status=c.getNext(key, data, null);
					}
				}
			drawer.finish();
			g.dispose();
			ImageIO.write(img, "JPG", imgFile);
			
			}
		c.close();
		
		//dumpGnuplot(new File("/commun/data/users/lindenb/_ignore.backup/jeter.zip"));
		
		}
	
	private void usage()
		{
		System.err.println(
			"Options:\n"+
			" --db-home (dir) berkeleydb-home\n"+
			" -f (fasta) reference fasta indexed with faidx\n"+
			" -q <int> quality treshold for BAM alignments\n"+
			" -L <Chrom:star-end> target interval (optional)\n"+
			" -w <int> window-size["+this.windowSize+"]\n"+
			" -s <int> step-size["+this.windowStep+"]\n"
			);
		}

	private void run(String[] args) throws Exception
		{
		
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				usage();
				return;
				}
			else if(args[optind].equals("--db-home") && optind+1 < args.length )
				{
				this.dbHome=new File(args[++optind]);
				}
			else if(args[optind].equals("-f") && optind+1 < args.length )
				{
				this.reference= ReferenceSequenceFileFactory.getReferenceSequenceFile(
						new File(args[++optind]),
						true //cut after first whitespace
						);
				}
			else if(args[optind].equals("-q") && optind+1 < args.length )
				{
				this.minQual=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-w") && optind+1 < args.length )
				{
				this.windowSize=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-s") && optind+1 < args.length )
				{
				this.windowStep=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-L") && optind+1 < args.length )
				{
				this.targetInterval=IntervalParser.parseOne(args[++optind]);
				if(this.targetInterval==null)
					{
					System.err.println("bad range:"+args[optind]);
					return;
					}
				//int len=this.targetInterval.getEnd()-this.targetInterval.getStart();
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return;
				}
			else
				{
				break;
				}
			++optind;
			}
		if(BUFFER_SIZE<=this.windowSize)
			{
			BUFFER_SIZE=(this.windowSize+1000);
			}
		
		if(this.dbHome==null)
			{
			System.err.println("--db-home missing");
			return;
			}
		
		if(this.reference==null)
			{
			System.err.println("Reference missing");
			return;
			}
		//this.referenceBuffer=new ReferenceBuffer(this.reference);
		//this.referenceBuffer.setBufferSize(BUFFER_SIZE);

		if(optind==args.length)
			{
			System.err.println("No bam defined");
			return;
			}
		SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

		while(optind!=args.length)
			{
			BamBuffer buf=new BamBuffer(new File(args[optind++]));
			//buf.setBufferSize(BUFFER_SIZE);
			//buf.setMinQual(this.minQual);
			this.bams.add(buf);
			}
		//this.test();
		this.run();
		}

	public static void main(String[] args)
		{
		try
			{
			new CNV20120927().run(args);
			}
		catch (Exception e)
			{
			e.printStackTrace();
			}
		}
}
