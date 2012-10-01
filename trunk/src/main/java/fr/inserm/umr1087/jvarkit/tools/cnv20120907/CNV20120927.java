package fr.inserm.umr1087.jvarkit.tools.cnv20120907;

import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import javax.imageio.ImageIO;

import org.apache.commons.math.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;

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
import fr.inserm.umr1087.jvarkit.tools.genomegraph.GenomeGraph;
import fr.inserm.umr1087.jvarkit.tools.genomegraph.Graphics2DGenomeGraphDrawer;
import fr.inserm.umr1087.jvarkit.util.bin.DefaultBinList;
import fr.inserm.umr1087.jvarkit.util.bin.Overlap;
import fr.inserm.umr1087.jvarkit.util.intervalparser.IntervalParser;
import fr.inserm.umr1087.jvarkit.util.picard.DepthBuffer;
import fr.inserm.umr1087.jvarkit.util.picard.ReferenceBuffer;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class CNV20120927
	{
	private int BUFFER_SIZE=1000000;//1E6
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
	private int minQual;
	private TidStartEndBinding tidStartEndBinding=new TidStartEndBinding();	
	
	
	private ReferenceSequenceFile reference=null;
	private int windowSize=100;
	private int windowStep=50;
	private List<BamBuffer> bams=new ArrayList<BamBuffer>();
	private Interval targetInterval=null;
	private File dbHome=null;
	private Environment environment=null;
	private Database region2depths;
	private Transaction txn=null;
	
	private class BamBuffer extends DepthBuffer
		{
		File file;
		BamBuffer(File file)
			{
			super(new SAMFileReader(file));
			getSamReader().setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
			this.file=file;
			}
		}
	

	private ReferenceBuffer referenceBuffer=null;

	private CNV20120927()
		{
		
		}

	@SuppressWarnings("unused")
	private void test() throws Exception
		{
		String chrom="chr1";
		int chromStart=12000;
		int chromEnd=chromStart+100;
		for(int i=chromStart;i<chromEnd;++i)
			{
			if((i-chromStart)%60==0) System.out.println();
			 System.out.print((char)referenceBuffer.getBaseAt(chrom, i));
			}
			
		System.out.println();
		System.exit(0);
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
	
	private void run() throws Exception
		{
		open();
		LoessInterpolator loessInterpolator=new LoessInterpolator();
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		SAMSequenceDictionary 	dict=this.reference.getSequenceDictionary();
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			if(targetInterval!=null && !targetInterval.getSequence().equals(chrom.getSequenceName()))
				{
				continue;
				}
			int  start=(targetInterval!=null?targetInterval.getStart():0);
			while(start+this.windowSize<=
					(targetInterval!=null?targetInterval.getEnd():chrom.getSequenceLength()))
				{
				int gc=0;
				int N=0;
				for(int i=start;N==0 && i<start+this.windowSize;++i)
					{
					byte base=referenceBuffer.getBaseAt(chrom.getSequenceName(), i);
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
					++start;
					continue;
					}
				TidStartEnd tidStartEnd=new TidStartEnd(
						chrom.getSequenceIndex(),
						start,
						start+windowSize
						);
				CNVRow row=new CNVRow(this.bams.size());
				row.setGcPercent(gc/(double)windowSize);
				
				for(int bx=0;
						bx< this.bams.size();
						++bx)
					{
					row.setDepth(bx,this.bams.get(bx).getMean(
						chrom.getSequenceName(),
						start,
						start+windowSize
						));
						
					}
				
				
				tidStartEndBinding.objectToEntry(tidStartEnd, key);
				CNVRow.BINDING.objectToEntry(row, data);
				
				if(this.region2depths.put(txn, key, data)!=OperationStatus.SUCCESS)
					{
					close();
					throw new RuntimeException("Cannot insert");
					}

				start+=this.windowStep;
				}
			}
		
		
		Cursor c=this.region2depths.openCursor(txn, null);
		List<Point2D.Double> points=new ArrayList<Point2D.Double>();
		OperationStatus status;
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			
			for(int bx=0;
					bx< this.bams.size();
					++bx)
				{
			
				PolynomialSplineFunction splineFun=null; 
				double minDepth= Double.MAX_VALUE;
				List<Double> vMedian=null;
				double medianValue=Double.NaN;
				
				for(int step=0;step<4;++step)
					{
					key=new DatabaseEntry();
					data=new DatabaseEntry();

					boolean first=true;
					points.clear();
					
					switch(step)
						{
						case 2:
							{
							vMedian=new ArrayList<Double>();
							break;
							}
						}
					
					
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
						
						
						CNVRow row=CNVRow.BINDING.entryToObject(data);
						switch(step)
							{
							case 0:
								{
								points.add(new Point2D.Double(
									row.getDepth(bx),
									row.getGcPercent()
									));
								break;
								}
							case 1:
							case 2:
							case 3:
								{
								double newDepth=0;
								switch(step)
									{
									case 1:
										{
										newDepth=splineFun.value(row.getDepth(bx));
										if(minDepth> newDepth) minDepth=newDepth;
										break;
										}
									case 2:
										{	
										newDepth= row.getDepth(bx)- minDepth;
										vMedian.add(newDepth);
										break;
										}
									case 3:
										{
										newDepth= row.getDepth(bx) / medianValue;
										break;
										}
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
							}
						}
					
					
					if(step==0 && points.isEmpty()) break;
					
					switch(step)
						{
						case 0:
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
							double x_val[]=new double[points.size()];
							double y_val[]=new double[points.size()];
							
							for(int i=0;i< points.size();++i)
								{
								x_val[i]=points.get(i).x;
								y_val[i]=points.get(i).y;
								
								if(i>0 && x_val[i]<=x_val[i-1])//loess plante si xi==x(i-1)
									{
									x_val[i]=x_val[i-1]+1E-6;
									}
								}
							points.clear();
							splineFun= loessInterpolator.interpolate(x_val,y_val);
							break;
							}
						case 1:
							{
							splineFun=null;
							break;
							}
						case 2:
							{
							Collections.sort(vMedian);
							medianValue = vMedian.get(vMedian.size()/2);
							vMedian=null;
							break;
							}
						}
						
					}
				}
			}
		
		for(int bx=0;
			bx< this.bams.size();
			++bx)
			{
			GenomeGraph grap=new GenomeGraph(this.reference.getSequenceDictionary());
			File imgFile=new File("/commun/data/users/lindenb/_ignore.backup/"+this.bams.get(bx).file.getName()+".jpg");
			key=new DatabaseEntry();
			status=c.getFirst(key, data, null);
			while(status==OperationStatus.SUCCESS)
				{
				
				TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
				int tid=seg.getChromId();
				int chromStart=seg.getStart();
				int chromEnd=seg.getEnd();
				List<Double> runmed=new ArrayList<Double>();
				CNVRow row=CNVRow.BINDING.entryToObject(data);
				runmed.add(row.getDepth(bx));
				
				
			
				
				/** 'runmed' forward/reverse */
				for(int side=0;side<2;++side)
					{
					Cursor cloneC=c.dup(true);
					DatabaseEntry key2=new DatabaseEntry();
					int runmedix=0;
					while(runmedix<2)
						{
						OperationStatus status2=(side==0?
								cloneC.getPrev(key2, data, null):
								cloneC.getNext(key2, data, null)
								);
						if(status2!=OperationStatus.SUCCESS) break;
						seg=this.tidStartEndBinding.entryToObject(key2);
						if(seg.getChromId()!=tid) break;
						row=CNVRow.BINDING.entryToObject(data);
						runmed.add(row.getDepth(bx));
						runmedix++;
						}
					cloneC.close();
					}
				Collections.sort(runmed);
				grap.put(tid,chromStart,chromEnd,Math.log(runmed.get(runmed.size()/2)));
				status=c.getNext(key, data, null);
				}
			if(grap.isEmpty()) continue;
			Graphics2DGenomeGraphDrawer drawer=new Graphics2DGenomeGraphDrawer();
			drawer.setMinY(-3.0);
			drawer.setMaxY(3.0);
			BufferedImage img=drawer.createImage(grap);
			ImageIO.write(img, "JPG", imgFile);
			}
		c.close();
		
		
		
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
		this.referenceBuffer=new ReferenceBuffer(this.reference);
		this.referenceBuffer.setBufferSize(BUFFER_SIZE);

		if(optind==args.length)
			{
			System.err.println("No bam defined");
			return;
			}
		SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

		while(optind!=args.length)
			{
			BamBuffer buf=new BamBuffer(new File(args[optind++]));
			LOG.info("opening "+buf.file);
			buf.setBufferSize(BUFFER_SIZE);
			buf.setMinQual(this.minQual);
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
