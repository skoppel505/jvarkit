package fr.inserm.umr1087.jvarkit.tools.cnv20120907;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.interpolation.LoessInterpolator;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import fr.inserm.umr1087.jvarkit.io.TupleBinding;
import fr.inserm.umr1087.jvarkit.util.bin.DefaultBinList;
import fr.inserm.umr1087.jvarkit.util.bin.Overlap;
import fr.inserm.umr1087.jvarkit.util.intervalparser.IntervalParser;

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

public class CNV20120907
	{
	private int BUFFER_SIZE=1000000;//1E6
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
	private static class QualCount
		{
		int qual;
		//int count;
		}
	
	
	private ReferenceSequenceFile reference=null;
	private int windowSize=100;
	private int windowStep=50;
	private List<BamBuffer> bams=new ArrayList<BamBuffer>();
	private List<WigBuffer> bbInput=new ArrayList<WigBuffer>();
	private List<QualCount> qual2count=new ArrayList<QualCount>();
	private Interval targetInterval=null;
	private PrintStream output=System.out;
	
	private class WigBuffer
		{
		BBFileReader bw;
		String chromName=null;
		int start0=0;
		int end0=0;
		DefaultBinList<Float> binIndex=null;
		double getMean(String chromId,int chromStart0,int chromEnd0)
			{
			if(!chromId.equals(this.chromName) ||
				this.binIndex==null ||
				this.start0> chromStart0 ||
				this.end0<=chromEnd0
				)
				{
				if(!this.bw.isBigWigFile()) return Double.NaN;
				this.binIndex=new DefaultBinList<Float>();
				this.chromName=chromId;
				this.start0=Math.max(chromStart0-CNV20120907.this.windowSize,0);
				this.end0=Math.max(chromEnd0,this.start0+CNV20120907.this.windowSize)+BUFFER_SIZE;
				LOG.info("Fill buffer for "+chromName+":"+start0+"-"+end0+" "+bw.getBBFilePath());
				
				
				BigWigIterator iter=this.bw.getBigWigIterator(
						this.chromName,
						this.start0,
						this.chromName,
						this.end0,
						false
						);
				while(iter.hasNext())
					{
					WigItem witem=iter.next();
					this.binIndex.put(witem.getStartBase(), witem.getEndBase(), witem.getWigValue());
					}
				}
			double count=0;
			double total=0;
			for(
					Iterator<Float> i=this.binIndex.iterator(chromStart0, chromEnd0,Overlap.QUERY_OVERLAP_FEATURE);
					i.hasNext();
					)
				{
				total+=i.next();
				count+=1;
				}
			if(count==0) return Double.NaN;
			return total/count;
			}

		}
	
	private class BamBuffer
		{
		File file;
		SAMFileReader samReader;
		String chromName=null;
		int start0=0;
		int end0=0;
		int[][] qual2depth=null;
		double getDepth(int qualityIndex,String chromId,int chromStart0,int chromEnd0)
			{
			if(!chromId.equals(this.chromName) ||
				this.qual2depth==null ||
				this.qual2depth[qualityIndex]==null ||
				this.start0> chromStart0 ||
				this.end0<chromEnd0
				)
				{
				
				this.qual2depth=new int[CNV20120907.this.qual2count.size()][];
                                this.chromName=chromId;
                                this.start0=Math.max(chromStart0-CNV20120907.this.windowSize,0);
                                this.end0=Math.max(chromEnd0,this.start0+CNV20120907.this.windowSize)+BUFFER_SIZE;
                                LOG.info("Fill buffer for "+chromName+":"+start0+"-"+end0+" "+file);
                                Interval interval=new Interval(chromId, this.start0+1, this.end0);//Coordinates are 1-based closed ended.
                                IntervalList iL=new IntervalList(this.samReader.getFileHeader());
                                iL.add(interval);
       
                                SamLocusIterator sli=new SamLocusIterator(this.samReader,iL,true);
                                if(CNV20120907.this.qual2count.size()==1)
                                	{
                                	sli.setMappingQualityScoreCutoff(CNV20120907.this.qual2count.get(0).qual);
                                	}
                                for(int i=0;i< CNV20120907.this.qual2count.size();++i)
                                        {
                                        qual2depth[i]=new int[this.end0-this.start0];
                                        Arrays.fill(this.qual2depth[i], 0);
                                        }
                                for(Iterator<SamLocusIterator.LocusInfo>  iter=sli.iterator();
                                                iter.hasNext();
                                                )
                                        {
                                		
                                        SamLocusIterator.LocusInfo locusInfo=iter.next();
                                        int pos0= locusInfo.getPosition() - 1;//1-based
                                        int offset=pos0-this.start0;
                                        
                                        if(CNV20120907.this.qual2count.size()==1)
                                        	{
                                        	//already filtered on quality, see above
                                        	if(offset>=this.qual2depth[0].length )
	                                             {
	                                             System.err.println("???");
	                                             continue;
	                                             }
                                        	this.qual2depth[0][offset]+=locusInfo.getRecordAndPositions().size();
                                        	}
                                        else
	                                        {
	                                        
	                                        for(RecordAndOffset rao:locusInfo.getRecordAndPositions())
	                                                {
	
	                                                for(int i=0;i< CNV20120907.this.qual2count.size();++i)
	                                                        {
	                                                        QualCount qc= CNV20120907.this.qual2count.get(i);
	                                                    	if(rao.getRecord().getMappingQuality()< qc.qual) continue;
	
	                                                        
	                                                        if(offset>=this.qual2depth[i].length )
	                                                                {
	                                                                System.err.println("???");
	                                                                continue;
	                                                                }
	                                                        this.qual2depth[i][offset]++;
	                                                        }
	                                                }
	                                        }
                                        }
                                //sli.close();
				}
			double count=0;
			double depth=0;
			for(int i=chromStart0;
					i< chromEnd0 && (i-this.start0) < this.qual2depth[qualityIndex].length ;
					++i)
				{
				count+=1;
				depth+= this.qual2depth[qualityIndex][i-this.start0];
				}
			return depth/count;
			}
		}
	
	private class ReferenceBuffer
		{
		String chromName=null;
		int start0=0;
		byte buffer[]=null;

		byte getBaseAt(String chromId,int index0)
			{
			if(!chromId.equals(this.chromName) ||
				buffer==null ||
				index0 < this.start0 ||
				index0 >= (this.start0+buffer.length))
				{
				this.chromName=chromId;
				this.start0 =Math.max(0,index0- CNV20120907.this.windowSize);
				int end=CNV20120907.this.reference.getSequenceDictionary().getSequence(chromId).getSequenceLength();
				end=Math.min(end, Math.max(this.start0+BUFFER_SIZE,this.start0+windowSize+1));
				LOG.info("refill DNA buffer  "+chromId+":"+start0);
				ReferenceSequence dna= CNV20120907.this.reference.getSubsequenceAt(
							chromId,
							this.start0+1,//1-based
							end//inclusive
							);
				this.buffer=dna.getBases();
				}
			return this.buffer[index0-this.start0];
			}
		}

	private ReferenceBuffer referenceBuffer=new ReferenceBuffer();
	
	private static class DataRow
		{
		String chrom;
		int chromStart=0;
		int chromEnd=0;
		double gc_percent=0;
		double wigs[];
		double depths[];
		
		void print(PrintStream out)
			{
			/** chromosome start-end*/
			out.print(this.chrom);
			out.print('\t');
			out.print(chromStart);
			out.print('\t');
			out.print(chromEnd);
			
			/** GC% */
			out.print('\t');
			out.print(gc_percent);
			
			/** loop over the big files */
			for(double wigV: this.wigs)
				{
				out.print('\t');
				out.printf("%.2f",wigV);
				}

			
			/** get coverage */
			for(double d:depths)
				{
				out.print('\t');
				out.printf("%.2f",d);
				}
			out.println();
			}
		}
	
	
	private class DataRowBinding implements TupleBinding<DataRow>
		{
		@Override
		public DataRow readObject(DataInputStream in) throws IOException
			{
			DataRow row=new DataRow();
			row.chrom=in.readUTF();
			row.chromStart=in.readInt();
			row.chromEnd=in.readInt();
			row.gc_percent=in.readDouble();
			int n=in.readInt();
			row.wigs=new double[n];
			for(int i=0;i< n;++i) row.wigs[i]=in.readDouble();
			n=in.readInt();
			row.depths=new double[n];
			for(int i=0;i< n;++i) row.depths[i]=in.readDouble();
			return row;
			}
		@Override
		public void writeObject(DataOutputStream out, DataRow row)
				throws IOException
			{
			out.writeUTF(row.chrom);
			out.writeInt(row.chromStart);
			out.writeInt(row.chromEnd);
			out.writeDouble(row.gc_percent);
			out.writeInt(row.wigs.length);
			for(int i=0;i< row.wigs.length ;++i)
				{
				out.writeDouble(row.wigs[i]);
				}
			out.writeInt(row.depths.length);
			for(int i=0;i< row.depths.length ;++i)
				{
				out.writeDouble(row.depths[i]);
				}
			}
		}
	private DataRowBinding dataRowBinding=new DataRowBinding();
	
	private CNV20120907()
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
		for(WigBuffer wigBuffer: this.bbInput)
			{
			for(int i=chromStart;i<chromEnd;++i)
				{
				System.out.print((i-chromStart)%10==0?'\n':' ');
				System.out.printf(
						"%d:%2.2f",i,wigBuffer.getMean(
						chrom,
						i,
						i+10
						));
				System.out.print(":"+wigBuffer.binIndex.get(i, i+10,Overlap.QUERY_OVERLAP_FEATURE));
				}
			System.out.println();
			
			System.out.println();
			}
			
			
		System.out.println();
		System.exit(0);
		}

	private void writeDataRow(DataRow row)
		{
		row.print(this.output);
		}
	
	
	
	private void redonRoutine() throws IOException,MathException
		{
		File dataFile=null;//TODO
		LoessInterpolator loessInterpolator=new LoessInterpolator();
		SAMSequenceDictionary 	dict=this.reference.getSequenceDictionary();
		String prevChrom=null;
		DataInputStream dataIn=new DataInputStream(new FileInputStream(dataFile));
		List<DataRow> rows=new ArrayList<CNV20120907.DataRow>();
		for(;;)
			{
			DataRow row=null;
			if(dataIn.available()>0)
				{
				row=this.dataRowBinding.readObject(dataIn);
				}
			if((row==null || !row.chrom.equals(prevChrom)) && !rows.isEmpty()) 
				{
				double xval[]=new double[rows.size()];
				double yval[]=new double[rows.size()];
				
				int col=0;
				for(BamBuffer r:bams)
					{
					for(int qualIndex=0;qualIndex< this.qual2count.size();++qualIndex)
						{
						for(int i=0;i< rows.size();++i)
							{
							xval[i]=rows.get(i).depths[col];
							yval[i]=rows.get(i).gc_percent;
							}
						
						col++;
						}
					}
				rows.clear();
				}
			if(row==null) break;
			rows.add(row);
			}
		dataIn.close();
		
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			Count c=chrom2count.get(chrom.getSequenceName());
			if(c.count==0) continue;
			double xval[]=new double[c.count];
			double yval[]=new double[c.count];
			dataIn=new DataInputStream(new FileInputStream(dataFile));
			int rowIndex=0;
			while(rowIndex<c.count)
				{
				DataRow row=this.dataRowBinding.readObject(dataIn);
				if(!row.chrom.equals(chrom.getSequenceName())) continue;
				
				xval[rowIndex]=(row.chromStart+row.chromEnd)/2.0;
				yval[rowIndex]=row.gc_percent;
				
				++rowIndex;
				}
			dataIn.close();
			if(rowIndex!=c.count) throw new RuntimeException();
			
			double y2val[]=loessInterpolator.smooth(xval, yval);
			}
		}
	
	private void run() throws Exception
		{
		this.output.print("#chrom\tchromStart\tchromEnd\tGC%");
		for(WigBuffer wigBuffer: this.bbInput)
			{
			this.output.print("\t"+wigBuffer.bw.getBBFilePath());
			}
		for(BamBuffer b: this.bams)
			{
			for(int qualIndex=0;qualIndex< this.qual2count.size();++qualIndex)
				{
				this.output.print("\t"+b.file+"[Q="+this.qual2count.get(qualIndex).qual+"]");
				}
			}
		this.output.println();
		
		DataRow dataRow=new DataRow();
		dataRow.wigs=new double[this.bbInput.size()];
		dataRow.depths=new double[this.bams.size()*this.qual2count.size()];
		
		SAMSequenceDictionary 	dict=this.reference.getSequenceDictionary();
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			if(targetInterval!=null && !targetInterval.getSequence().equals(chrom.getSequenceName()))
				{
				continue;
				}
			LOG.info("chrom:="+chrom.getSequenceName()+"/"+chrom.getSequenceLength());
			dataRow.chrom=chrom.getSequenceName();
			dataRow.chromStart=(targetInterval!=null?targetInterval.getStart():0);
			while(dataRow.chromStart+this.windowSize<=
					(targetInterval!=null?targetInterval.getEnd():chrom.getSequenceLength()))
				{
				int gc=0;
				int N=0;
				for(int i=dataRow.chromStart;N==0 && i<dataRow.chromStart+this.windowSize;++i)
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
					++dataRow.chromStart;
					continue;
					}
				dataRow.chromEnd=dataRow.chromStart+windowSize;
				dataRow.gc_percent=gc/(double)windowSize;
				
				
				int loopIdx=0;
				/** loop over the big files */
				for(WigBuffer wigBuffer: this.bbInput)
					{
					dataRow.wigs[loopIdx++]= wigBuffer.getMean(
							dataRow.chrom,
							dataRow.chromStart,
							dataRow.chromEnd
							);
						
					}
				
				loopIdx=0;
				/** get coverage */
				for(BamBuffer r:bams)
					{
					for(int qualIndex=0;qualIndex< this.qual2count.size();++qualIndex)
						{
						dataRow.wigs[loopIdx++]=r.getDepth(
								qualIndex,
								dataRow.chrom,
								dataRow.chromStart,
								dataRow.chromEnd
								);
						}
					}				
				writeDataRow(dataRow);
				
				dataRow.chromStart+=this.windowStep;
				}
			}
		}
	
	private void usage()
		{
		System.err.println(
			"Options:\n"+
			" -f (fasta) reference fasta indexed with faidx\n"+
			" -bw (file) add this bigwig file\n"+
			" -q <int> add this quality treshold for BAM alignments\n"+
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
			else if(args[optind].equals("-f") && optind+1 < args.length )
				{
				this.reference= ReferenceSequenceFileFactory.getReferenceSequenceFile(
						new File(args[++optind]),
						true //cut after first whitespace
						);
				}
			else if(args[optind].equals("-bw") && optind+1 < args.length )
				{
				WigBuffer wb=new WigBuffer();
				wb.bw=new BBFileReader(args[++optind]);
				this.bbInput.add(wb);
				}
			else if(args[optind].equals("-q") && optind+1 < args.length )
				{
				QualCount qc=new QualCount();
				qc.qual=Integer.parseInt(args[++optind]);
				this.qual2count.add(qc);
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
		if(this.qual2count.isEmpty())
			{
			QualCount qc=new QualCount();
			qc.qual=0;
			this.qual2count.add(qc);
			}
		if(this.reference==null)
			{
			System.err.println("Reference missing");
			return;
			}
		if(optind==args.length)
			{
			System.err.println("No bam defined");
			return;
			}
		SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

		while(optind!=args.length)
			{
			BamBuffer buf=new BamBuffer();
			buf.file=new File(args[optind++]);
			LOG.info("opening "+buf.file);
			buf.samReader=new SAMFileReader(buf.file);
			buf.samReader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
			this.bams.add(buf);
			}
		//this.test();
		this.run();
		}

	public static void main(String[] args)
		{
		try
			{
			new CNV20120907().run(args);
			}
		catch (Exception e)
			{
			e.printStackTrace();
			}
		}
}
