package fr.inserm.umr1087.jvarkit.tools.cnv20120907;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import fr.inserm.umr1087.jvarkit.util.bin.DefaultBinList;
import fr.inserm.umr1087.jvarkit.util.bin.Overlap;

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
	private static final int BUFFER_SIZE=1000000;//1E6
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
	List<QualCount> qual2count=new ArrayList<QualCount>();
	
	
	private class WigBuffer
		{
		BBFileReader bw;
		String chromName=null;
		int start0=0;
		int end0=0;
		DefaultBinList<Float> binIndex=null;
		double getDepth(String chromId,int chromStart0,int chromEnd0)
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
	
				SamLocusIterator sli=new SamLocusIterator(this.samReader,iL);
				
				for(int i=0;i< CNV20120907.this.qual2count.size();++i)
					{
					qual2depth[i]=new int[end0-start0];
					Arrays.fill(this.qual2depth[i], 0);
					}
				for(Iterator<SamLocusIterator.LocusInfo>  iter=sli.iterator();
						iter.hasNext();
						)
					{
					SamLocusIterator.LocusInfo locusInfo=iter.next();
					int pos0= locusInfo.getPosition() - 1;//1-based
					int offset=pos0-this.start0;
					
					for(RecordAndOffset rao:locusInfo.getRecordAndPositions())
						{
						for(int i=0;i< CNV20120907.this.qual2count.size();++i)
							{
							QualCount qc= CNV20120907.this.qual2count.get(i);
							if(rao.getBaseQuality()< qc.qual) continue;
							
							if(offset>=this.qual2depth[i].length ) continue;
							this.qual2depth[i][offset]++;
							}
						}
					}
				sli.close();
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
				
				this.start0 =Math.max(0,index0- CNV20120907.this.windowSize);
				LOG.info("refill DNA buffer  "+chromId+":"+start0);
				ReferenceSequence dna= CNV20120907.this.reference.getSubsequenceAt(
							chromId,
							this.start0+1,//1-based
							Math.max(this.start0+BUFFER_SIZE,this.start0+windowSize+1)//inclusive
							);
				this.buffer=dna.getBases();
				this.chromName=chromId;
				}
			return this.buffer[index0-this.start0];
			}
		}

	private ReferenceBuffer referenceBuffer=new ReferenceBuffer();

	private CNV20120907()
		{
		
		}

	

	
	private void run() throws Exception
		{
		SAMSequenceDictionary 	dict=this.reference.getSequenceDictionary();
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			LOG.info("chrom:="+chrom.getSequenceName()+"/"+chrom.getSequenceLength());
			int start=0;
			while(start+this.windowSize<= chrom.getSequenceLength())
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
				
				/** chromosome start-end*/
				System.out.print(chrom.getSequenceName());
				System.out.print('\t');
				System.out.print(start);
				System.out.print('\t');
				System.out.print(start+windowSize);
				
				/** GC% */
				System.out.print('\t');
				System.out.print(gc/(double)windowSize);
				
				/** loop over the big files */
				for(WigBuffer wigBuffer: this.bbInput)
					{
					
					System.out.print('\t');
					System.out.printf("%.2f",wigBuffer.getDepth(
							chrom.getSequenceName(),
							start,
							start+windowSize
							));
						
					}

				
				/** get coverage */
				for(BamBuffer r:bams)
					{
					for(int qualIndex=0;qualIndex< this.qual2count.size();++qualIndex)
						{
						System.out.print('\t');
						System.out.printf("%.2f",r.getDepth(
								qualIndex,
								chrom.getSequenceName(),
								start,
								start+windowSize
								));
						}
					}
				
				System.out.println();
				
				start+=this.windowStep;
				}
			}
		}
	
	private void usage()
		{
		System.err.println(
			"Options:\n"+
			" -f (fasta) reference fasta indexed with faidx\n"+
			" -bw (file) add this bigwig file\n"+
			" -q <int> add this quality treshold for BAM alignments\n"
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
		while(optind!=args.length)
			{
			BamBuffer buf=new BamBuffer();
			buf.file=new File(args[optind++]);
			LOG.info("opening "+buf.file);
			buf.samReader=new SAMFileReader(buf.file);
			this.bams.add(buf);
			}
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
