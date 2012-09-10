package fr.inserm.umr1087.jvarkit.tools.cnv20120907;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

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
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
	private static class QualCount
		{
		int qual;
		int count;
		}
	
	
	private ReferenceSequenceFile reference=null;
	private int windowSize=100;
	private int windowStep=50;
	private int minBaseQual=30;
	private List<SAMFileReader> bams=new ArrayList<SAMFileReader>();
	private List<BBFileReader> bbInput=new ArrayList<BBFileReader>();
	List<QualCount> qual2count=new ArrayList<QualCount>();
	
	private class ReferenceBuffer
		{
		String name=null;
		int start0=0;
		byte buffer[]=null;

		byte getBaseAt(String chromId,int index0)
			{
			if(!chromId.equals(this.name) || buffer==null || index0 < this.start0 || index0 >= (this.start0+buffer.length))
				{
				this.start0 =Math.max(0,this.start0- CNV20120907.this.windowSize);
				ReferenceSequence dna= CNV20120907.this.reference.getSubsequenceAt(
							chromId,
							this.start0+1,//1-based
							(this.start0+windowSize)//inclusive
							);
				this.buffer=dna.getBases();
				this.name=chromId;
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
				Interval interval=new Interval(
						chrom.getSequenceName(),
						start+1,
						(start+windowSize)
						);
				ReferenceSequence dna=this.reference.getSubsequenceAt(
						chrom.getSequenceName(),
						start+1,//1-based
						(start+windowSize)//inclusive
						);
				int gc=0;
				int N=0;
				for(byte base:dna.getBases())
					{
					switch(base)
						{
						case 'n': case 'N': ++N;break;
						case 's': case 'S':
						case 'c': case 'C':
						case 'g': case 'G': ++gc; break;
						default:break;
						}
					if(N>0) break;
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
				for(BBFileReader bbReader: this.bbInput)
					{
					if(bbReader.isBigWigFile())
						{
						double wig_total=0;
						int wig_count=0;
						BigWigIterator iter=bbReader.getBigWigIterator(
								chrom.getSequenceName(),
								start,
								chrom.getSequenceName(),
								start+windowSize,
								false
								);
						while(iter.hasNext())
							{
							WigItem witem=iter.next();
							wig_total+=witem.getWigValue();
							++wig_count;
							}
						System.out.print('\t');
						System.out.printf("%.2f",wig_total/(double)wig_count);
						}
					}

				
				/** get coverage */
				for(SAMFileReader r:bams)
					{
					IntervalList iL=new IntervalList(r.getFileHeader());
					iL.add(interval);//Coordinates are 1-based closed ended. 

					SamLocusIterator sli=new SamLocusIterator(r,iL);
					
					for(QualCount qc:qual2count)
						{
						qc.count=0;
						}
					for(Iterator<SamLocusIterator.LocusInfo>  i=sli.iterator();
							i.hasNext();
							)
						{
						
						SamLocusIterator.LocusInfo locusInfo=i.next();


						for(RecordAndOffset rao:locusInfo.getRecordAndPositions())
							{
							for(QualCount qc:qual2count)
								{
								if(rao.getBaseQuality()< qc.qual) continue;
								qc.count++;
								}
							}
						}
					sli.close();
					for(QualCount qc:qual2count)
						{
						System.out.print('\t');
						System.out.printf("%.2f",qc.count/(double)windowSize);
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
				BBFileReader bbReader=new BBFileReader(args[++optind]);
				this.bbInput.add(bbReader);
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
			File f=new File(args[optind++]);
			LOG.info("opening "+f);
			SAMFileReader sfr=new SAMFileReader(f);
			this.bams.add(sfr);
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
