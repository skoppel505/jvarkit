package fr.inserm.umr1087.jvarkit.tools.bamstats04;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class BamStats04
	{
	private File bedFile=null;
	private boolean skipDuplicates=true;
	private int minQual=0;
	private int basesperbin=10;
	private int num_bin=20;
	private BamStats04()
		{
		
		
		}
	
	private void scan(File bam) throws Exception
		{
		long bases2count[]=new long[num_bin];
		Arrays.fill(bases2count, 0L);
		Pattern tab=Pattern.compile("[\t]");
		String tokens[];
		BufferedReader in=new BufferedReader(new FileReader(this.bedFile));
		SAMFileReader samReader = new SAMFileReader(bam);
		long total=0L;
		String line=null;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			tokens=tab.split(line,5);
			if(tokens.length<3) throw new IOException("bad bed line in "+line+" "+this.bedFile);
			String chrom=tokens[0];
			int chromStart=Integer.parseInt(tokens[1]);
			int chromEnd=Integer.parseInt(tokens[2]);
			IntervalList L=new IntervalList(samReader.getFileHeader());
			/* picard javadoc:  - Sequence name - Start position (1-based) - End position (1-based, end inclusive)  */
			Interval interval=new Interval(chrom, chromStart+1, chromEnd);
			L.add(interval);
			int counts[]=new int[chromEnd-chromStart];
			Arrays.fill(counts, 0);
			SamLocusIterator iter=new SamLocusIterator(samReader, L, true);
			iter.setMappingQualityScoreCutoff(this.minQual);
			iter.setEmitUncoveredLoci(true);
			while(iter.hasNext())
				{
				SamLocusIterator.LocusInfo li=iter.next();
				if(li.getPosition()< interval.getStart()) continue;
				if(li.getPosition()> interval.getEnd()) continue;
				int count=0;
				for(SamLocusIterator.RecordAndOffset sr: li.getRecordAndPositions())
					{
					if(skipDuplicates && sr.getRecord().getDuplicateReadFlag() ) continue;
					count++;
					}
				counts[li.getPosition()-interval.getStart()]=count;
				}
			iter.close();
			
			for(int depth:counts)
				{
				int cat=depth/this.basesperbin;
				if(cat>=bases2count.length) cat=bases2count.length-1;
				bases2count[cat]++;
				++total;
				}
			
			}
		in.close();
		samReader.close();
		System.out.print(bam.toString()+"\t"+total);
		for(long bases:bases2count)
			{
			System.out.print("\t");
			System.out.print(bases);
			}
		System.out.println();
		}

	
	private void run(String args[]) throws Exception
		{
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.LENIENT);
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				System.out.println("Pierre Lindenbaum PhD. 2013.");
				System.out.println(" -b bedfile (required).");
				System.out.println(" -q (int) min-qual.");
				System.out.println(" -D do NOT ignore duplicates.");
				return;
				}
			else if(args[optind].equals("-b") && optind+1< args.length)
				{
				this.bedFile=new File(args[++optind]);
				}
			else if(args[optind].equals("-q") && optind+1< args.length)
				{
				this.minQual=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-D"))
				{
				this.skipDuplicates=false;
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
		if(optind==args.length)
				{
				System.err.println("Bam input missing");
				System.exit(-1);
				}
		
		System.out.print("#filename\ttotal_bases");
		

		for(int i=0;i< this.num_bin;++i)
			{
			System.out.print("\t[" + (i*this.basesperbin));
			if(i+1==this.num_bin)
				{
				System.out.print("-all[");
				}
			else
				{
				System.out.print("-" + ((i+1)*this.basesperbin) + "[");
				}
			}
		System.out.println();

		while(optind< args.length)
			{
			File bam=new File(args[optind++]);
			scan(bam);
			}
				
		
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception
		{
		BamStats04 app=new BamStats04();
		app.run(args);
		}

}
