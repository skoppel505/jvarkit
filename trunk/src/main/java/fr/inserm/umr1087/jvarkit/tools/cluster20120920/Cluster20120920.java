package fr.inserm.umr1087.jvarkit.tools.cluster20120920;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;




public class Cluster20120920
	{
	@SuppressWarnings("unused")
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");

	private int margin=10;
	private int minQual=15;
	private List<Interval> intervalList=new ArrayList<Interval>();
	private Cluster20120920()
		{
		
		}
	private void dump(
			File bamFile,
			List<SAMRecord> cluster
			)
		{
		if(cluster.size()> 1)
			{
			SAMRecord front=cluster.get(0);
			int start1=Math.min(front.getAlignmentStart(),front.getMateAlignmentStart());
			int end1=start1+Math.abs(front.getInferredInsertSize());
			System.out.println(
				front.getReferenceName()+"\t"+
				start1+"\t"+
				end1+"\t"+
				cluster.size()+"\t"+
				bamFile
				);
			}
		
		cluster.clear();
		}
	private void scan(File bamFile,SAMFileReader samReader,SAMRecordIterator iter)
		{
		List<SAMRecord> cluster=new ArrayList<SAMRecord>();
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.getMateUnmappedFlag()) continue;
			if(rec.getNotPrimaryAlignmentFlag()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(rec.getProperPairFlag()) continue;
			if(rec.getMappingQuality()< this.minQual) continue;
			if(!rec.getFirstOfPairFlag()) continue;
			LOG.info("ok");
			if(!cluster.isEmpty())
				{
				SAMRecord front=cluster.get(0);
				if(front.getAlignmentStart()>front.getAlignmentEnd()) throw new IllegalStateException();
				int start1=Math.min(front.getAlignmentStart(),front.getMateAlignmentStart())-this.margin;
				int end1=start1+Math.abs(front.getInferredInsertSize())+2*this.margin;
				
				int start2=Math.min(rec.getAlignmentStart(),rec.getMateAlignmentStart())-this.margin;
				int end2=start2+Math.abs(rec.getInferredInsertSize())+2*this.margin;
				if(!front.getReferenceName().equals(rec.getReferenceName()) || end2<start1 || start2>end1)
					{
					dump(bamFile,cluster);
					}
				}
			cluster.add(rec);
			}
		dump(bamFile,cluster);
		}
	
	private void scanBam(File bamFile)
		{
		SAMFileReader samReader=new SAMFileReader(bamFile);
		samReader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
		if(this.intervalList.isEmpty())
			{
			SAMRecordIterator iter=samReader.iterator();
			scan(bamFile,samReader,iter);
			iter.close();
			}
		else
			{
			for(Interval interval:this.intervalList)
				{
				SAMRecordIterator iter=samReader.queryOverlapping(
						interval.getSequence(),
						interval.getStart()+1,
						interval.getEnd()
						);
				scan(bamFile,samReader,iter);
				iter.close();
				}
			}
		samReader.close();
		}
	
	private void usage()
		{
		System.err.println(
			"  [options] bam1 bam2 ... bamN"+
			"Options:\n"+
			" -b bed file (optional)\n"+
			" -m <int> margin-tolerance-size["+this.margin+"]\n"+
			" -q <int> add this quality treshold for BAM alignments\n"
			);
		}
	
	private void loadBed(File file) throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		BufferedReader in=new BufferedReader(new FileReader(file));
		String line;
		while((line=in.readLine())!=null)
			{
			String tokens[]=tab.split(line);
			this.intervalList.add(new Interval(
				tokens[0],
				Integer.parseInt(tokens[1]),
				Integer.parseInt(tokens[2])
				));
			}
		in.close();
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
			else if(args[optind].equals("-q") && optind+1 < args.length )
				{
				this.minQual=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-m") && optind+1 < args.length )
				{
				this.margin=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-b") && optind+1 < args.length )
				{
				this.loadBed(new File(args[++optind]));
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
			System.err.println("No bam defined");
			return;
			}
		SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
	
		while(optind!=args.length)
			{
			File file=new File(args[optind++]);
			scanBam(file);
			
			}
		//this.test();
		}
	public static void main(String[] args)
		{
		try
			{
			LOG.setLevel(Level.OFF);
			new Cluster20120920().run(args);
			}
		catch (Exception e)
			{
			e.printStackTrace();
			}
		}
}
