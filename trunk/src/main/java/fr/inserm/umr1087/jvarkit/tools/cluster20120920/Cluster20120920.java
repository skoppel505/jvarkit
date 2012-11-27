package fr.inserm.umr1087.jvarkit.tools.cluster20120920;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordCoordinateComparator;
import net.sf.samtools.SAMRecordIterator;



/**
 * Cluster20120920
 * @author lindenb
 *
 */
public class Cluster20120920
	{
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");

	private int margin=10;
	private int minQual=15;
	private List<Interval> intervalList=new ArrayList<Interval>();
	private List<SamRecordAndFile> buffer=new LinkedList<Cluster20120920.SamRecordAndFile>();
	private List<SamRecordSource> samSources=new ArrayList<SamRecordSource>();
	private SAMRecordCoordinateComparator COMPARATOR=new SAMRecordCoordinateComparator();
	
	private class SamRecordSource
		{
		File bamFile;
		SAMRecord current=null;
		SAMFileReader samReader=null;
		SAMRecordIterator iter=null;
		
		SamRecordSource(File bamFile)
			{
			this.bamFile=bamFile;
			this.samReader=new SAMFileReader(bamFile);
			this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
			this.iter=null;
			}
		
		public SamRecordAndFile peek()
			{
			while(current==null && iter!=null )
				{
				if(iter.hasNext())
					{
					SAMRecord sm=iter.next();
					if(sm.getReadUnmappedFlag()) continue;
					if(sm.getMateUnmappedFlag()) continue;
					if(sm.getNotPrimaryAlignmentFlag()) continue;
					if(sm.getReadFailsVendorQualityCheckFlag()) continue;
					if(!sm.getProperPairFlag()) continue;
					if(sm.getMappingQuality()< Cluster20120920.this.minQual) continue;
					current=sm;
					break;
					}
				else
					{
					iter.close();
					iter=null;
					}
				
				}
			if(current==null) return null;
			return new SamRecordAndFile(current,this.bamFile);
			}
		
		public SamRecordAndFile consumme()
			{
			if(current==null) throw new IllegalStateException();
			SAMRecord  rec=current;
			current=null;
			return new SamRecordAndFile(rec,this.bamFile);
			}
		}

	
	private static class SamRecordAndFile
		{
		SAMRecord record;
		File bamFile;
		SamRecordAndFile(SAMRecord record,File bamFile)
			{
			if(record==null) throw new NullPointerException();
			this.record=record;
			this.bamFile=bamFile;
			}
		@Override
		public String toString() {
			return record.getReferenceIndex()+":"+record.getAlignmentStart();
			}
		}
	
	
	private SamRecordAndFile next()
		{
		boolean flags[]=new boolean[this.samSources.size()];
		for(;;)
			{
			if(!buffer.isEmpty())
				{
				return buffer.remove(0);
				}
			Arrays.fill(flags, false);
			for(int i=0;i< this.samSources.size();++i)
				{
				SamRecordAndFile rec= this.samSources.get(i).peek();
				if(rec==null) continue;
				if(buffer.isEmpty())
					{
					buffer.add(rec);
					flags[i]=true;
					}
				/*  If the two records
				 * are equal enough that their ordering in a
				 * sorted SAM file would be arbitrary,
				 * this method returns 0.
				 */
				else if(Cluster20120920.this.COMPARATOR.fileOrderCompare(
						buffer.get(0).record,
						rec.record) >0)
					{
					buffer.clear();
					Arrays.fill(flags, false);
					buffer.add(rec);
					flags[i]=true;
					}
				else if(Cluster20120920.this.COMPARATOR.fileOrderCompare(
						buffer.get(0).record,
						rec.record) <0)
					{
					//ignore record
					}
				else
					{
					buffer.add(rec);
					flags[i]=true;
					}
				}
			if(buffer.isEmpty()) return null;
			for(int i=0;i< flags.length;++i)
				{
				if(!flags[i]) continue;
				this.samSources.get(i).consumme();
				}
			}
		}
	
	
	private Cluster20120920()
		{
		
		}
	
	
	
	boolean first_dump=true;
	private void dump(
			List<SamRecordAndFile> cluster
			)
		{
		if(cluster.size()> 1)
			{
			SAMRecord front=cluster.get(0).record;
			int start1=Math.min(front.getAlignmentStart(),front.getMateAlignmentStart());
			int end1=start1+Math.abs(front.getInferredInsertSize());
			System.out.print(
				front.getReferenceName()+"\t"+
				start1+"\t"+
				end1+"\t"+
				cluster.size()
				);
			for(SamRecordSource srcFile:this.samSources)
				{
				int c=0;
				for(SamRecordAndFile srf:cluster)
					{
					if(srf.bamFile==srcFile.bamFile) 
						{
						++c;
						}
					}
				System.out.print("\t"+c);
				}	
			System.out.println();
			if(first_dump)
				{
				int c=0;
				for(SamRecordAndFile srf:cluster)
					{
					System.err.println(""+(++c)+" "+srf.bamFile+" "+
							front.getReferenceName()+"\t"+
							start1+"\t"+
							end1+"\t"+
							srf.record.getReadName()+"\t"+
							srf.record.getReferenceName()+"\t"+
							Math.min(srf.record.getAlignmentStart(),srf.record.getMateAlignmentStart())							
							);
					}
				first_dump=false;
				}
			}
		
		cluster.clear();
		}
	
	/** scan */
	private void scanAll()
		{
		List<SamRecordAndFile> cluster=new ArrayList<SamRecordAndFile>();
		SamRecordAndFile rec=null;
		while((rec=next())!=null)
			{
			if(!rec.record.getFirstOfPairFlag()) continue;
			LOG.info("ok");
			if(!cluster.isEmpty())
				{
				SAMRecord front=cluster.get(0).record;
				if(front.getAlignmentStart()>front.getAlignmentEnd()) throw new IllegalStateException();
				int start1=Math.min(front.getAlignmentStart(),front.getMateAlignmentStart())-this.margin;
				int end1=start1+Math.abs(front.getInferredInsertSize())+2*this.margin;
				
				int start2=Math.min(rec.record.getAlignmentStart(),rec.record.getMateAlignmentStart())-this.margin;
				int end2=start2+Math.abs(rec.record.getInferredInsertSize())+2*this.margin;
				if(!front.getReferenceName().equals(rec.record.getReferenceName()) || end2<start1 || start2>end1)
					{
					dump(cluster);
					}
				}
			cluster.add(rec);
			}
		dump(cluster);
		}
	
	private void scan()
		{
		if(this.intervalList.isEmpty())
			{
			for(SamRecordSource src:this.samSources)
				{
				src.iter=src.samReader.iterator();
				}
			scanAll();
			}
		else
			{
			for(Interval interval:this.intervalList)
				{
				for(SamRecordSource src:this.samSources)
					{
					src.iter=src.samReader.queryOverlapping(
							interval.getSequence(),
							interval.getStart()+1,
							interval.getEnd()
							);
					}
				scanAll();
				}
			}
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
			this.samSources.add(new SamRecordSource(new File(args[optind++])));
			}
		System.out.print(
				"#chrom\t"+
				"start\t"+
				"end\t"+
				"total"
				);
			for(SamRecordSource srcFile:this.samSources)
				{
				System.out.print("\t"+srcFile.bamFile);
				}	
			System.out.println();

		this.scan();
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
