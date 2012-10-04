package fr.inserm.umr1087.jvarkit.util.picard;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

public class AbstractTView
	{
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");

	private int columns=80;
	private int max_rows=25;
	private SAMFileReader samReader;
	private String chromosome;
	private int position1;
	private ReferenceSequenceFile referenceSequenceFile;
	private int min_dist=2;
	private Map<Integer,Integer> insertionsInRef=new TreeMap<Integer,Integer>();
	
	
	
	
	protected boolean overlap(SAMRecord r1,SAMRecord r2)
		{
		return !(
				r1.getAlignmentEnd()+this.min_dist < r2.getAlignmentStart() ||
				r1.getAlignmentStart() > r2.getAlignmentEnd()+this.min_dist
				);
		}
	
	public void build()
		{
		if(this.referenceSequenceFile!=null)
			{
			SAMSequenceRecord ssr=this.referenceSequenceFile.getSequenceDictionary().getSequence(this.chromosome);
			if(ssr!=null)
				{
				int L=ssr.getSequenceLength();
				ReferenceSequence refseq=this.referenceSequenceFile.getSubsequenceAt(
						this.chromosome,
						this.position1,
						Math.min(position1+this.columns, L)
						);
				byte bases[]=refseq.getBases();
				//todo
				}
			}
		List<SAMRecord> records=new ArrayList<SAMRecord>();

		SAMRecordIterator iter=this.samReader.queryOverlapping(
				this.chromosome,
				this.position1,
				this.position1+this.columns
				);
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			Iterator<CigarAlignment> c=CigarAlignment.iterator(rec);
			while(c.hasNext())
				{
				CigarAlignment caln=c.next();
				
				if(	caln.isInsertRef() && caln.getIndexInCigarElement()==0)
					{
					}
				}
			
			records.add(rec);
			}
		iter.close();
		LOG.info("records:"+records.size()+" "+insertionsInRef);
		/* sort reads */
		Collections.sort(records,new Comparator<SAMRecord>()
				{
				@Override
				public int compare(SAMRecord o1, SAMRecord o2)
					{
					int i=o1.getAlignmentStart()-o2.getAlignmentStart();
					if(i!=0) return i;
					i=o1.getAlignmentEnd()-o2.getAlignmentEnd();
					if(i!=0) return i;
					i=o1.getReadName().compareTo(o2.getReadName());
					return i;
					}
				});
		/* pack reads */
		List<List<SAMRecord>> lines=new ArrayList<List<SAMRecord>>();
		while(!records.isEmpty())
			{
			SAMRecord first=records.remove(0);
			for(int y=0;y< lines.size() && first!=null;++y)
				{
				List<SAMRecord> rowxy=lines.get(y);
				boolean ok=true;
				for(SAMRecord xy:rowxy)
					{
					if(overlap(first,xy))
						{
						ok=false;
						break;
						}
					}
				if(!ok) continue;
				rowxy.add(first);
				first=null;
				}
			if(first==null) continue;
			List<SAMRecord> newrow=new ArrayList<SAMRecord>();
			newrow.add(first);
			lines.add(newrow);
			}
				
		for(int y=0;y< lines.size();++y)
			{
			List<SAMRecord> row=lines.get(y);
			int prev_x=0;
			for(SAMRecord rec:row)
				{
				Iterator<CigarAlignment> c=CigarAlignment.iterator(rec);
				while(c.hasNext())
					{
					CigarAlignment caln=c.next();
					int x=caln.getReferencePosition1()-this.position1;
					if(x<0) continue;
					if(x>=this.columns) continue;
					while(prev_x<x)
						{
						System.out.print(':');
						++prev_x;
						}
					prev_x=x+1;
					if(caln.isInsertRef())
						{
						System.out.print("!"+caln.getReadBase());
						}
					else if(caln.isDeletionRef())
						{
						System.out.print('-');
						}
					else
						{
						System.out.print(caln.getReadBase());
						}
					}
				}
			while(prev_x< this.columns)
				{
				System.out.print('~');
				++prev_x;
				}
			
			
			for(SAMRecord rec:row)
				{
				System.out.print(" "+rec.getCigarString());
				}
			
			System.out.println();
			}
		
		}
	

	public static void main(String[] args)
		{
		args=new String[]{"/commun/data/users/cfaucheron/aln_20120329/S0529/data_S0529/S0529_sort.nodup.bam"};
		
		AbstractTView tview=new AbstractTView();
		tview.referenceSequenceFile=ReferenceSequenceFileFactory.getReferenceSequenceFile(new File("/commun/data/pubdb/ucsc/hg19/chromosomes/hg19.fa"));
		tview.chromosome="chr1";
		tview.columns=50;
		tview.position1=9988;
		tview.samReader=new SAMFileReader(new File(args[0]));
		tview.build();
		LOG.info("Done");
		}
}
