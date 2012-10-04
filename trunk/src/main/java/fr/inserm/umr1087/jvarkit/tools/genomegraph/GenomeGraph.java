package fr.inserm.umr1087.jvarkit.tools.genomegraph;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class GenomeGraph
	{
	@SuppressWarnings("unused")
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");

	public static class Data
		{
		int start;
		int end;
		double value;
		}
	
	static class ChromInfo
		{
		String name;
		long length;
		Integer minPos=null;
		Integer maxPos=null;
		List<Data> data=new ArrayList<GenomeGraph.Data>();
		@Override
		public String toString() {
			return name+" L="+length+" ("+minPos+"/"+this.maxPos+") data-size="+data.size();
			}
		}
	Double minValue=null;
	Double maxValue=null;
	private SAMSequenceDictionary samSequenceDictionary;
	private List<ChromInfo> index2chrominfo=new ArrayList<ChromInfo>(30);
	
	public GenomeGraph(SAMSequenceDictionary samSequenceDictionary)
		{
		this.samSequenceDictionary=samSequenceDictionary;
		}
	
	public boolean isEmpty()
		{
		return minValue==null || maxValue==null;
		}
	
	List<ChromInfo> getIndex2chrominfo()
		{
		return index2chrominfo;
		}
	
	public SAMSequenceDictionary getSamSequenceDictionary()
		{
		return samSequenceDictionary;
		}
	
	public Data put(String chromName,int start,int end,double value)
		{
		SAMSequenceRecord samSeqRec=getSamSequenceDictionary().getSequence(chromName);
		if(samSeqRec==null)
			{
			System.err.println("Unknown Chromosome "+chromName);
			return null;
			}
		int seqIdx=samSeqRec.getSequenceIndex();
		return put(seqIdx,start,end,value);
		}
	
	public Data put(int chromIndex,int start,int end,double value)
		{
		if(Double.isNaN(value)) return null;
		while(index2chrominfo.size()<=chromIndex)
			{
			index2chrominfo.add(null);
			}
		ChromInfo ci=index2chrominfo.get(chromIndex);
		if(ci==null)
			{
			if(chromIndex >= getSamSequenceDictionary().size())
				{
				throw new IndexOutOfBoundsException();
				}
			SAMSequenceRecord samSeqRec=getSamSequenceDictionary().getSequence(chromIndex);
			ci=new ChromInfo();
			ci.name=samSeqRec.getSequenceName();
			ci.length=samSeqRec.getSequenceLength();
			this.index2chrominfo.set(chromIndex,ci);
			}
		if(minValue==null || minValue>value) this.minValue=value;
		if(maxValue==null || maxValue<value) this.maxValue=value;
		if(ci.minPos==null || ci.minPos>start) ci.minPos=start;
		if(ci.maxPos==null || ci.maxPos<end) ci.maxPos=end;

		Data d=new Data();
		d.start=start;
		d.end=end;
		d.value=value;
		ci.data.add(d);
		return d;
		}

	
	protected void compile()
		{
		
		for(ChromInfo ci:this.index2chrominfo)
			{
			if(ci==null) continue;
			}
		}
	
	}
