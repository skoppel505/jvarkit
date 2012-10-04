package fr.inserm.umr1087.jvarkit.segments;

import fr.inserm.umr1087.jvarkit.util.bin.Bin;

public class TidStartEnd implements Comparable<TidStartEnd>
	{
	private int tid;
	private int start;
	private int end;
	
	public TidStartEnd(int tid,int start,int end)
		{
		this.tid=tid;
		this.start=start;
		this.end=end;
		}
	
	public int getChromId()
		{
		return tid;
		}

	public int getStart()
		{
		return start;
		}
	
	public int getEnd()
		{
		return end;
		}

	public int getBin()
		{
		return Bin.bin(start, end);
		}
	
	@Override
	public int hashCode()
		{
		final int prime = 31;
		int result = 1;
		result = prime * result + tid;
		result = prime * result + start;
		result = prime * result + end;
		return result;
		}


	@Override
	public boolean equals(Object obj)
		{
		if (this == obj) { return true; }
		if (obj == null) { return false; }
		if (getClass() != obj.getClass()) { return false; }
		TidStartEnd other = (TidStartEnd) obj;
		if (start != other.start) { return false; }
		if (end != other.end) { return false; }
		if (tid != other.tid) { return false; }
		return true;
		}


	@Override
	public int compareTo(TidStartEnd o)
		{
		int i=tid-o.tid;
		if(i!=0) return i;
		i=start-o.start;
		if(i!=0) return i;
		return end-o.end;
		}
	
	@Override
	public String toString()
		{
		return String.valueOf(tid)+":"+start+"-"+end;
		}
	
	}
