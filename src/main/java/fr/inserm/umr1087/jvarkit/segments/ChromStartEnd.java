package fr.inserm.umr1087.jvarkit.segments;

import fr.inserm.umr1087.jvarkit.util.bin.Bin;

public class ChromStartEnd implements Comparable<ChromStartEnd>
	{
	private String chrom;
	private int start;
	private int end;
	
	public ChromStartEnd(String chrom,int start,int end)
		{
		this.chrom=chrom;
		this.start=start;
		this.end=end;
		}
	
	public String getChromosome()
		{
		return chrom;
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
		result = prime * result + chrom.hashCode();
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
		ChromStartEnd other = (ChromStartEnd) obj;
		if (start != other.start) { return false; }
		if (end != other.end) { return false; }
		if (!chrom.equals(other.chrom)) { return false; }
		return true;
		}


	@Override
	public int compareTo(ChromStartEnd o)
		{
		int i=this.chrom.compareTo(o.chrom);
		if(i!=0) return i;
		i=start-o.start;
		if(i!=0) return i;
		return end-o.end;
		}
	
	@Override
	public String toString()
		{
		return chrom+":"+start+"-"+end;
		}
	
	}
