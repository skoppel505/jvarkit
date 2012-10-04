package fr.inserm.umr1087.jvarkit.segments;

public class ChromPosition implements Comparable<ChromPosition>
	{
	private String chrom;
	private int position;
	
	public ChromPosition(String chrom,int position)
		{
		this.chrom=chrom;
		this.position=position;
		}
	
	public String getChromosome()
		{
		return chrom;
		}


	public int getPosition()
		{
		return position;
		}


	@Override
	public int hashCode()
		{
		final int prime = 31;
		int result = 1;
		result = prime * result + chrom.hashCode();
		result = prime * result + position;
		return result;
		}


	@Override
	public boolean equals(Object obj)
		{
		if (this == obj) { return true; }
		if (obj == null) { return false; }
		if (getClass() != obj.getClass()) { return false; }
		ChromPosition other = (ChromPosition) obj;
		if (position != other.position) { return false; }
		if (!chrom.equals(other.chrom)) { return false; }
		return true;
		}


	@Override
	public int compareTo(ChromPosition o)
		{
		int i=this.chrom.compareTo(o.chrom);
		if(i!=0) return i;
		return position-o.position;
		}
	
	@Override
	public String toString()
		{
		return chrom+":"+position;
		}
	
	}
