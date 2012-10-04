package fr.inserm.umr1087.jvarkit.segments;

public class TidPosition implements Comparable<TidPosition>
	{
	private int tid;
	private int position;
	
	public TidPosition(int tid,int position)
		{
		this.tid=tid;
		this.position=position;
		}
	
	public int getChromId()
		{
		return tid;
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
		result = prime * result + tid;
		result = prime * result + position;
		return result;
		}


	@Override
	public boolean equals(Object obj)
		{
		if (this == obj) { return true; }
		if (obj == null) { return false; }
		if (getClass() != obj.getClass()) { return false; }
		TidPosition other = (TidPosition) obj;
		if (position != other.position) { return false; }
		if (tid!=other.tid) { return false; }
		return true;
		}


	@Override
	public int compareTo(TidPosition o)
		{
		int i=tid-o.tid;
		if(i!=0) return i;
		return position-o.position;
		}
	
	@Override
	public String toString()
		{
		return String.valueOf(tid)+":"+position;
		}
	
	}
