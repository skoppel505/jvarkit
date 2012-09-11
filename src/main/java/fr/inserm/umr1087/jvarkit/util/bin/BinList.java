package fr.inserm.umr1087.jvarkit.util.bin;

import java.util.Iterator;
import java.util.List;

public interface BinList<T>
	{
	public List<T> get(int chromStart,int chromEnd,Overlap overlap);
	public int count(int chromStart,int chromEnd,Overlap overlap);
	public T put(int chromStart,int chromEnd,T object);
	public Iterator<T> iterator(int chromStart,int chromEnd,Overlap overlap);
	}
