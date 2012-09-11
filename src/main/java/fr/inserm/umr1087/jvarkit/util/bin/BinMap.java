package fr.inserm.umr1087.jvarkit.util.bin;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

public interface BinMap<T>
	{
	public Set<String> getChromosomes();
	public List<T> get(String chromId,int chromStart,int chromEnd,Overlap overlap);
	public int count(String chromId,int chromStart,int chromEnd,Overlap overlap);
	public T put(String chromId,int chromStart,int chromEnd,T object);
	public Iterator<T> iterator(String chromId,int chromStart,int chromEnd,Overlap overlap);
	}
