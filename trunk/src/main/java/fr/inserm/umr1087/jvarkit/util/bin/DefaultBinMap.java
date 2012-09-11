package fr.inserm.umr1087.jvarkit.util.bin;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class DefaultBinMap<T> implements BinMap<T>
	{
	private Map<String, DefaultBinList<T>> chrom2bin=new TreeMap<String, DefaultBinList<T>>();
	@Override
	public Set<String> getChromosomes()
		{
		return Collections.unmodifiableSet(chrom2bin.keySet());
		}

	@Override
	public List<T> get(String chromId, int chromStart, int chromEnd,
			Overlap overlap
			)
		{
		DefaultBinList<T> L=chrom2bin.get(chromId);
		if(L==null) return Collections.emptyList();
		return L.get(chromStart, chromEnd, overlap);
		}

	@Override
	public int count(String chromId, int chromStart, int chromEnd,
			Overlap overlap)
		{
		DefaultBinList<T> L=chrom2bin.get(chromId);
		if(L==null) return 0;
		return L.count(chromStart, chromEnd, overlap);
		}

	@Override
	public T put(String chromId, int chromStart, int chromEnd, T object) {
		DefaultBinList<T> L=chrom2bin.get(chromId);
		if(L==null)
			{
			L=new DefaultBinList<T>();
			chrom2bin.put(chromId, L);
			}
		return L.put(chromStart, chromEnd, object);
		}

	@Override
	public Iterator<T> iterator(String chromId, int chromStart, int chromEnd,
			Overlap overlap)
		{
		DefaultBinList<T> L=chrom2bin.get(chromId);
		if(L==null)
			{
			List<T> L2=Collections.emptyList();
			return L2.iterator();
			}
		return L.iterator(chromStart, chromEnd, overlap);
		}

}
