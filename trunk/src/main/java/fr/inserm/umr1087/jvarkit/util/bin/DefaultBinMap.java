package fr.inserm.umr1087.jvarkit.util.bin;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;

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
	
	public boolean contains(String chromId, int chromStart, int chromEnd,
			Overlap overlap)
		{
		return iterator(chromId,chromStart, chromEnd, overlap).hasNext();
		}
	
	public static DefaultBinMap<String[]> readBedFile(File bedFile) throws IOException
		{
		DefaultBinMap<String[]> map=new DefaultBinMap<String[]>();
		Pattern tab=Pattern.compile("[\t]");
		BufferedReader in=new BufferedReader(new FileReader(bedFile));
		String line;
		while((line=in.readLine())!=null)
			{
			String tokens[]=tab.split(line);
			if(tokens.length<3) throw new IOException("Bad BED row "+line+" in "+bedFile);
			map.put(tokens[0],Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2]), tokens);
			}
		
		in.close();
		return map;
		}
}
