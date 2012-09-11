 package fr.inserm.umr1087.jvarkit.util.bin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;


public class DefaultBinList<T>
	implements BinList<T>
	{
	private static class Carrier<T>
		{
		int chromStart;
		int chromEnd;
		T object;
		
		Carrier(int chromStart,int chromEnd,T object)
			{
			this.chromEnd=chromEnd;
			this.chromStart=chromStart;
			this.object=object;
			}
		
		int getBin()
			{
			return Bin.bin(chromStart, chromEnd);
			}

		}
	private SortedMap<Integer, List<Carrier<T> > > bin2list;
	
	public DefaultBinList()
		{
		this.bin2list=new TreeMap<Integer, List<Carrier<T>>>();
		}
		
	@Override
	public int count(int chromStart, int chromEnd,Overlap overlap)
		{
		int c=0;
		Iterator<T> i=iterator(chromStart,chromEnd, overlap);
		while(i.hasNext())
			{
			i.next();
			c++;
			}
		return c;
		}
	
	@Override
	public List<T> get(int chromStart, int chromEnd,Overlap overlap)
		{
		List<T> L=new ArrayList<T>();
		Iterator<T> i=iterator(chromStart,chromEnd,overlap);
		while(i.hasNext())
			{
			L.add(i.next());
			}
		return L;
		}
	
	@Override
	public Iterator<T> iterator(int chromStart, int chromEnd,Overlap overlap) {
		return new MyIterator(chromStart,chromEnd,overlap);
		}
	
	public T put(int chromStart, int chromEnd, T object)
		{
		Carrier<T> c=new Carrier<T>(chromStart,chromEnd,object);
		int bin=c.getBin();
		List<Carrier<T>> L= this.bin2list.get(bin);
		if(L==null)
			{
			L=new ArrayList<Carrier<T>>();
			this.bin2list.put(bin,L);
			}
		L.add(c);
		return object;
		}
	
	private class MyIterator implements Iterator<T>
		{
		private Overlap overlap;
		private int chromStart;
		private int chromEnd;
		private List<Integer> bins;
		private int binIndex=0;
		private Iterator<Carrier<T>> listiter=null;
		private T _next;
		private boolean _hasNextCalled=false;
		private boolean _hasNext=false;

		MyIterator(int chromStart,int chromEnd,Overlap overlap)
			{
			this.overlap=(overlap==null?Overlap.QUERY_OVERLAP_FEATURE:overlap);
			this.chromStart=chromStart;
			this.chromEnd=chromEnd;
			this.bins=Bin.bins(chromStart, chromEnd);
			}
		@Override
		public boolean hasNext()
			{
			if(_hasNextCalled) return _hasNext;
			_hasNextCalled=true;
			_hasNext=false;
			_next=null;
			
			for(;;)
				{
				/* never called */
				if(listiter==null)
					{
					binIndex=0;
					List<Carrier<T>>  L= DefaultBinList.this.bin2list.get(this.bins.get(binIndex));
					if(L==null) L=Collections.emptyList(); 
					listiter=L.iterator();
					continue;
					}
				/* iterator is eof */
				if(!listiter.hasNext())
					{
					++binIndex;
					/* no more bin */
					if(binIndex>=this.bins.size())
						{
						return false;
						}
					
					List<Carrier<T>>  L= DefaultBinList.this.bin2list.get(this.bins.get(binIndex));
					if(L==null) L=Collections.emptyList(); 
					listiter=L.iterator();
					continue;
					}
				
				Carrier<T> c=listiter.next();
				if(!this.overlap.match(
						this.chromStart, this.chromEnd,
						c.chromStart,c.chromEnd
						))
					{
					continue;
					}
				_next=c.object;
				_hasNext=true;
				return true;
				}
			}
		@Override
		public T next()
			{
			if(!_hasNextCalled) hasNext();
			if(!hasNext()) throw new IllegalStateException();
			T o=_next;
			_hasNext=false;
			_hasNextCalled=false;
			_next=null;
			return o;
			}
		@Override
		public void remove()
			{
			throw new UnsupportedOperationException();
			}
		}
	
	}
