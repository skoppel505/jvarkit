package fr.inserm.umr1087.jvarkit.collections;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

public class Algorithms<T>
	{
	private Comparator<T> comparator=null;
	private static final int DEFAULT_MERGE_SORT_SIZE=100000;
	private int mergeSortSize=DEFAULT_MERGE_SORT_SIZE;
	
	public Algorithms(Comparator<T> comparator)
		{
		this.comparator=comparator;
		}
	
	public Algorithms()
		{
		}
	
	public void setComparator(Comparator<T> comparator)
		{
		this.comparator = comparator;
		}
	
	
	@SuppressWarnings("unchecked")
	public int compare(T o1, T o2)
		{
		return comparator==null? ((Comparable<T>)o1).compareTo(o2):
				comparator.compare(o1, o2);
		}
	
	private void swap(List<T> L,int i,int j)
		{
		T o1=L.get(i);
		T o2=L.get(j);
		L.set(i,o2);
		L.set(j,o1);
		}
	
	public void shuffle(List<T> L)
		{
		Random rand=new Random();
		for(int i=0;i<L.size();++i)
			{
			int j=rand.nextInt(L.size());
			if(i!=j) swap(L,i,j);
			}
		}
	
	public void quickSort(List<T> L)
		{
		if(L.size()<2) return;
		quickSort(L,0,L.size()-1);
		if(!isSorted(L)) throw new RuntimeException();
		}
	
	//http://www.vogella.com/articles/JavaAlgorithmsQuicksort/article.html
	private void quickSort(List<T> L,int low, int high)
  		{
	    int i = low, j = high;
	    // Get the pivot element from the middle of the list
	    T pivot = L.get(low + (high-low)/2);

	    // Divide into two lists
	    while (i <= j) {
	      // If the current value from the left list is smaller then the pivot
	      // element then get the next element from the left list
	      while (compare(L.get(i) , pivot)< 0)
	      	{
	        i++;
	      	}
	      // If the current value from the right list is larger then the pivot
	      // element then get the next element from the right list
	      while (compare(L.get(j), pivot)>0) {
	        j--;
	      }

	      // If we have found a values in the left list which is larger then
	      // the pivot element and if we have found a value in the right list
	      // which is smaller then the pivot element then we exchange the
	      // values.
	      // As we are done we can increase i and j
	      if (i <= j) {
	        swap(L,i, j);
	        i++;
	        j--;
	      }
	    }
	    // Recursion
	    if (low < j)
	      quickSort(L,low, j);
	    if (i < high)
	    	quickSort(L,i, high);
	  }
	
	
	public boolean isSorted(List<T> L)
		{
		for(int i=0;i+1< L.size();++i)
			{
			if(compare(L.get(i), L.get(i+1))>0)
				{
				return false;
				}
			}
		return true;
		}
	
	public void mergeSort(List<T> L)
		{
		int mergeLength=mergeSortSize;
		int start=0;
		while(start<L.size())
			{	
			int end=Math.min(L.size(),start+mergeLength);
			quickSort(L,start,end-1);
			int index0=0;
			int index1=start;
			T o0=null;
			T o1=null;
			int index=0;

			while(index0< start || index1 < end)
				{
				if(o0==null && index0<start)
					{
					o0=L.get(index0);
					System.err.println("Pop 0:"+index0+"="+o0);
					index0++;
					}
				if(o1==null && index1 < end)
					{
					o1=L.get(index1);
					System.err.println("Pop 1:"+index0+"="+o1);
					index1++;
					}
				
				if(o0!=null && o1!=null)
					{
					int i=compare(o0,o1);
					if(i==0)
						{
						L.set(index++, o0);
						L.set(index++, o1);
						o0=null;
						o1=null;
						}
					else if(i<0)
						{
						L.set(index++, o0);
						o0=null;
						}
					else if(i>0)
						{
						L.set(index++, o1);
						o1=null;
						}
					}
				else if(o0==null && o1==null)
					{
					throw new RuntimeException("Boum");
					}
				else if(o0!=null)
					{
					L.set(index++, o0);
					o0=null;
					}
				else if(o1!=null)
					{
					L.set(index++, o1);
					o1=null;
					}
				}
			
			//System.err.println("i0:"+index0+" i1:"+index1+" start:"+start+" end:"+end+" L.size:"+L.size()+" index:"+index);

			start=end;
			}
		
		
		if(!isSorted(L))
			{
			for(int i=0;i< L.size();++i)
				{
				System.err.println("["+i+"/"+L.size()+"]"+L.get(i));
				if(i>0 && compare(L.get(i-1),L.get(i))>0) break;
				}		
			throw new RuntimeException();
			}
		}
	public static void main(String[] args)
		{
		Random r=new Random(0L);
		Algorithms<Integer> a=new Algorithms<Integer>();
		for(;;)
			{
			a.mergeSortSize=2+r.nextInt(2);
			List<Integer> L=new ArrayList<Integer>();
			int N=r.nextInt(10)+10;
			while(L.size()<N)
				{
				int i=r.nextInt(10);
				System.err.println("["+L.size()+"]"+i);
				L.add(i);
				}
			L.add(-1);
			System.out.println(a.isSorted(L));
			a.quickSort(L);
			System.out.println(a.isSorted(L));
			a.shuffle(L);
			System.out.println(a.isSorted(L));
			a.mergeSort(L);
			System.out.println(a.isSorted(L));
			}
		}
	}
