package fr.inserm.umr1087.jvarkit.collections;

import java.util.ArrayList;
import java.util.List;

public class DefaultListFactory<T> implements Listfactory<T> {
	@Override
	public List<T> createList(int capacity)
		{
		return new ArrayList<T>(capacity);
		}
	@Override
	public void disposeList(List<T> list)
		{
		if(list==null) return;
		list.clear();
		}
	@Override
	public void close() {
		
		}
	}
