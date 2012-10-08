package fr.inserm.umr1087.jvarkit.collections;

import java.util.List;

public interface Listfactory<T>
	{
	public List<T> createList(int capacity);
	public void disposeList(List<T> list);
	public void close();
	}
