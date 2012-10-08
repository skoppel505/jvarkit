package fr.inserm.umr1087.jvarkit.bdb.algorithm;

import java.util.AbstractList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import com.sleepycat.bind.tuple.IntegerBinding;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

import fr.inserm.umr1087.jvarkit.collections.Listfactory;

public class StoredListFactory<T>
	implements Listfactory<T>
	{
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
	private static long ID_GENERATOR=0L;
	private TupleBinding<T> binding;
	private Environment environment;
	private Transaction txn;
	private Map<String, MyStoredList> id2storedList=new HashMap<String,MyStoredList>();
	
	private class MyStoredList
		extends AbstractList<T>
		{
		Database database;
		String id;
		int _size=0;
		public MyStoredList()
			{
			DatabaseConfig cfg=new DatabaseConfig();
			cfg.setAllowCreate(true);
			cfg.setReadOnly(false);
			cfg.setTransactional(true);
			this.id="stored.list."+(++ID_GENERATOR);
			this.database=environment.openDatabase(StoredListFactory.this.txn, id, cfg);
			}
		
		@Override
		public void clear()
			{
			DatabaseEntry key=new DatabaseEntry();
			DatabaseEntry data=new DatabaseEntry();
			Cursor c=database.openCursor(txn, null);
			while(c.getNext(key, data, null)==OperationStatus.SUCCESS)
				{
				c.delete();
				}
			c.close();
			_size=0;
			}
		
		@Override
		public boolean add(T e)
			{
			DatabaseEntry key=new DatabaseEntry();
			DatabaseEntry data=new DatabaseEntry();
			IntegerBinding.intToEntry(_size, key);
			getBinding().objectToEntry(e, data);
			if(this.database.put(txn, key, data)!=OperationStatus.SUCCESS)
				{
				throw new IndexOutOfBoundsException("cannot insert");
				}
			++_size;
			return true;
			}
		@Override
		public T set(int index, T e)
			{
			if(index<0 || index>=size())
				{
				throw new IndexOutOfBoundsException("0<="+index+"<"+size());
				}
			DatabaseEntry key=new DatabaseEntry();
			DatabaseEntry data=new DatabaseEntry();
			IntegerBinding.intToEntry(index, key);
			getBinding().objectToEntry(e, data);
			if(this.database.put(txn, key, data)!=OperationStatus.SUCCESS)
				{
				throw new IndexOutOfBoundsException("cannot insert");
				}
			return null;//todo
			}

		@Override
		public T get(int index)
			{
			DatabaseEntry key=new DatabaseEntry();
			DatabaseEntry data=new DatabaseEntry();
			IntegerBinding.intToEntry(index, key);
			if(this.database.get(txn, key, data, null)!=OperationStatus.SUCCESS)
				{
				throw new IndexOutOfBoundsException("0<="+index+"<"+size());
				}
			return getBinding().entryToObject(data);
			}

		@Override
		public int size()
			{
			return _size;
			}
		
		
		}
	
	public StoredListFactory( Environment environment,Transaction txn, TupleBinding<T> binding)
		{
		this.environment=environment;
		this.binding=binding;
		this.txn=txn;
		}
	
	public TupleBinding<T> getBinding() {
		return binding;
		}
	@Override
	public List<T> createList(int capacity)
		{
		MyStoredList L=new MyStoredList();
		id2storedList.put(L.id, L);
		return L;
		}
	@Override
	public void disposeList(List<T> list)
		{
		MyStoredList L=(MyStoredList)list;
		L.database.close();
		this.environment.removeDatabase(this.txn, L.id);
		id2storedList.remove(L.id);
		this.environment.cleanLog();
		}
	@Override
	public void close()
		{
		for(String id:id2storedList.keySet())
			{
			LOG.info("database "+id+" wasn't disposed");
			}
		}
	
	}
