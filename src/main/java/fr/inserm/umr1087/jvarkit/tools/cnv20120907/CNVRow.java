package fr.inserm.umr1087.jvarkit.tools.cnv20120907;

import java.util.Arrays;

import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

public class CNVRow
	{
	private double gcPercent;
	private double depths[];
	public static class Binding
		extends TupleBinding<CNVRow>
		{
		@Override
		public CNVRow entryToObject(TupleInput in)
			{
			int nSamples=in.readInt();
			CNVRow o=new CNVRow(nSamples);
			o.gcPercent=in.readDouble();
			for(int i=0;i< nSamples;++i)
				{
				o.depths[i]=in.readDouble();
				}
			return o;
			}
		@Override
		public void objectToEntry(CNVRow o, TupleOutput out)
			{
			out.writeInt(o.depths.length);
			out.writeDouble(o.gcPercent);
			for(double d: o.depths)
				{
				out.writeDouble(d);
				}
			}
		}
	
	public static final TupleBinding<CNVRow> BINDING=new Binding();
	
	public CNVRow(int nSamples)
		{
		this.depths=new double[nSamples];
		}
	public void setGcPercent(double gcPercent)
		{
		this.gcPercent = gcPercent;
		}
	public double getGcPercent()
		{
		return gcPercent;
		}
	
	public void setDepth(int index,double d)
		{
		this.depths[index] = d;
		}
	
	public double getDepth(int index)
		{
		return depths[index];
		}
	
	public boolean containsNanDepth()
		{
		for(double d: this.depths) if(Double.isNaN(d)) return true;
		return false;
		}
	
	public int getSampleCount()
		{
		return this.depths.length;
		}
	
	public double getMedianDepth()
		{
		double v[]=new double[this.depths.length];
		System.arraycopy(v, 0, this.depths,0,this.depths.length);
		Arrays.sort(v);
		return v[v.length/2];
		}
	
	@Override
	public String toString() {
		return "gc:"+gcPercent+" depths:"+Arrays.toString(this.depths);
		}
	}
