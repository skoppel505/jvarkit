package fr.inserm.umr1087.jvarkit.tools.cnv20120907;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;


import com.sleepycat.bind.tuple.DoubleBinding;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

import fr.inserm.umr1087.jvarkit.bdb.algorithm.StoredListFactory;
import fr.inserm.umr1087.jvarkit.collections.Algorithms;
import fr.inserm.umr1087.jvarkit.collections.DefaultListFactory;
import fr.inserm.umr1087.jvarkit.collections.Listfactory;
import fr.inserm.umr1087.jvarkit.math.Lowess;
import fr.inserm.umr1087.jvarkit.r.RLoess;
import fr.inserm.umr1087.jvarkit.segments.TidStartEnd;
import fr.inserm.umr1087.jvarkit.segments.bdb.TidStartEndBinding;
import fr.inserm.umr1087.jvarkit.segments.bdb.TidStartEndSorter;
import fr.inserm.umr1087.jvarkit.tools.genomegraph.Graphics2DGenomeGraphDrawer;
import fr.inserm.umr1087.jvarkit.util.intervalparser.IntervalParser;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class CNV20120927
	{
	private int BUFFER_SIZE=1000000;//1E6
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
	private int minQual=30;
	private TidStartEndBinding tidStartEndBinding=new TidStartEndBinding();	
	
	
	private ReferenceSequenceFile reference=null;
	private int windowSize=100;
	private int windowStep=50;
	private List<BamBuffer> bams=new ArrayList<BamBuffer>();
	private Interval targetInterval=null;
	private File dbHome=null;
	private Environment environment=null;
	private Database region2depths;
	private Transaction txn=null;
	private Listfactory<Double> listOfDoubleFactory=new DefaultListFactory<Double>();
	
	private class BamBuffer
		{
		File file;
		BamBuffer(File file)
			{
			//super(new SAMFileReader(file));
			
			this.file=file;
			}
		}
	
	private static class XY
		{
		double x;
		double y;
		int index;
		}
	
	private static class XYBinding
		extends TupleBinding<XY>
		{
		@Override
		public XY entryToObject(TupleInput input)
			{
			XY xy=new XY();
			xy.index=input.readInt();
			xy.x=input.readDouble();
			xy.y=input.readDouble();
			return xy;
			}
		@Override
		public void objectToEntry(XY xy, TupleOutput output) {
			output.writeInt(xy.index);
			output.writeDouble(xy.x);
			output.writeDouble(xy.y);
			}
		}
	
	
	private static abstract class AbstractExtractList
		extends AbstractList<Double>
		{
		private List<XY> delegate;
		AbstractExtractList(List<XY> delegate)
			{
			this.delegate=delegate;
			}
		
		protected abstract double extract(XY xy);
		protected abstract double set(XY xy,double v);
		
		@Override
		public Double get(int index) {
			return extract(delegate.get(index));
			}
		
		@Override
		public Double set(int index, Double v)
			{
			XY xy= delegate.get(index);
			double v2=set(xy, v);
			delegate.set(index, xy);
			return v2;
			}
		
		@Override
		public int size() {
			return delegate.size();
			}	
		}
	
	private static class XYListOfX extends AbstractExtractList
		{
		XYListOfX(List<XY> delegate) { super(delegate);}
		@Override
		protected  double extract(XY xy) { return xy.x;}
		@Override
		protected  double set(XY xy,double v) { double old=xy.x; xy.x=v; return old;}
		}
	
	private static class XYListOfY extends AbstractExtractList
		{
		XYListOfY(List<XY> delegate) { super(delegate);}
		@Override
		protected  double extract(XY xy) { return xy.y;}
		@Override
		protected  double set(XY xy,double v) { double old=xy.y; xy.y=v; return old;}
		}



	private CNV20120927()
		{
		
		}

	
	private static double log2(double num)
		{
		return (Math.log(num)/Math.log(2));
		} 
	
	private static double calcMediane(List<Double> L)
		{
		if(L.size()==1) return L.get(0);
		int mid=L.size()/2;
		if(L.size()%2==0)
			{
			return (L.get(mid-1)+L.get(mid))/2.0;
			}
		return L.get(mid);
		}
	
	private void open()
		{
		String dbName1="region2depths";
		EnvironmentConfig envConfig= new EnvironmentConfig();
		envConfig.setAllowCreate(true);
		envConfig.setReadOnly(false);
		envConfig.setConfigParam(EnvironmentConfig.LOG_FILE_MAX,"250000000");
		envConfig.setTransactional(true);
		this.environment= new Environment(dbHome, envConfig);
		this.txn=this.environment.beginTransaction(null, null);
		for(String db:this.environment.getDatabaseNames())
			{
			if(db.equals(dbName1)) this.environment.removeDatabase(txn, db);
			}
		DatabaseConfig cfg= new DatabaseConfig();
		cfg.setAllowCreate(true);
		cfg.setReadOnly(false);
		cfg.setTransactional(true);
		cfg.setBtreeComparator(TidStartEndSorter.class);
		this.region2depths= this.environment.openDatabase(txn,dbName1,cfg);
		
		/*this.listOfDoubleFactory=new StoredListFactory<Double>(
				this.environment,
				this.txn,
				new DoubleBinding()
				);*/
		this.listOfDoubleFactory=new DefaultListFactory<Double>();
		}
	
	
	
	private void close()
		{
		if(listOfDoubleFactory!=null)
			{
			try{this.listOfDoubleFactory.close();} catch(Exception err) {}
			this.listOfDoubleFactory=null;
			}
		if(this.region2depths!=null)
			{
			try{region2depths.close();} catch(Exception err) {}
			this.region2depths=null;
			}
		if(this.txn!=null)
			{
			try{txn.commit();} catch(Exception err) {}
			this.txn=null;
			}
		if(this.environment!=null)
			{
			try{environment.close();} catch(Exception err) {}
			this.environment=null;
			}
		}
	
	
	private void dump(File file) throws IOException
		{
		LOG.info("dump data to "+file);
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		Cursor c=this.region2depths.openCursor(txn, null);
		PrintWriter pw=new PrintWriter(file);
		pw.print("CHROM\tSTART\tEND\tGC");
		for(int bx=0;
		bx< this.bams.size();
		++bx)
			{
			pw.print("\t"+this.bams.get(bx).file.getName().replace(".bam", ""));
			}
		pw.println();
		while(c.getNext(key, data, null)==OperationStatus.SUCCESS)
			{
			TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
			CNVRow row=CNVRow.BINDING.entryToObject(data);
			
			pw.print(
					this.reference.getSequenceDictionary().getSequence(seg.getChromId()).getSequenceName()+
					"\t"+seg.getStart()+
					"\t"+seg.getEnd()+
					"\t"+row.getGcPercent()
					);
			for(int i=0;i< this.bams.size();++i)
				{
				pw.print("\t"+row.getDepth(i));
				}
			pw.println();
				
			}
		c.close();
		pw.flush();
		pw.close();
		}
	
	private void dump(File file,int bamIndex) throws IOException
		{
		LOG.info("dump data to "+file);
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		int nRows=0;
		Cursor c=this.region2depths.openCursor(txn, null);
		PrintWriter pw=new PrintWriter(file);
		pw.println();
		while(c.getNext(key, data, null)==OperationStatus.SUCCESS && nRows++<20)
			{
			TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
			CNVRow row=CNVRow.BINDING.entryToObject(data);
			
			pw.print(
					this.reference.getSequenceDictionary().getSequence(seg.getChromId()).getSequenceName()+
					"\t"+seg.getStart()
					);
		
			pw.print("\t"+row.getDepth(bamIndex));
				
			pw.println();
			}
		c.close();
		pw.flush();
		pw.close();
		}
	
	
	@SuppressWarnings("unused")
	private void dumpGnuplot(File file) throws IOException
		{
		
		ZipOutputStream zout=new ZipOutputStream(new FileOutputStream(file));
		
		for(int bx=0;
			bx< this.bams.size();
			++bx)
			{
			long shift=0L;
			ZipEntry zentry=new ZipEntry(this.bams.get(bx).file.getName()+".tsv");
			zout.putNextEntry(zentry);
			PrintStream pout=new PrintStream(zout);
			for(SAMSequenceRecord chrom: this.reference.getSequenceDictionary().getSequences())
				{
				if(targetInterval!=null && !chrom.getSequenceName().equals(targetInterval.getSequence())) continue;
				long chromMax=0L;
				DatabaseEntry key=new DatabaseEntry();
				DatabaseEntry data =new DatabaseEntry();
				Cursor c=this.region2depths.openCursor(txn, null);
				while(c.getNext(key, data, null)==OperationStatus.SUCCESS)
					{
					
					TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
					if(seg.getChromId()!=chrom.getSequenceIndex())
						{
						continue;
						}
					int chromStart=seg.getStart();				
					CNVRow row=CNVRow.BINDING.entryToObject(data);
					if(row.containsNanDepth()) continue;
					chromMax=Math.max(chromMax, seg.getEnd());
				    double y=row.getDepth(bx);
					pout.println(""+(shift+chromStart)+"\t"+y);
					}
				c.close();
				shift+=chromMax;
				}
			pout.flush();
			zout.closeEntry();
			}
		
		
		ZipEntry zentry=new ZipEntry("jeter.gnuplot");
		zout.putNextEntry(zentry);
		PrintStream gnuplotout=new PrintStream(zout);
		gnuplotout.println("set term postscript");
		gnuplotout.println("set output jeter.ps");
		
		for(int bx=0;
		bx< this.bams.size();
		++bx)
			{
			gnuplotout.println("set title \"SAMPLE-"+(1+bx)+" "+this.bams.get(bx).file.getName()+"\"");
			gnuplotout.println("set xlabel \"Position\"");
			gnuplotout.println("set ylabel \"Depth\"");
			gnuplotout.println("set yrange [-3:3]");
			gnuplotout.println("plot \""+ this.bams.get(bx).file.getName()+".tsv\" using 1:2 notitle\n");
			}
		gnuplotout.flush();
		zout.closeEntry();
		
		zout.finish();
		zout.close();
		}
	
	private abstract class ForEachRow
		{
		private Cursor c;
		DatabaseEntry key;
		DatabaseEntry data;
		int bamIndex=-1;
		TidStartEnd currSeg;
		CNVRow currRow;
		int rowIndex;

		protected ForEachRow(int bamIndex)
			{
			this.bamIndex=bamIndex;
			}
		
		public void scanChromosome(int tid) throws Exception
			{
			LOG.info("Scan bam:"+this.bamIndex+" tid:"+tid+" class:"+getClass());
			key=new DatabaseEntry();
			data=new DatabaseEntry();

			this.c=CNV20120927.this.region2depths.openCursor(txn, null);
			boolean first=true;
			this.rowIndex=-1;
			for(;;)
				{
				OperationStatus status;
				if(first) 
					{
					first=false;
					if(tid==-1)
						{
						status=c.getFirst(key, data,null);
						}
					else
						{
						TidStartEnd init=new TidStartEnd(tid,
								0,
								0
								);
					
						CNV20120927.this.tidStartEndBinding.objectToEntry(init, key);
						status=c.getSearchKeyRange(key, data,null);
						
						}
					}
				else
					{
					status=c.getNext(key, data,null);
					}
				if(status!=OperationStatus.SUCCESS) break;
				this.currSeg=CNV20120927.this.tidStartEndBinding.entryToObject(key);
				if(tid!=-1)
					{
					if(this.currSeg.getChromId()<tid) continue;
					if(this.currSeg.getChromId()>tid) break;
					}
				++rowIndex;

				this.currRow=CNVRow.BINDING.entryToObject(data);
				apply();
				}
			c.close();
			}
		public void update()
			{
			CNVRow.BINDING.objectToEntry(this.currRow,data);
			if(c.putCurrent(data)!=OperationStatus.SUCCESS)
				{
				c.close();
				close();
				throw new RuntimeException("Cannot update "+this.currRow);
				}
			}
		public abstract void apply() throws Exception;
		}
	
	private class ForEachRowAccumulateData
		extends ForEachRow
		{
		List<Double> vData;
		ForEachRowAccumulateData(int bamIndex)
			{
			super(bamIndex);
			vData=CNV20120927.this.listOfDoubleFactory.createList(10000);
			}
		@Override
		public void scanChromosome(int tid) throws Exception
			{
			vData.clear();
			super.scanChromosome(tid);
			}
		@Override
		public void apply() throws Exception {
			vData.add(currRow.getDepth(this.bamIndex));
			}
		}
	
	private class ForEachRowCalcDepth
		extends ForEachRow
		{
		short depth0[];
		ForEachRowCalcDepth(int bamIndex,short depth0[])
			{
			super(bamIndex);
			this.depth0=depth0;
			}
		@Override
		public void apply() throws Exception
			{
			int count=0;
			double sum=0.0;
			for(int i=this.currSeg.getStart();i<this.currSeg.getEnd() &&
					i < depth0.length;
					++i)
				{
				sum+=depth0[i];
				++count;
				}
			
			
			currRow.setDepth(this.bamIndex, (count==0?Double.NaN:sum/count));
			update();
			}
		}
	
	/****************************************************************************************
	 * 
	 * AbstractForEachRowLoess
	 *
	 */
	private abstract class AbstractForEachRowLoess
		extends ForEachRow
		{
		List<Double> vData=null;
		AbstractForEachRowLoess(int bamIndex)
			{
			super(bamIndex);
			}
		
		}
	
	
	/****************************************************************************************
	 * 
	 * ForEachRowCalcLoess
	 *
	 */
	private class ForEachRowCalcLoess
		extends AbstractForEachRowLoess
		{
		private RLoess rLoessProc;
		ForEachRowCalcLoess(int bamIndex)
			{
			super(bamIndex);
			}
		@Override
		public void scanChromosome(int tid) throws Exception
			{
			LOG.info("scanning loess for bamidx="+this.bamIndex+" tid:"+tid);
			this.rLoessProc=new RLoess();
			this.rLoessProc.setListOfDoubleFactory(CNV20120927.this.listOfDoubleFactory);
			super.scanChromosome(tid);
			super.vData=rLoessProc.smooth();
			LOG.info("CalcLoess: vData.size="+vData.size()+" rowIndex:"+this.rowIndex);
			}
		@Override
		public void apply() throws Exception
			{
			rLoessProc.add(
					this.currRow.getGcPercent(),
					this.currRow.getDepth(super.bamIndex)
					);
			}
		}
	
	/****************************************************************************************
	 * 
	 * ForEachRowCalcLoess
	 *
	 */
	private class ForEachRowCalcLowess
		extends AbstractForEachRowLoess
		{
		private Lowess lowess;
		private List<XY> points;
		ForEachRowCalcLowess(int bamIndex)
			{
			super(bamIndex);
			}
		private void sortPoints(Comparator<XY> cmp)
			{
			LOG.info("sort "+points.size());
			Algorithms<XY> algos=new Algorithms<CNV20120927.XY>(cmp);
			algos.quickSort(this.points);
			LOG.info("sorted");
			}
		@Override
		public void scanChromosome(int tid) throws Exception
			{
			LOG.info("scanning loess for bamidx="+this.bamIndex+" tid:"+tid);
			
			Listfactory<XY> listOfPointsfactory=new DefaultListFactory<XY>();
			/*
			XYBinding binding=new XYBinding();
			new StoredListFactory<CNV20120927.XY>(
					CNV20120927.this.environment,
					CNV20120927.this.txn,
					binding
					);*/
			
			this.lowess=new Lowess();
			this.lowess.setListFactory(CNV20120927.this.listOfDoubleFactory);
			
			points=listOfPointsfactory.createList(10000);
			
			
			super.scanChromosome(tid);
			
			
			sortPoints(new Comparator<XY>()
						{
						@Override
						public int compare(XY o1, XY o2)
							{
							if(o1.x==o2.x)
								{
								if(o1.y==o2.y) return 0;
								return o1.y < o2.y?-1:1;
								}
							return o1.x < o2.x?-1:1;
							}
						});

			List<Double> var_x= new XYListOfX(this.points);
			List<Double> var_y= new XYListOfY(this.points);

			for(int i=0;i< var_y.size();++i)
				{
				if(i<10) System.err.println(""+i+":"+var_y.get(i));
				}
			
			super.vData=this.lowess.lowess(
				var_x,
				var_y,
				var_x.size()
				);
			
			for(int i=0;i< super.vData.size();++i)
				{
				XY xy=this.points.get(i);
				xy.y=super.vData.get(i);
				if(i<10) System.err.println(""+i+":"+xy.y);
				this.points.set(i,xy);
				}
			
			
			sortPoints(new Comparator<XY>()
					{
					@Override
					public int compare(XY o1, XY o2)
						{
						return o1.index -o2.index;
						}
					});
			
			for(int i=0;i< points.size();++i)
				{
				XY xy=this.points.get(i);
			
				super.vData.set(i,xy.y);
				}
			
			listOfPointsfactory.disposeList(points);
			LOG.info("CalcLoess: vData.size="+vData.size()+" rowIndex:"+this.rowIndex);
			}
		@Override
		public void apply() throws Exception
			{
			XY xy=new XY();
			xy.x=this.currRow.getGcPercent();
			xy.y=this.currRow.getDepth(super.bamIndex);
			xy.index=this.rowIndex;
			this.points.add(xy);
			}
		}
	
	
	/****************************************************************************************
	 * 
	 * ForEachRowApplyLoessAndGetMin
	 *
	 */
	private class ForEachRowApplyLoessAndGetMin
	extends ForEachRow
		{
		private List<Double> vData;
		double minDepth;
		ForEachRowApplyLoessAndGetMin(int bamIndex,List<Double> vData)
			{
			super(bamIndex);
			this.vData=vData;
			}
		@Override
		public void scanChromosome(int tid) throws Exception
			{
			this.minDepth=Double.MAX_VALUE;
			super.scanChromosome(tid);
			LOG.info("min depth:"+minDepth+" tid:"+tid+" bamIdx:"+this.bamIndex+" vsdata.size:"+this.vData.size());
			CNV20120927.this.listOfDoubleFactory.disposeList(this.vData);
			}
		@Override
		public void apply() throws Exception
			{
			double newDepth=vData.get(this.rowIndex);
			if(minDepth> newDepth) minDepth=newDepth;
			currRow.setDepth(this.bamIndex, newDepth);
			update();
			}
		}
	
	/****************************************************************************************
	 * 
	 * ForEachRowSubstractMinAndGetMedian
	 *
	 */
	private class ForEachRowSubstractMinAndGetMedian
	extends ForEachRow
		{
		private double minDepth;
		double medianValue;
		private List<Double> _vData;
		ForEachRowSubstractMinAndGetMedian(int bamIndex,double minDepth)
			{
			super(bamIndex);
			this.minDepth=minDepth;
			}
		@Override
		public void scanChromosome(int tid) throws Exception
			{
			_vData=CNV20120927.this.listOfDoubleFactory.createList(10000);
			super.scanChromosome(tid);
			Collections.sort(_vData);
			medianValue = calcMediane(this._vData);
			CNV20120927.this.listOfDoubleFactory.disposeList(_vData);
			LOG.info("medianValue depth:"+medianValue+" tid:"+tid+" bamIdx:"+this.bamIndex+" vdata.size:"+
					_vData.size()+" rowIndex:"+this.rowIndex
					);
			}
		@Override
		public void apply() throws Exception
			{
			double depth=currRow.getDepth(this.bamIndex);
			double newDepth=depth-this.minDepth;
			this._vData.add(newDepth);
			currRow.setDepth(this.bamIndex, newDepth);
			update();
			}
		}
	
	
	/****************************************************************************************
	 * 
	 * ForEachRowSubstractMinAndGetMedian
	 *
	 */
	private class ForEachRowDivMedian
	extends ForEachRow
		{
		private double medianValue;
		ForEachRowDivMedian(int bamIndex,double medianValue)
			{
			super(bamIndex);
			this.medianValue=medianValue;
			}
		
		@Override
		public void apply() throws Exception
			{
			double depth=currRow.getDepth(this.bamIndex);
			double newDepth=depth/medianValue;
			currRow.setDepth(this.bamIndex, newDepth);
			update();
			}
		}

	
	/****************************************************************************************
	 * 
	 * ForEachRowHorizontalMedian
	 *
	 */
	private class ForEachRowHorizontalMedian
		extends ForEachRow
		{
		double cutOff;
		ForEachRowHorizontalMedian(double cutOff)
			{
			super(-1);
			this.cutOff=cutOff;
			}
		
		@Override
		public void apply() throws Exception
			{
			List<Double> v=new ArrayList<Double>(currRow.getSampleCount()-1);
			for(int i=1;i<currRow.getSampleCount();++i)
				{
				v.add(currRow.getDepth(i));
				if(this.rowIndex==0)
					{
					System.err.println("###["+i+"]="+v.get(v.size()-1));
					}
				}
			Collections.sort(v);
			double mediane=calcMediane(v);
			if(this.rowIndex==0)
				{
				System.err.println("###mediane="+mediane);
				}
			if(mediane>this.cutOff)
				{
				for(int i=0;i< currRow.getSampleCount();++i)
					{
					if(this.rowIndex==0)
						{
						System.err.println("###new["+i+"]="+currRow.getDepth(i)/mediane);
						}
					
					currRow.setDepth(i, currRow.getDepth(i)/mediane);
					}
				update();
				}
			else
				{
				super.c.delete();
				}
			}
		}
	
	private class ForEachRowMed
	extends ForEachRow
		{
		private List<Double> vData;
		private int margin;
		ForEachRowMed(int bamIndex,List<Double> vData,int margin)
			{
			super(bamIndex);
			this.vData=vData;
			this.margin=margin;
			}
		@Override
		public void scanChromosome(int tid) throws Exception {
			super.scanChromosome(tid);
			CNV20120927.this.listOfDoubleFactory.disposeList(vData);
			}
		@Override
		public void apply() throws Exception
			{
			double newDepth=0;
			
			if(vData.get(rowIndex)!=currRow.getDepth(bamIndex))
				{
				throw new RuntimeException("Err "+vData.get(rowIndex)+"/"+currRow.getDepth(bamIndex));
				}
			List<Double> runmed=new ArrayList<Double>(margin*2+1);
			for(int i=Math.max(0,rowIndex-margin);
					i<= rowIndex+margin && i< this.vData.size();
					++i)
				{
				runmed.add(vData.get(i));
				}
			Collections.sort(runmed);
			newDepth=log2(calcMediane(runmed));
			if(Double.isNaN(newDepth)) newDepth=0;
			
			currRow.setDepth(this.bamIndex, newDepth);
			update();
			}
		}
	
	private void run() throws Exception
		{
		File baseDir=new File("/commun/data/users/lindenb/_ignore.backup/");
		open();
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		SAMSequenceDictionary 	dict=this.reference.getSequenceDictionary();
		/** loop over each chromosome */
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			if(targetInterval!=null && !targetInterval.getSequence().equals(chrom.getSequenceName()))
				{
				continue;
				}
			/** loop over the bases of the chromosomes */
			byte bases[]=this.reference.getSequence(chrom.getSequenceName()).getBases();
			int  start0=(targetInterval!=null?targetInterval.getStart():0);
			while(start0+this.windowSize<
					(targetInterval!=null?targetInterval.getEnd():bases.length))
				{
				int gc=0;
				int N=0;
				for(int i=start0;N==0 && i<start0+this.windowSize;++i)
					{
					byte base=bases[i];
					switch(base)
						{
						case 'n': case 'N': ++N;break;
						case 's': case 'S':
						case 'c': case 'C':
						case 'g': case 'G': ++gc; break;
						default:break;
						}
					}
				if(N>0)
					{
					/* if there is a 'N' in the window, skip to next base */
					++start0;
					continue;
					}
				TidStartEnd tidStartEnd=new TidStartEnd(
						chrom.getSequenceIndex(),
						start0,
						start0+windowSize
						);
				CNVRow row=new CNVRow(this.bams.size());
				row.setGcPercent(gc/(double)windowSize);
				
				
				
				
				tidStartEndBinding.objectToEntry(tidStartEnd, key);
				CNVRow.BINDING.objectToEntry(row, data);
				
				if(this.region2depths.put(txn, key, data)!=OperationStatus.SUCCESS)
					{
					close();
					throw new RuntimeException("Cannot insert");
					}

				start0+=this.windowStep;
				}
			bases=null;
			
			/* get the depth for all bams */
			short depth0[]=new short[chrom.getSequenceLength()];
			for(int bx=0;
				bx< this.bams.size();
				++bx)
					{
					LOG.info("Getting depth for "+chrom.getSequenceName()+" "+this.bams.get(bx).file);
					Arrays.fill(depth0, (short)0);
					SAMFileReader sf=new SAMFileReader(this.bams.get(bx).file);
					sf.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
					
					/* create a SamLocusIterator to get the DEPTH at each position */
					Interval interval;
					if(targetInterval==null)
						{
						interval=new Interval(
			            		chrom.getSequenceName(),
			            		1,
			            		depth0.length
			            		);//Coordinates are 1-based closed ended.
						}
					else
						{
						interval=new Interval(
			            		chrom.getSequenceName(),
			            		Math.max(1, this.targetInterval.getStart()),
			            		Math.min(this.targetInterval.getEnd(),depth0.length)
			            		);//Coordinates are 1-based closed ended.
						}
					IntervalList iL=new IntervalList(sf.getFileHeader());
		            iL.add(interval);
		            SamLocusIterator sli=new SamLocusIterator(sf,iL,true);           
		            sli.setMappingQualityScoreCutoff(this.minQual);
		            for(Iterator<SamLocusIterator.LocusInfo>  iter=sli.iterator();
		                iter.hasNext();
		                )
		                {
		                SamLocusIterator.LocusInfo locusInfo=iter.next();
		                int heredepth=locusInfo.getRecordAndPositions().size();
		                if(heredepth>Short.MAX_VALUE)
		                	{
		                	LOG.info("WARNING depth >"+Short.MAX_VALUE+" at "+locusInfo.getPosition());
		                	}
		                depth0[locusInfo.getPosition()-1]=(short)Math.min(heredepth, Short.MAX_VALUE);
		                }
		            
		            sf.close();
		            ForEachRow forEach=new ForEachRowCalcDepth(bx,depth0);
		            forEach.scanChromosome(chrom.getSequenceIndex());
		            forEach=null;
		            
		           
					}
			depth0=null;
			}
		
		dump(new File(baseDir,"jeter.tsv"));
		
		for(SAMSequenceRecord chrom: dict.getSequences())
			{
			if(targetInterval!=null && !chrom.getSequenceName().equals(targetInterval.getSequence())) continue;
			int tid=chrom.getSequenceIndex();
			for(int bx=0;
					bx< this.bams.size();
					++bx)
				{
				boolean uselowess=false;
				AbstractForEachRowLoess calcLoess=null;
				calcLoess=(uselowess?
						new ForEachRowCalcLowess(bx):
						new ForEachRowCalcLoess(bx));
				calcLoess.scanChromosome(tid);
								
				ForEachRowApplyLoessAndGetMin applyLoess=new ForEachRowApplyLoessAndGetMin(bx, calcLoess.vData);
				
				calcLoess=null;
				applyLoess.scanChromosome(tid);
				
				
				//dump(new File(baseDir,"jeter.loess."+bx+".tsv"),bx);
				
				LOG.info("min depth:"+applyLoess.minDepth);
				double minDepth=applyLoess.minDepth;
				applyLoess=null;
				
				ForEachRowSubstractMinAndGetMedian subMinMedian=new ForEachRowSubstractMinAndGetMedian(bx, minDepth);
				subMinMedian.scanChromosome(tid);
				double median=subMinMedian.medianValue;
				subMinMedian=null;
				LOG.info("median:"+subMinMedian);
				//dump(new File(baseDir,"jeter.min."+bx+".tsv"),bx);
				
				ForEachRowDivMedian divMedian=new ForEachRowDivMedian(bx,median);
				divMedian.scanChromosome(tid);
				
				//dump(new File(baseDir,"jeter.divmedian."+bx+".tsv"),bx);
				}
			
			ForEachRowHorizontalMedian forEachHorizontalMedian=
				new ForEachRowHorizontalMedian(
						0.2
						);
			forEachHorizontalMedian.scanChromosome(tid);
			forEachHorizontalMedian=null;
			
			for(int bx=0;
			bx< this.bams.size();
			++bx)
				{
				//dump(new File(baseDir,"jeter.hmedian."+bx+".tsv"),bx);
				}
			
			for(int bx=0;
			bx< this.bams.size();
			++bx)
				{
				ForEachRowAccumulateData forEachVData=new ForEachRowAccumulateData(bx);
				forEachVData.scanChromosome(tid);
				List<Double> vData=forEachVData.vData;
				forEachVData=null;
				ForEachRowMed forEachMed=new ForEachRowMed(bx, vData, 3);
				forEachMed.scanChromosome(tid);
				//dump(new File(baseDir,"jeter.med."+bx+".tsv"),bx);
				}
			}
		
		PrintWriter html=new PrintWriter(new File("/commun/data/users/lindenb/_ignore.backup/jeter.html"));
		html.print("<html><body>");
		
		Cursor c=this.region2depths.openCursor(txn, null);
		OperationStatus status;

		for(int bx=0;
			bx< this.bams.size();
			++bx)
			{
			Graphics2DGenomeGraphDrawer drawer=new Graphics2DGenomeGraphDrawer(this.reference.getSequenceDictionary(),this.targetInterval);
			drawer.setMinY(-3.0);
			drawer.setMaxY(3.0);
			BufferedImage img=drawer.createImage(BufferedImage.TYPE_INT_RGB);
			Graphics2D g=(Graphics2D)img.getGraphics();
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			drawer.init(g);
			
			File imgFile=new File("/commun/data/users/lindenb/_ignore.backup/"+this.bams.get(bx).file.getName()+".jpg");
			html.println("<div>");
			html.println("<h1>"+imgFile.getName()+"</h1>");

			key=new DatabaseEntry();
			status=c.getFirst(key, data, null);
			
			while(status==OperationStatus.SUCCESS)
				{
				
				TidStartEnd seg=this.tidStartEndBinding.entryToObject(key);
				int tid=seg.getChromId();
				int chromStart=seg.getStart();
				int chromEnd=seg.getEnd();
				
				CNVRow row=CNVRow.BINDING.entryToObject(data);
				if(!row.containsNanDepth())
					{
				    double y=row.getDepth(bx);
				    double pix_y= drawer.convertValueToPixelY(y);
				    double pix_x= drawer.convertPositionToPixelX(tid, (chromStart+chromEnd)/2);
				    Composite oldComposite=g.getComposite();
				    g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,0.2f));
				    g.setColor(Color.RED);
				    g.fillOval((int)pix_x-2,(int)pix_y-2,5,5);
				    g.setComposite(oldComposite);
					status=c.getNext(key, data, null);
					}
				}
			
			drawer.finish();
			g.dispose();
			ImageIO.write(img, "JPG", imgFile);
			html.println("<img src='"+imgFile.getName()+"'/>");
			html.println("</div>");
			}
		c.close();
		
		html.print("</body></html>");
		html.flush();
		html.close();
		dump(new File(baseDir,"jeter.end.tsv"));
		
		}
	
	private void usage()
		{
		System.err.println(
			"Options:\n"+
			" --db-home (dir) berkeleydb-home\n"+
			" -f (fasta) reference fasta indexed with faidx\n"+
			" -q <int> quality treshold for BAM alignments\n"+
			" -L <Chrom:star-end> target interval (optional)\n"+
			" -w <int> window-size["+this.windowSize+"]\n"+
			" -s <int> step-size["+this.windowStep+"]\n"
			);
		}

	private void run(String[] args) throws Exception
		{
		
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				usage();
				return;
				}
			else if(args[optind].equals("--db-home") && optind+1 < args.length )
				{
				this.dbHome=new File(args[++optind]);
				}
			else if(args[optind].equals("-f") && optind+1 < args.length )
				{
				this.reference= ReferenceSequenceFileFactory.getReferenceSequenceFile(
						new File(args[++optind]),
						true //cut after first whitespace
						);
				}
			else if(args[optind].equals("-q") && optind+1 < args.length )
				{
				this.minQual=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-w") && optind+1 < args.length )
				{
				this.windowSize=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-s") && optind+1 < args.length )
				{
				this.windowStep=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("-L") && optind+1 < args.length )
				{
				this.targetInterval=IntervalParser.parseOne(args[++optind]);
				if(this.targetInterval==null)
					{
					System.err.println("bad range:"+args[optind]);
					return;
					}
				//int len=this.targetInterval.getEnd()-this.targetInterval.getStart();
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return;
				}
			else
				{
				break;
				}
			++optind;
			}
		if(BUFFER_SIZE<=this.windowSize)
			{
			BUFFER_SIZE=(this.windowSize+1000);
			}
		
		if(this.dbHome==null)
			{
			System.err.println("--db-home missing");
			return;
			}
		
		if(this.reference==null)
			{
			System.err.println("Reference missing");
			return;
			}
		//this.referenceBuffer=new ReferenceBuffer(this.reference);
		//this.referenceBuffer.setBufferSize(BUFFER_SIZE);

		if(optind==args.length)
			{
			System.err.println("No bam defined");
			return;
			}
		SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

		while(optind!=args.length)
			{
			BamBuffer buf=new BamBuffer(new File(args[optind++]));
			//buf.setBufferSize(BUFFER_SIZE);
			//buf.setMinQual(this.minQual);
			this.bams.add(buf);
			}
		//this.test();
		this.run();
		}

	public static void main(String[] args)
		{
		try
			{
			new CNV20120927().run(args);
			}
		catch (Exception e)
			{
			e.printStackTrace();
			}
		}
}
