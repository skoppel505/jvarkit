package fr.inserm.umr1087.jvarkit.tools.genomegraph;

import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import net.sf.picard.util.Interval;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;



public abstract class AbstractGenomeGraphDrawer
	{
	protected static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");

	protected class Bound
		{
		SAMSequenceRecord samSeqRec;
		int chromStart;
		int chromEnd;
		double pix_x;
		double pix_length;
		int index=0;
		public Rectangle2D.Double getSize()
			{
			return new Rectangle2D.Double(
				this.getX(),
				this.getY(),
				this.getWidth(),
				this.getHeight()
				);
			}
		double getX()
			{
			return pix_x;
			}
		double getWidth()
			{
			return pix_length;
			}
		double getY()
			{
			return getDrawingBounds().getY();
			}
		double getHeight()
			{
			return getDrawingBounds().getHeight();
			}
		
		public long getMinBase()
			{
			return  chromStart;
			}
		
		public long getLengthBp()
			{
			
			return  chromEnd-chromStart;
			}
		
		double convertBaseToPixelX(int position)
			{
			return getX() + ((position-this.getMinBase())/(double)getLengthBp())*getWidth();
			}
		@Override
		public String toString() {
			return samSeqRec.getSequenceName()+" "+getSize();
			}
		}
	protected SAMSequenceDictionary samSequenceDictionary;
	private Interval interval;
	protected List<Bound> boundaries=new ArrayList<Bound>();
	private Map<String,Bound> name2bound=new HashMap<String, Bound>();
	private Map<Integer,Bound> tidbound=new HashMap<Integer, Bound>();
	protected long total_bases=0L;
	private Insets margin=new Insets(50, 200, 100, 100);
	private Dimension dimension=new Dimension(1000,400);
	private double distBetweenChrom=10;
	private Double userMinY=null;
	private Double userMaxY=null;
	
	
	protected AbstractGenomeGraphDrawer(SAMSequenceDictionary samSequenceDictionary)
		{
		this(samSequenceDictionary,null);
		}
	protected AbstractGenomeGraphDrawer(SAMSequenceDictionary samSequenceDictionary,Interval interval)
		{
		this.samSequenceDictionary=samSequenceDictionary;
		this.interval=interval;
		}
	
	public Dimension getSize()
		{
		return dimension;
		}
	
	public int getWidth()
		{
		return getSize().width;
		}
	
	public int getHeight()
		{
		return getSize().height;
		}

	public double getMinY()
		{
		return this.userMinY==null?-1.0:this.userMinY;
		}
	
	public double getMaxY()
		{
		return this.userMaxY==null?1.0:this.userMaxY;
		}
	
	public void setMaxY(Double userMaxY)
		{
		this.userMaxY = userMaxY;
		}
	
	public void setMinY(Double userMinY)
		{
		this.userMinY = userMinY;
		}
	
	public Insets getInsets()
		{
		return margin;
		}
	
	protected Rectangle getDrawingBounds()
		{
		return new Rectangle(
			this.margin.left,
			this.margin.top,
			getWidth()-(this.margin.left+this.margin.right),
			getHeight()-(this.margin.top+this.margin.bottom)
			);
		}
	
	public double convertValueToPixelY(double value)
		{
		Rectangle r=getDrawingBounds();
		return r.y+r.height*(1.0 -(value-getMinY())/(getMaxY()-getMinY()));
		}
	
	public double convertPositionToPixelX(String chrom,int position)
		{
		Bound b=this.name2bound.get(chrom);
		if(b==null)
			{
			LOG.info("unknown chromosome "+chrom);
			return -1;
			}
		return b.convertBaseToPixelX(position);
		}

	public double convertPositionToPixelX(int tid,int position)
		{
		Bound b=this.tidbound.get(tid);
		if(b==null)
			{
			LOG.info("unknown chromosome "+tid);
			return -1;
			}
		return b.convertBaseToPixelX(position);
		}
	
	public void init()
		{
		int countChromosomes=0;
		this.boundaries.clear();
		this.total_bases=0L;
		this.tidbound.clear();
		this.name2bound.clear();
		for(SAMSequenceRecord ci: this.samSequenceDictionary.getSequences())
			{
			if(this.interval!=null && !ci.getSequenceName().equals(interval.getSequence())) continue;
			++countChromosomes;
			if(this.interval!=null)
				{
				this.total_bases+= 1+(this.interval.getEnd()-this.interval.getStart());
				}
			else
				{
				this.total_bases+= ci.getSequenceLength();
				}
			}
		double pixel_width_for_chromosomes= getDrawingBounds().getWidth() - (
				countChromosomes==0?
				0:
				(countChromosomes-1)*this.distBetweenChrom
				);
		double x=margin.left;
		for(SAMSequenceRecord ci: this.samSequenceDictionary.getSequences())
			{
			if(this.interval!=null && !ci.getSequenceName().equals(interval.getSequence())) continue;
			Bound b=new Bound();
			b.samSeqRec=ci;
			b.pix_x=x;
			
			double L;
			
			if(this.interval!=null)
				{
				b.chromStart=this.interval.getStart();
				b.chromEnd=this.interval.getEnd();
				L = 1+(b.chromEnd-b.chromStart);
				}
			else
				{
				b.chromStart=0;
				b.chromEnd=ci.getSequenceLength();
				L = ci.getSequenceLength();
				}
			b.pix_length = pixel_width_for_chromosomes*(L/this.total_bases);
			b.index=this.boundaries.size();
			this.boundaries.add(b);
			this.tidbound.put(ci.getSequenceIndex(),b);
			this.name2bound.put(ci.getSequenceName(),b);
			x += b.pix_length;
			x += this.distBetweenChrom;
			}
		
		}
	
	
	
	}
