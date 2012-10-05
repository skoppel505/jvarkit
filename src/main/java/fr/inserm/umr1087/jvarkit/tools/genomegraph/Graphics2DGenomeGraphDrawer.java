package fr.inserm.umr1087.jvarkit.tools.genomegraph;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMSequenceDictionary;




public class Graphics2DGenomeGraphDrawer extends AbstractGenomeGraphDrawer
	{
	private Graphics2D gCurr;
	
	
	public Graphics2DGenomeGraphDrawer(
			SAMSequenceDictionary samSequenceDictionary, Interval interval) {
		super(samSequenceDictionary, interval);
	}

	public Graphics2DGenomeGraphDrawer(
			SAMSequenceDictionary samSequenceDictionary) {
		super(samSequenceDictionary);
	}

	protected void paintData(AbstractGenomeGraphDrawer.Bound bound,GenomeGraph.Data data)
		{
		double x0= bound.convertBaseToPixelX(data.start);
		double x1= bound.convertBaseToPixelX(data.end);
		double xmid=(x0+x1)/2.0;
		double y=convertValueToPixelY(data.value);
		gCurr.fill(new Ellipse2D.Double(xmid-2, y-2, 5, 5));
		//gCurr.fill(new Rectangle2D.Double(x0, y-3, x1==x0?1:x1-x0,6));
		}
	
	protected void paintBound(AbstractGenomeGraphDrawer.Bound bound)
		{
		gCurr.setColor(Color.DARK_GRAY);
		Rectangle2D r=bound.getSize();
		AffineTransform oldtr=this.gCurr.getTransform();
		AffineTransform tr=new AffineTransform(oldtr);
		tr.translate(r.getCenterX(), r.getMinY()-3);
		tr.rotate(-Math.PI/2.0);
		this.gCurr.setTransform(tr);
		gCurr.drawString(bound.samSeqRec.getSequenceName(), 0, 0);
		this.gCurr.setTransform(oldtr);
		
		
		gCurr.setColor(bound.index%2==0?Color.WHITE:Color.LIGHT_GRAY);
		gCurr.fill(r);
		
		gCurr.setColor(Color.BLACK);
		gCurr.draw(new Line2D.Double(
			r.getMinX(),
			r.getMaxY(),
			r.getMaxX(),
			r.getMaxY()
			));
		
	
		}
	
	protected void paintBackground()
		{
		gCurr.setColor(Color.WHITE);
		gCurr.fillRect(0, 0, this.getWidth(), this.getHeight());
		Paint oldPaint=this.gCurr.getPaint();
		
		
		this.gCurr.setPaint(oldPaint);
		}
	
	protected void paintForeground()
		{
		
		Rectangle r=getDrawingBounds();
		

		
		gCurr.setColor(Color.RED);
		gCurr.draw(r);
		
		gCurr.setColor(Color.BLACK);
		gCurr.drawLine(r.x,r.y+r.height,r.x,r.y);
		
		int nStep=10;
		for(int step=0;step<=nStep;++step)
			{
			double v= getMinY()+ step*((getMaxY()-getMinY()))/(double)nStep;
			double y= r.getMaxY()-step*(r.getHeight()/(double)nStep);
			gCurr.draw(new Line2D.Double(r.getX(), y, r.getX()-5,y));
			gCurr.drawString(String.valueOf(v),(int)r.getX()-200, (int)y);
			}
		
		}
	
	public Graphics2D getGraphics()
		{
		return gCurr;
		}
	
	private void _priv_init()
		{
		paintBackground();
		for(AbstractGenomeGraphDrawer.Bound bound:super.boundaries)
			{
			LOG.info(bound.toString());
			paintBound(bound);
			}

		}
	
	public void init(Graphics2D g)
		{
		super.init();
		this.gCurr=g;
		_priv_init();
		}
	
	
	
	public void finish()
		{
		paintForeground();
		}
	
	public BufferedImage createImage(int imageType)
		{
		BufferedImage img= new BufferedImage(this.getWidth(),this.getHeight(),imageType);
		return img;
		}
	
	
	}
