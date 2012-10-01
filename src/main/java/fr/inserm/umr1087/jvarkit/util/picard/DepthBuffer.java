package fr.inserm.umr1087.jvarkit.util.picard;

import java.util.Arrays;
import java.util.Iterator;
import java.util.logging.Logger;

import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.samtools.SAMFileReader;

public class DepthBuffer
	{
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");

	public static final int DEFAULT_BUFFER_SIZE=1000000;
	private int buffer_size=DEFAULT_BUFFER_SIZE;
	private SAMFileReader samReader;
	private String chromName=null;
	private int start0=0;
	private int[] depth=null;
	private int minQual=0;
	public DepthBuffer(SAMFileReader samReader)
		{
		this.samReader=samReader;
		}
	
	public SAMFileReader getSamReader()
		{
		return samReader;
		}
	
	public void setBufferSize(int buffer_size) {
		this.buffer_size = Math.max(10000,buffer_size);
		}
	
	public int getMinQual()
		{
		return minQual;
		}
	
	public void setMinQual(int minQual)
		{
		this.minQual = minQual;
		}
	
	public void dispose()
		{
		this.chromName=null;
		this.depth=null;
		this.start0=0;
		}
	
	public double getMean(String chromId,int start0Incl,int end0Excl)
		{
		int L=end0Excl-start0Incl;
		if(L<=0) return Double.NaN;
		double sum=0;
		while(start0Incl< end0Excl)
			{
			sum+=getDepth(chromId, start0Incl);
			++start0Incl;
			}
		return sum/L;
		}
	
	public int getDepth(String chromId,int position0)
		{
		if(!chromId.equals(this.chromName) ||
			this.depth==null ||
			this.start0> position0 ||
			this.start0+this.depth.length<=position0
			)
			{				
            this.chromName=chromId;
            this.start0=position0;
			this.depth=new int[this.buffer_size];

            
            LOG.info("Fill buffer for "+chromName+":"+start0+"-"+(this.start0+this.depth.length));
            Interval interval=new Interval(
            		chromId,
            		this.start0+1,
            		this.start0+this.depth.length
            		);//Coordinates are 1-based closed ended.
            IntervalList iL=new IntervalList(this.samReader.getFileHeader());
            iL.add(interval);

            SamLocusIterator sli=new SamLocusIterator(this.samReader,iL,true);           
            sli.setMappingQualityScoreCutoff(this.minQual);
                       
             for(Iterator<SamLocusIterator.LocusInfo>  iter=sli.iterator();
	                iter.hasNext();
	                )
                {
                SamLocusIterator.LocusInfo locusInfo=iter.next();
                int pos0= locusInfo.getPosition() - 1;//1-based
                int offset=pos0-this.start0;
                if(offset>=this.depth.length )
			         {
			         System.err.println("???");
			         continue;
			         }
                this.depth[offset]=locusInfo.getRecordAndPositions().size();
                }
			//sli.close();
			}
		return this.depth[position0-this.start0];
		}
	
	}
