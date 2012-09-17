package fr.inserm.umr1087.jvarkit.tools.invert20120917;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import net.sf.picard.util.*;
import net.sf.picard.util.SamLocusIterator.LocusInfo;

import net.sf.samtools.*;



public class Invert20120917 {
    private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
    private int minCountReadPerCluster=10;
    private File inputBam;
    private SAMFileReader samFileReader1;
    private SAMFileReader samFileReader2;
    private int mappingQuality=30;
    private int alignId=0;
   

    private void dump(List<SAMRecord> records)
		{
    	int beg1=Integer.MAX_VALUE;
    	int end1=Integer.MIN_VALUE;
    	for(SAMRecord rec2: records)
        	{
    		beg1=Math.min(beg1,rec2.getAlignmentEnd());
    		beg1=Math.min(beg1,rec2.getAlignmentStart());
    		end1=Math.max(end1,rec2.getAlignmentEnd());
    		end1=Math.max(end1,rec2.getAlignmentStart());
        	}
    	int countReads=0;
    	Interval interval=new Interval(records.get(0).getReferenceName(), beg1, end1);
    	LOG.info("interval:"+interval);
        /*IntervalList L=new IntervalList(this.samFileReader2.getFileHeader());
        L.add(interval);
    	SamLocusIterator iter=new SamLocusIterator(this.samFileReader2,L);
    	while(iter.hasNext())
    		{
    		LocusInfo info=iter.next();
    		countReads+=info.getRecordAndPositions().size();
    		}*/
    	++alignId;
    	for(SAMRecord rec2: records)
            {
            System.out.println(
            			
            			rec2.getReferenceName()+"\t"+
            			rec2.getAlignmentStart()+"\t"+
            			rec2.getAlignmentEnd()+"\t"+
            			alignId+"\t"+
            			records.size()+"\t"+
            			rec2.getReadName() +"\t"+
            			this.inputBam//+"\t"
            			//+(countReads/(double)(interval.getEnd()-interval.getStart()))
            			);
            }
		}

    private boolean overlap(SAMRecord r1,SAMRecord r2)
        {
        if(!r1.getReferenceName().equals(r2.getReferenceName())) return false;
        int beg1=Math.min(r1.getAlignmentStart(), r1.getAlignmentEnd());
        int end1=Math.max(r1.getAlignmentStart(), r1.getAlignmentEnd());
        int beg2=Math.min(r2.getAlignmentStart(), r2.getAlignmentEnd());
        int end2=Math.max(r2.getAlignmentStart(), r2.getAlignmentEnd());
        return !(end1<beg2 || end2<beg1);
        }
   
    private void scan()
        {
        LOG.info("reading..."+this.inputBam);
        ArrayList<SAMRecord> records=new ArrayList<SAMRecord>();
        SAMRecordIterator iter=samFileReader1.iterator();
      
        while(iter.hasNext())
            {
            SAMRecord rec=iter.next();
            //if(++nRead%100000==0) LOG.info("reads:"+nRead);
            if(!rec.getReadPairedFlag()) continue;
            if(rec.getProperPairFlag()) continue;
            if(rec.getReadUnmappedFlag()) continue;
            if(rec.getMateUnmappedFlag()) continue;
            if(!rec.getFirstOfPairFlag()) continue;
            if(rec.getReadFailsVendorQualityCheckFlag()) continue;
            if(rec.getDuplicateReadFlag()) continue;
            if(!rec.getMateReferenceName().equals(rec.getReferenceName())) continue;
            if( rec.getReadNegativeStrandFlag() !=
                rec.getMateNegativeStrandFlag())
                {
                continue;
                }
            if(rec.getMappingQuality()<this.mappingQuality) continue;
            if(!records.isEmpty())
                {
                SAMRecord last=records.get(records.size()-1);
                if(!overlap(last,rec))
                    {
                    if(records.size()>=this.minCountReadPerCluster)
                    	{
                    	dump(records);
                    	}
                    records.clear();
                    }
                }
           
            records.add(rec);
            }
        }
    public int run(String[] args)
            throws Exception
        {
    	int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				return 0;
				}
			else if(args[optind].equals("-Q"))
				{
				this.mappingQuality=Integer.parseInt(args[++optind]);
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return 0;
				}
			else
				{
				break;
				}
			++optind;
			}
        //String url= "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA12878/exome_alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.exome.20120522.bam";
        //String url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA12878/exome_alignment/NA12878.mapped.illumina.mosaik.CEU.exome.20110411.bam";
        //InputStream in=new URL(url).openStream();
		while(optind < args.length)
			{
			this.inputBam=new File(args[optind++]);
	        this.samFileReader1=new SAMFileReader(this.inputBam);
	        this.samFileReader2=new SAMFileReader(this.inputBam);
	        this.samFileReader1.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	        this.samFileReader2.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	        scan();
	        this.samFileReader1.close();
	        this.samFileReader2.close();
			}
        return 0;
        }

    /**
     * @param args
     */
    public static void main(String[] args)
        throws Exception
        {
    	Invert20120917 app=new Invert20120917();
        app.run(args);
        }
}
