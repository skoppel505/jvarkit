package fr.inserm.umr1087.jvarkit.util.picard;

import java.util.logging.Logger;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMSequenceRecord;

public class ReferenceBuffer
	{
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");

	public static final int DEFAULT_BUFFER_SIZE=1000000;
	private int buffer_size=DEFAULT_BUFFER_SIZE;
	private SAMSequenceRecord samSeqRec=null;
	private int start0=0;
	private byte buffer[]=null;
	private ReferenceSequenceFile refSeqFile;
	
	public ReferenceBuffer(ReferenceSequenceFile refSeqFile)
		{
		this.refSeqFile=refSeqFile;
		}

	public void setBufferSize(int buffer_size)
		{
		this.buffer_size = Math.max(10000,buffer_size);
		}

	
	public void dispose()
		{
		this.samSeqRec=null;
		this.buffer=null;
		this.start0=0;
		}

	
	public byte getBaseAt(String chromId,int index0)
		{
		if( this.samSeqRec==null ||
			!chromId.equals(this.samSeqRec.getSequenceName()) ||
			buffer==null ||
			index0 < this.start0 ||
			index0 >= (this.start0+buffer.length))
			{
			this.samSeqRec=this.refSeqFile.getSequenceDictionary().getSequence(chromId);
			if(this.samSeqRec==null) throw new IllegalArgumentException("unknown chromosome:"+chromId);
			this.start0 =index0;
			int end=Math.min(this.start0+this.buffer_size,samSeqRec.getSequenceLength());
			LOG.info("refill DNA buffer  "+chromId+":"+start0);
			ReferenceSequence dna= this.refSeqFile.getSubsequenceAt(
						chromId,
						this.start0+1,//1-based
						end//inclusive
						);
			this.buffer=dna.getBases();
			}
		return this.buffer[index0-this.start0];
		}

	}
