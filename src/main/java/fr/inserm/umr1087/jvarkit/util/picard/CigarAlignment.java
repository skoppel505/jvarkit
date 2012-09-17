package fr.inserm.umr1087.jvarkit.util.picard;

import java.io.File;
import java.util.Iterator;
import java.util.logging.Logger;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class CigarAlignment
	{
	@SuppressWarnings("unused")
	private static final Logger LOG=Logger.getLogger("fr.inserm.umr1087.jvarkit");
	
	
	private int genomePos1;
	private int readPos0;
	private CigarElement cigarElement;
	private SAMRecord record;
	private byte bases[]=null;
	private byte qualities[]=null;
	private int repeat_index;
	
	private CigarAlignment(MyIterator iter)
		{
		this.genomePos1=iter.genomePos1;
		this.readPos0=iter.readPos0;
		this.cigarElement=iter.cigarElement;
		this.record=iter.record;
		this.repeat_index=iter.repeat_index;
		}
	
	public SAMRecord getSAMRecord()
		{
		return record;
		}
	
	public int getReferencePosition1()
		{
		return genomePos1;
		}
	
	public int getReadPosition0()
		{
		return readPos0;
		}
	
	public int getIndexInCigarElement()
		{
		return this.repeat_index;
		}
	
	public CigarElement getCigarElement()
		{
		return cigarElement;
		}
	public CigarOperator getCigarOperator()
		{
		return getCigarElement().getOperator();
		}
	
	public int getCigarCount()
		{
		return getCigarElement().getLength();
		}
	/** base in the read doesn't exists in the REFERENCE. Extra bases in the read */
	public boolean isInsertRef()
		{
		switch(getCigarOperator())
			{
			case I: case S: return true;
			default:return false;
			}
		}
	
	/**  base in the ref doesn't exists in the READ. Extra bases in the REF */
	public boolean isDeletionRef()
		{
		switch(getCigarOperator())
			{
			case D: return true;
			default:return false;
			}
		}
	
	public char getReadBase()
		{
		if(isDeletionRef()) throw new IllegalStateException("can't get base in reference insertion");
		if(this.bases==null)
			{
			bases=getSAMRecord().getReadBases();
			}
		return (char)this.bases[getReadPosition0()];
		}
	
	public int getReadQuality()
		{
		if(isDeletionRef()) throw new IllegalStateException("can't get QUAL in reference insertion");
		if(this.qualities==null)
			{
			qualities=getSAMRecord().getBaseQualities();
			}
		return this.qualities[getReadPosition0()];
		}
	

	public static Iterator<CigarAlignment> iterator(SAMRecord record)
		{
		return new MyIterator(record);
		}

	private static class MyIterator implements Iterator<CigarAlignment>
			{
			private SAMRecord record;
			private int genomePos1;
			private int readPos0;
			private int _next_genomePos1;
			private int _next_readPos0;
		
			private Cigar cigar;
			private int cigar_element_index=-1;
			private int repeat_index=-1;
			private boolean _has_next_called;
			private boolean _has_next;
			private CigarElement cigarElement;
			private CigarOperator cigarOperator;
			
			public MyIterator(SAMRecord record)
				{
				this.record=record;
				if(record.getReadUnmappedFlag()) throw new IllegalArgumentException("unmapped read");
				this.readPos0=0;
				this.genomePos1 = record.getAlignmentStart();
				this._next_genomePos1=this.genomePos1;
				this._next_readPos0=this.readPos0;
				if(this.genomePos1<1) throw new IllegalArgumentException("record");
				this.cigar= record.getCigar();
				}
			
			@Override
			public boolean hasNext()
				{
				if(_has_next_called) return _has_next;
				_has_next_called=true;
				if(this.cigar_element_index==-1)
					{
					this.cigar_element_index=0;
					this.repeat_index=0;
					}
				else
					{
					this.repeat_index++;
					}
				if( this.cigar_element_index < this.cigar.numCigarElements() &&
					this.repeat_index >= this.cigar.getCigarElement(this.cigar_element_index).getLength()
					)
					{
					this.repeat_index=0;
					this.cigar_element_index++;
					}
				//System.err.println("cigarElement:"+this.cigar_element_index+" repeat:"+this.repeat_index);

				if(cigar_element_index>=this.cigar.numCigarElements()) 
					{
					_has_next=false;
					return false;
					}
					
				_has_next=true;
				this.cigarElement = this.cigar.getCigarElement(cigar_element_index);
				this.cigarOperator =  this.cigarElement.getOperator();
				//System.out.println("cigar:"+this.cigarElement.getOperator()+":"+this.cigar_element_index+" "+this.cigarElement.getLength());
				//System.err.println("read:"+this.readPos0+" genom:"+this.genomePos1);
				switch(this.cigarOperator)
					{

					case M:
						{
						this._next_genomePos1 = this.genomePos1+1;
						this._next_readPos0 = this.readPos0+1;
						break;
						}
					//insertion from the reference
					case I:
					case S:
						{
						//this._next_genomePos1 = this.genomePos1+1;
						this._next_readPos0 = this.readPos0+1;
						break;	
						}
					//deletion from the reference
					case D:
						{
						this._next_genomePos1 = this.genomePos1+1;
						//this._next_readPos1 = this.readPos1+1;
						break;	
						}
					case P:
					case H:
					case N:
						{
							
						break;
						}
			
					default: throw new IllegalStateException("Cannot handle operator:"+this.cigarOperator);
					}
				
				return true;
				}
			@Override
			public CigarAlignment next()
				{
				if(!_has_next_called) hasNext();
				if(!_has_next) throw new IllegalStateException();
				_has_next_called=false;
				_has_next=false;
				CigarAlignment next=new CigarAlignment(this);
				this.genomePos1 = this._next_genomePos1 ;
				this.readPos0 = this._next_readPos0 ;
				
				return next;
				}
			@Override
			public void remove() {
				throw new UnsupportedOperationException();
				}
			}
	
	public static void main(String[] args)
		{
		args=new String[]{"/commun/data/users/cfaucheron/aln_20120329/S0529/data_S0529/S0529_sort.nodup.bam"};
		ReferenceSequenceFile rsf=ReferenceSequenceFileFactory.getReferenceSequenceFile(new File("/commun/data/pubdb/ucsc/hg19/chromosomes/hg19.fa"));
		int count=0;
		for(String filename:args)
			{
			File file=new File(filename);
			SAMFileReader samIn=new SAMFileReader(file);
			SAMRecordIterator r=samIn.iterator();
			while(r.hasNext())
				{

				SAMRecord rec=r.next();
				if(rec.getReadUnmappedFlag()) continue;

				if(++count>10000) break;
				
				if(rec.getAlignmentStart()> rec.getAlignmentEnd()) throw new IllegalStateException();
				byte bases[]=rsf.getSubsequenceAt(
						rec.getReferenceName(),
						rec.getAlignmentStart(),
						Math.max(rec.getAlignmentEnd(),
								rec.getAlignmentStart()+rec.getCigar().getPaddedReferenceLength()
								)
						).getBases();
				Iterator<CigarAlignment> i=CigarAlignment.iterator(rec);
				/*System.err.println(rec.getCigarString());
				System.err.println(bases.length);
				System.err.println("start:"+rec.getAlignmentStart());*/
				StringBuilder s1=new StringBuilder();
				StringBuilder s2=new StringBuilder();
				
				while(i.hasNext())
					{
					CigarAlignment caln=i.next();
					/*
					System.err.println(rec.getCigarString());
					
					System.err.println("bases.length:"+bases.length);
					System.err.println("refpos:"+caln.getReferencePosition1());
					System.err.println("readpos:"+rec.getAlignmentStart());
					*/
					if(caln.getReferencePosition1()-rec.getAlignmentStart()>=bases.length)
						{
						System.out.println("SHORT!");
						System.out.println("op:"+caln.getCigarOperator());
						System.out.println("read start:"+rec.getAlignmentStart());
						System.out.println("clan.pos1:"+caln.getReferencePosition1());
						System.out.println("read end:"+rec.getAlignmentEnd());
						System.out.println("bases.length:"+bases.length);
						System.out.println("getPaddedReferenceLength:"+rec.getCigar().getPaddedReferenceLength());
						System.out.println("getReferenceLength:"+rec.getCigar().getReferenceLength());
						System.out.println("getReadLength:"+rec.getCigar().getReadLength());
						System.out.println("cigar.read.length:"+Cigar.getReadLength(rec.getCigar().getCigarElements()));
						count=2000;
						break;
						}
					if(caln.isInsertRef())
						{
						s2.append("-");
						s1.append(caln.getReadBase());
						}
					else if(caln.isDeletionRef())
						{	
						s2.append((char)bases[caln.getReferencePosition1()-rec.getAlignmentStart()]);
						s1.append("-");
						}
					else
						{
						s2.append((char)bases[caln.getReferencePosition1()-rec.getAlignmentStart()]);
						s1.append(caln.getReadBase());
						}
					//System.out.println(s1);
					//System.out.println(s2);
					//System.out.println();
					}
				System.out.println(rec.getCigarString()+" "+rec.getReferenceName()+":"+rec.getAlignmentStart());
				System.out.println("ref :"+new String(bases));
				System.out.println("read:"+new String(rec.getReadBases()));
				System.out.println();
				System.out.println(s1);
				System.out.println(s2);
				System.out.println();
				}
			samIn.close();
			}
		}
	}
