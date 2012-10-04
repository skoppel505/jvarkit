package fr.inserm.umr1087.jvarkit.util.picard;


import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

public class CigarUtils extends net.sf.picard.util.CigarUtil
	{
	public int getCigarAlignmentLength(SAMRecord rec)
		{
		for(CigarElement e:rec.getCigar().getCigarElements())
			{
			
			}
		return 0;
		}
	public static String toString(CigarElement e)
		{
		return String.valueOf(e.getLength())+e.getOperator();
		}
}
