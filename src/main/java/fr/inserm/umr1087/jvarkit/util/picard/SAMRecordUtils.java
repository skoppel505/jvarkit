package fr.inserm.umr1087.jvarkit.util.picard;

import net.sf.picard.util.Interval;
import net.sf.samtools.SAMRecord;

public class SAMRecordUtils extends net.sf.samtools.SAMRecordUtil
{

public static Interval toInterval(SAMRecord record)
	{
	return new Interval(record.getReferenceName(),
			record.getAlignmentStart(),
			record.getAlignmentEnd(),
			record.getReadNegativeStrandFlag(),
			record.getReadName()
			);
	}

}
