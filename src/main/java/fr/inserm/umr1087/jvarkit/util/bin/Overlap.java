package fr.inserm.umr1087.jvarkit.util.bin;

public enum Overlap
	{
	QUERY_CONTAINS_FEATURE
		{
		@Override
		public boolean match(int qStart, int qEnd, int fStart, int fEnd) {
			return qStart<=fStart && qEnd>=fEnd;
			}
		},
	FEATURE_CONTAINS_QUERY
		{
		@Override
		public boolean match(int qStart, int qEnd, int fStart, int fEnd) {
			return fStart<=qStart && fEnd>=qEnd;
			}
		},
	QUERY_OVERLAP_FEATURE
		{
		@Override
		public boolean match(int qStart, int qEnd, int fStart, int fEnd) {
			return !(fEnd<=qStart || qEnd<=fStart);
			}
		};
	public abstract boolean match(
		int qStart,
		int qEnd,
		int fStart,
		int fEnd
		);
	}
