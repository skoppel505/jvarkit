package fr.inserm.umr1087.jvarkit.util.bin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Bin {
private static final int MAX_BIN = 37450;

static private final int binOffsets[] = {512+64+8+1, 64+8+1, 8+1, 1, 0};

static public int bin(int start,int end)
	{
	final int _binFirstShift=17;
	final int _binNextShift=3;
	
	int startBin = start, endBin = end-1, i;
	startBin >>= _binFirstShift;
	endBin >>= _binFirstShift;
	for (i=0; i< binOffsets.length; ++i)
	    {
	    if (startBin == endBin)
	        return binOffsets[i] + startBin;
	    startBin >>= _binNextShift;
	    endBin >>= _binNextShift;
	    }
	throw new IllegalArgumentException("out of range in findBin (max is 512M):"+ start+"-"+end);
	}

static public List<Integer> bins(int chromStart,int chromEnd)
	{
	int k, end = chromEnd;
	if (chromStart >= end) return Collections.emptyList();
	if (end >= 1<<29) end = 1<<29;
	--end;
	List<Integer> L=new ArrayList<Integer>(MAX_BIN);
	L.add(0);
	for (k =    1 + (chromStart>>26); k <=    1 + (end>>26); ++k) L.add(k);
	for (k =    9 + (chromStart>>23); k <=    9 + (end>>23); ++k) L.add(k);
	for (k =   73 + (chromStart>>20); k <=   73 + (end>>20); ++k) L.add(k);
	for (k =  585 + (chromStart>>17); k <=  585 + (end>>17); ++k) L.add(k);
	for (k = 4681 + (chromStart>>14); k <= 4681 + (end>>14); ++k) L.add(k);
	return L;
	}

public static void main(String[] args) {
	System.err.println(bin(10004,10005));
	}

}
