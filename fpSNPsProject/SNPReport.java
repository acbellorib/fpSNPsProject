package fps;

import java.io.*;

public class SNPReport
{
	public String transcriptName;
	public int transcriptLength;
	public int snpPos;
	public char refAlleleAtLast;
	public int accumBaseA;
	public int accumBaseC;
	public int accumBaseT;
	public int accumBaseG;
	public int accumUndetermined;
	public int acummRefAllele;
	public int totalDiffAlleles;
	public int numberOfOverlappingReads;
	public int numTotalAvailReads;
	public int numUsedReads;
	public float roundedPercentOfMappedUsedReads;
	public boolean presenceOf3OrMoreAlleles;
	public boolean presenceOfBothAlleles;
	public boolean presenceOnlyAlternateAlleles;
	public boolean presenceOnlyReferenceAllele;
	public boolean presenceOfOddOccurrence;
	public boolean hasBLASTResultsBeforeFilter;
	public boolean hasBLASTResultsAfterFilter;
	public String paralogReport;
	public void outputToFile(BufferedWriter writer) throws IOException
	{
		writer.write(transcriptName+"\t");
		writer.write(transcriptLength+"\t");
		writer.write(snpPos+"\t");
		writer.write(refAlleleAtLast+"\t");
		writer.write(accumBaseA+"\t");
		writer.write(accumBaseC+"\t");
		writer.write(accumBaseT+"\t");
		writer.write(accumBaseG+"\t");
		writer.write(accumUndetermined+"\t");
		writer.write(acummRefAllele+"\t");
		writer.write(totalDiffAlleles+"\t");
		writer.write(numberOfOverlappingReads+"\t");
		writer.write(numTotalAvailReads+"\t");
		writer.write(numUsedReads+"\t");
		writer.write(roundedPercentOfMappedUsedReads+"\t");
		writer.write(presenceOf3OrMoreAlleles+"\t");
		writer.write(presenceOfBothAlleles+"\t");
		writer.write(presenceOnlyAlternateAlleles+"\t");
		writer.write(presenceOnlyReferenceAllele+"\t");
		writer.write(presenceOfOddOccurrence+"\t");
		writer.write(hasBLASTResultsBeforeFilter+"\t");
		writer.write(hasBLASTResultsAfterFilter+"\t");
		writer.write(paralogReport+"\t");
		writer.newLine();
	}
	public static void printHeader(BufferedWriter writer) throws IOException
	{
		writer.write("transcript"+"\t");
		writer.write("length"+"\t");
		writer.write("snpPos"+"\t");
		writer.write("refAllele"+"\t");
		writer.write("Qty 'A'"+"\t");
		writer.write("Qty 'C'"+"\t");
		writer.write("Qty 'T'"+"\t");
		writer.write("Qty 'G'"+"\t");
		writer.write("Qty 'Not Det'"+"\t");
		writer.write("Tot RefAllele"+"\t");
		writer.write("Tot DiffAlleles"+"\t");
		writer.write("Qty OverReads"+"\t");
		writer.write("Qty AvailReads"+"\t");
		writer.write("Qty UsedReads"+"\t");
		writer.write("%MappedReads"+"\t");
		writer.write("3OrMoreAlleles?"+"\t");
		writer.write("BothAlleles?"+"\t");
		writer.write("OnlyAltAlleles?"+"\t");
		writer.write("OnlyRefAllele?"+"\t");
		writer.write("AnyOddOccur?"+"\t");
		writer.write("BLAST hits (BEFORE filtering)?"+"\t");
		writer.write("BLAST hits (AFTER filtering)?"+"\t");
		writer.write("BLAST hit:HasRefAllele?:Qty 'A':Qty 'C':Qty 'T':Qty 'G':Qty 'Not Det';"+"\t");
		writer.newLine();
	}
	
	public void outputToStdout()
	{	
		/*
		//print header
		System.out.print("transcript"+"\t");
		System.out.print("length"+"\t");
		System.out.print("snpPos"+"\t");
		System.out.print("refAllele"+"\t");
		System.out.print("Qty 'A'"+"\t");
		System.out.print("Qty 'C'"+"\t");
		System.out.print("Qty 'T'"+"\t");
		System.out.print("Qty 'G'"+"\t");
		System.out.print("Qty 'Not Det'"+"\t");
		System.out.print("Tot RefAllele"+"\t");
		System.out.print("Tot DiffAlleles"+"\t");
		System.out.print("Qty OverReads"+"\t");
		System.out.print("Qty AvailReads"+"\t");
		System.out.print("Qty UsedReads"+"\t");
		System.out.print("%MappedReads"+"\t");
		System.out.print("3OrMoreAlleles?"+"\t");
		System.out.print("BothAlleles?"+"\t");
		System.out.print("OnlyAltAlleles?"+"\t");
		System.out.print("OnlyRefAllele?"+"\t");
		System.out.print("AnyOddOccur?"+"\t");
		System.out.println();
		
		//print data
		System.out.print(transcriptName+"\t");
		System.out.print(transcriptLength+"\t");
		System.out.print(snpPos+"\t");
		System.out.print(refAlleleAtLast+"\t");
		System.out.print(accumBaseA+"\t");
		System.out.print(accumBaseC+"\t");
		System.out.print(accumBaseT+"\t");
		System.out.print(accumBaseG+"\t");
		System.out.print(accumUndetermined+"\t");
		System.out.print(acummRefAllele+"\t");
		System.out.print(totalDiffAlleles+"\t");
		System.out.print(numberOfOverlappingReads+"\t");
		System.out.print(numTotalAvailReads+"\t");
		System.out.print(numUsedReads+"\t");
		System.out.print(roundedPercentOfMappedUsedReads+"\t");
		System.out.print(presenceOf3OrMoreAlleles+"\t");
		System.out.print(presenceOfBothAlleles+"\t");
		System.out.print(presenceOnlyAlternateAlleles+"\t");
		System.out.print(presenceOnlyReferenceAllele+"\t");
		System.out.print(presenceOfOddOccurrence+"\t");
		System.out.println();
		*/
		System.out.println("\n------------------------------------ START: Reporting block for transcript "+transcriptName+" at position "+snpPos+" -------------------------------------");
		System.out.println("\nTranscript length: "+transcriptLength+" bases.");
		System.out.println("SNP position: "+snpPos+".");
		System.out.println("Reference allele at that position: "+refAlleleAtLast+".");
		System.out.println("Number of reads overlapping the SNP position: "+numberOfOverlappingReads+".");
		System.out.println("Quantity of reads containing the 'A'/'a' allele at that position: "+accumBaseA+".");
		System.out.println("Quantity of reads containing the 'C'/'c' allele at that position: "+accumBaseC+".");
		System.out.println("Quantity of reads containing the 'T'/'t' allele at that position: "+accumBaseT+".");
		System.out.println("Quantity of reads containing the 'G'/'g' allele at that position: "+accumBaseG+".");
		System.out.println("Quantity of reads containing an UNDETERMINED allele at that position: "+accumUndetermined+".");
		System.out.println("Summary of reads containing the REFERENCE ALLELE '"+refAlleleAtLast+"' at that position: "+acummRefAllele+".");
		System.out.println("Summary of reads containing an ALTERNATE ALLELE at that position: "+totalDiffAlleles+".");
		System.out.println("\nTOTAL NUMBER of AVAILABLE READS on transcript assembly time for the mapping related to the "+transcriptName+".sorted.bam file: "+numTotalAvailReads+".");
		System.out.println("ACTUAL NUMBER of MAPPED/USED READS in the mapping related to the "+transcriptName+".sorted.bam file: "+numUsedReads+".");
		System.out.println("PERCENTAGE of MAPPED/USED READS from the initial available set: "+roundedPercentOfMappedUsedReads+"%.");
		System.out.println("\nBoth the REFERENCE ALLELE and some type of ALTERNATE ALLELE are present: "+presenceOfBothAlleles+".");
		System.out.println("Apart from the REFERENCE ALLELE, there is MORE THAN ONE TYPE of ALTERNATE ALLELE present: "+presenceOf3OrMoreAlleles+".");
		System.out.println("ONLY ALTERNATE ALLELE(s) is/are present and there is NO occurrence of the REFERENCE ALLELE: "+presenceOnlyAlternateAlleles+".");
		System.out.println("ONLY THE REFERENCE ALLELE(s) is present: "+presenceOnlyReferenceAllele+".");
		System.out.println("Houston, we have a problem...: "+presenceOfOddOccurrence+".");
		System.out.println("Any kind of BLAST hits found? "+hasBLASTResultsBeforeFilter+".");
		System.out.println("Perfect matches for BLAST hits found? "+hasBLASTResultsAfterFilter+".");
		System.out.println("List of BLAST hits found in the format 'BLAST hit:HasRefAllele?:Qty 'A':Qty 'C':Qty 'T':Qty 'G':Qty 'Not Det';': "+paralogReport+".");
		System.out.println("\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- END: Reporting block for transcript "+transcriptName+" at position "+snpPos+" -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*");
	}
}
