package fps;

import java.io.*;
//import java.nio.file.*;
import java.util.*;
import utils.fasta.*;
import utils.validation.transcripts.*;
import net.sf.picard.reference.*;
//import net.sf.samtools.*;
//import java.math.*;

public class AnalyzeHomozygousSNPs
{
	static final String USAGE = "java fps.AnalyzeHomozygousSNPs <FASTA file> <Truncate name at space character? true | false> <List of Homozygous SNPs file> <Absolute path of Trinity .reads directory> <Number of mismatches for BWA alignment> <Number of BWA running threads> <outputFile> <Run BLAST against any specific database? true | false> <Target BLAST database file | none>";
	public static int totalUsedReads = 0;
	public static HashMap<String, ArrayList<Integer>> usedReadsLookup = new HashMap<String, ArrayList<Integer>>();
	//this keeps a tab of all our SNPReport objects
	static ArrayList<SNPReport> snpReports = new ArrayList<SNPReport>();
	public static HashMap<Integer, Float> bitScoreLookup = new HashMap<Integer, Float>();
	public static String databaseFile = null;
	
	static boolean debug = true;
	//static boolean debug = false;
	
	public static void main(String[] args)
	{
		if (args.length != 9) {
			System.out.println(USAGE);
			System.exit(0);
		}
		try
		{
			boolean runBLAST = Boolean.parseBoolean(args[7]);
			if (runBLAST == true)
			{
				databaseFile = new String(args[8]);
			}
			File fastaFile = new File(args[0]);
			boolean truncateNameAtSpaceChar = Boolean.parseBoolean(args[1]);
			File listHomSNPsFile = new File(args[2]);
			String path2TrinityDir = new String(args[3]);
			int numberOfMismatches4BWA = Integer.parseInt(args[4]);
			int numberOfThreads4BWA = Integer.parseInt(args[5]);
			analyze(fastaFile, truncateNameAtSpaceChar, listHomSNPsFile, path2TrinityDir, numberOfMismatches4BWA, numberOfThreads4BWA, databaseFile);
			outputResults(new File(args[6]));
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static void outputResults(File outputFile) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		//print a header 
		SNPReport.printHeader(writer);
		//for each SNP, print a row with the SNP report
		for(SNPReport snpReport : snpReports)
			snpReport.outputToFile(writer);
		//don't forget to close the output stream
		writer.close();
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static void analyze(File fastaFile, boolean truncateNameAtSpaceChar, File listHomSNPsFile, String path2TrinityDir, int numberOfMismatches4BWA, int numberOfThreads4BWA, String databaseFile) throws IOException
	{
		System.out.println("\nParsing a given transcripts assembly file and storing each entry...");
		// read all transcripts into sequence lookup (key = transcript name, value = sequence)
		LinkedHashMap<String, String> sequenceLookup = FastaParser.parseFile(fastaFile, truncateNameAtSpaceChar);
		// creates the list of entries of the kind comp33_c0_seq1_266 (transcript and SNP position) from a given list of Homozygous SNP file
		ArrayList<String> snpList = parseSNPList(listHomSNPsFile);
		// make lookup of transcripts (keys) versus lists of SNPs within these (values)
		HashMap<String, TreeSet<Integer>> snpLookup = makeSNPLookup(snpList);
		// for each transcript in keyset of lookup
		for (String transcName : snpLookup.keySet())
		{
			// look up its sequence in the sequence lookup
			String sequence = sequenceLookup.get(transcName);
			// output it in FASTA format
			File transcFastaFile = outputSeq2Fasta(transcName, sequence);
			// index this file for BWA and by SAMTOOLS
			if(debug) System.out.println("\nProcessing the file: "+transcFastaFile.getAbsolutePath());
			String com4BWAIndex = "bwa index "+transcFastaFile.getAbsolutePath();
			RunJob.runJob(com4BWAIndex);
			// index this file by SAMTOOLS
			String com4Faidx = "samtools faidx "+transcFastaFile.getAbsolutePath();
			RunJob.runJob(com4Faidx);
			// resultFaidxFile(com4Faidx, transcName);
			// get the proper reads file from Trinity output and converting it to an equivalent FASTA file of the reads only
			if (debug) System.out.println("\n\nTrying to retrieve the used reads in the original transcript assembly...");
			// find the proper directory for the specific transcript compXXXX in the trinity_out_dir path
			String shortTranscName = transcName.substring(0, transcName.indexOf("_"));
			// Line from when I was using the find command approach: String comFindRawCompReadsFile = "find "+path2TrinityDir+" -name "+shortTranscName+".reads";
			// Line from when I was using the find command approach: System.out.println("\nExecuting find command...\n");
			if (debug) System.out.println("\nTrying to fetch the "+path2TrinityDir+shortTranscName+".reads file");
			// Line from when I was using the find command approach: System.out.println(comFindRawCompReadsFile);
			// Line from when I was using the find command approach: String comFindRawCompReadsFile = path2TrinityDir+shortTranscName;
			// Line from when I was using the find command approach: String resultFindRawCompReadsFile = tailoredFind(comFindRawCompReadsFile);
			String resultFindRawCompReadsFile = path2TrinityDir+shortTranscName+".reads";
			// Line from when I was using the find command approach: String resultFindRawCompReadsFile = comFindRawCompReadsFile.
			File handledRawCompReadsFile = trinityCompReadsFile2Fasta(resultFindRawCompReadsFile, shortTranscName);
			// map all reads with BWA to the transcript itself ... All this block will be refactored later in order to call a specific helper method for each task.........
			if (debug) System.out.println("\nExecuting the mappings of each set of reads related to the respective transcript assembly event...");
			// aligning the reads
			System.out.println("\nGenerating the required file: "+transcName+".sai");
			String com4BWAAlignment = "bwa aln -n "+numberOfMismatches4BWA+" -l 100 -t "+numberOfThreads4BWA+" -o 0 -f "+transcName+".sai "+transcFastaFile.getAbsolutePath() +" "+handledRawCompReadsFile.getAbsolutePath();
			RunJob.runJob(com4BWAAlignment);
			// converting SAM2BAM
			if (debug) System.out.println("\nCreating the required file: "+transcName+".bam");
			File resultFromBWAAlignment = new File(transcName+".sai");
			String com4BWASAM = "bwa samse -f "+transcName+".sam "+transcFastaFile.getAbsolutePath()+" "+resultFromBWAAlignment.getAbsolutePath()+" "+handledRawCompReadsFile.getAbsolutePath();
			RunJob.runJob(com4BWASAM);
			File bwaSAMFile = new File(transcName+".sam");
			String comSAM2BAM = "samtools view -S -b -o "+transcName+".bam "+bwaSAMFile.getAbsolutePath();
			RunJob.runJob(comSAM2BAM);
			// sorting BAM
			File convBAMFile = new File(transcName+".bam");
			if (debug) System.out.println("\nCreating the required file: "+transcName+".sorted.bam");
			String comSortBAM = "samtools sort "+convBAMFile.getAbsolutePath()+" "+transcName+".sorted";
			RunJob.runJob(comSortBAM);
			// indexing BAM
			if (debug) System.out.println("\nCreating the required file: "+transcName+".sorted.bam.bai");
			File sortedBAMFile = new File(transcName+".sorted.bam");
			String comIndexBAM = "samtools index "+sortedBAMFile.getAbsolutePath();
			RunJob.runJob(comIndexBAM);
			new File(transcName+".sorted.bam.bai");
			//System.out.println("\n -------- CHECKPOINT -------- : Confirming that the totalUsedReads variable has value = 0 before being used for another transcript event: "+ totalUsedReads);
			//System.out.println("\n -------- CHECKPOINT -------- : Confirming the values of the usedReadsLookup HASHMAP....");
			//System.out.println("\n -------- CHECKPOINT -------- : Getting the values for the given transcName "+transcName+" in usedReadsLookup HASHMAP: "+usedReadsLookup.get(transcName));
			
			// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			// Cleaning up stage (to be coded after the proper inventory of all created files....)
			// transcFastaFile.delete();
			// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			
			// for each SNP in the transcript
			for (Integer snpPos : snpLookup.get(transcName))
			{
				//a new report object for this SNP only
				if (debug) System.out.println("\nEntering at the loop step for "+snpPos+".");
				SNPReport snpReport = new SNPReport();
				snpReports.add(snpReport);
				// run samtools flagstat to get the percentage/count of reads mapped
				if (debug) System.out.println("\nGetting the proper statistics related to the file: "+transcName+".sorted.bam");
				String comFlagstat = "samtools flagstat "+sortedBAMFile.getAbsolutePath()+" > "+transcName+".flagstat.txt";
				resultFlagstatFile(comFlagstat, transcName, snpReport);
				// get the reference allele at the SNP position
				int finalPosOfTranscript = sequence.length();
				//System.out.println("\nThe transcript "+transcName+" has a length of "+finalPosOfTranscript+" bases.");
				char refAlleleAtPos = sequence.charAt((snpPos-1));
				// ------------------------------------------------------------------------------- NOT USED ANYMORE BLOCK - kept here for further documentation -----------------------------------------------------------------------------------------------------------------------------
				//System.out.println("\nThe found reference allele for the transcript "+transcName+" at position "+position+" is: "+refAlleleAtPos+".");
				// get all reads that intersect with the SNP position (it will be handled by AlleleExtractor class)
				// extract all the alleles at the SNP position (it will be handled by AlleleExtractor class)
				// count how many alternate alleles we have (it will be handled by AlleleExtractor class)
				// for each alternate allele (it will be handled by AlleleExtractor class)
				// record its count (it will be handled by AlleleExtractor class)
				// ------------------------------------------------------------------------------- END OF NOT USED ANYMORE BLOCK ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				IndexedFastaSequenceFile faidx = new IndexedFastaSequenceFile(new File(transcFastaFile.getAbsolutePath()));
				//extract allele counts
				AlleleExtractor.extractAlleles(new File(sortedBAMFile.getAbsolutePath()), snpPos.intValue(), faidx, transcName, refAlleleAtPos,finalPosOfTranscript, snpReport, databaseFile);
				// report/output this appropriately
				if(debug)
					snpReport.outputToStdout();
			}
		}
		System.out.println("\ndone - Processed and stored each entry of the given transcripts assembly file...");
		System.out.println("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =");
		/*System.out.println("\nNow entering in the stage 2 of the analysis...");
		System.out.println("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =");
		System.out.println();
		System.out.println("\nDoing lots of BLAST operations to compute the best alignments and associated bit scores from all the previous BLAST operations executed for each set of reads for each transcript...");
		String com4OverallBLAST = "/home/ar41690/phd_libraries/scan.sh _OverlappingReads.fasta";
		RunJob.runJob(com4OverallBLAST);
		String com4appendingBLASTResults = "/home/ar41690/phd_libraries/appendingBlasts.sh _OverReads_vs_FLcDNA.txt";
		RunJob.runJob(com4appendingBLASTResults);
		System.out.println("\ndone - Processed and BLASTed all the necessary files for the best alignments calculation stage...");
		System.out.println("\nNow populating the best 'bit scores' lookup table...");
		File resultOfAllBlasts = new File("resultOfAllBLASTs.txt");
		// make lookup of whole lenghts of reads (keys) versus better BLAST bit scores obtained (values)
		makeBitScoreLookup(resultOfAllBlasts);*/
		System.out.println("\nGood job!!! Tool execution FINISHED!");
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static ArrayList<String> parseSNPList(File listFile)
	{
		System.out.println("\nReading a given Homs. SNP file and populating a list of its entries...");
		// reads the Homs. SNPs file to create an array list (of Strings) for each entry of the file
		BufferedReader br = null;
		ArrayList<String> listTest = new ArrayList<String>();
		try
		{
			// reads each line of the file and populate the array list named listTest
			String line;
			br = new BufferedReader(new FileReader(listFile));
			while ((line = br.readLine()) != null)
			{
				listTest.add(line.trim());
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		finally
		{
			try
			{
				if (br != null)
					br.close();
			}
			catch (IOException ex)
			{
				ex.printStackTrace();
			}
		}
		System.out.println("\ndone -- The processed file had "+listTest.size()+" lines.");
		return listTest;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	// make lookup of transcripts (keys) versus lists of SNPs within these (values) from a given snpList produced by the parseSNPList method
	private static HashMap<String, TreeSet<Integer>> makeSNPLookup(ArrayList<String> snpList)
	{
		System.out.println("\nMaking SNP Lookup 'data structure'...");
		// Stores transcripts names as keys and lists of SNPs positions as values
		HashMap<String, TreeSet<Integer>> snpLookup = new HashMap<String, TreeSet<Integer>>();
		// comp33_c0_seq1_266
		// or
		// comp1308_c0_seq1_64
		// comp1308_c0_seq1_110
		// iterates over list of SNPs
		for (String snpName : snpList)
		{
			// extracting transcript names and SNPs positions list
			String transcName = snpName.substring(0, snpName.lastIndexOf("_"));
			String position = snpName.substring(snpName.lastIndexOf("_")+1);
			// get the list of SNP positions for this transcript and check if it exists
			TreeSet<Integer> snpPositions = snpLookup.get(transcName);
			// if the transcript doesn't exist, creates a new entry for this transcript and its snpPositions list
			if (snpPositions == null)
			{
				// new instance of snpPositions object
				snpPositions = new TreeSet<Integer>();
				// adding the list of transcript and respective snpPosition
				snpLookup.put(transcName, snpPositions);
			}
			// adding the position to list of snpPositions for this transcript
			snpPositions.add(Integer.parseInt(position));
		}
		System.out.println("\ndone -- Entries processed: "+snpLookup.size()+".");
		return snpLookup;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static File outputSeq2Fasta(String sequenceName, String sequence) throws IOException
	{
		// naming the file with the sequence name and adding the .fasta extension
		File fastaFile = new File(sequenceName+".fasta");
		// opening an output Stream to that file
		BufferedWriter bw = new BufferedWriter(new FileWriter(fastaFile));
		// writing sequence name FASTA style
		bw.write(">"+sequenceName);
		// newline method putting the plataform specific new line char
		bw.newLine();
		// writes the sequence itself
		bw.write(sequence);
		// newline method putting the plataform specific new line char
		bw.newLine();
		// closing the output Stream
		bw.close();
		return fastaFile;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static File trinityCompReadsFile2Fasta(String workdir, String shortTranscName) throws IOException
	{
		// Manipulation of the proper compXXXX.reads file in the Trinity chrysalis/RawComps.YYYY/ directory structure
		// naming the file with the sequence name and adding the .fasta extension
		File fastaFile = new File(shortTranscName+"_UsedReads.fasta");
		// opening an output Stream to that file
		BufferedWriter bw = new BufferedWriter(new FileWriter(fastaFile));
		Integer readCount = 0;
		// System.out.println("\nReading the RawComp reads file for a given transcript...");
		try
		{
			// reading the RawComp reads file
			// System.out.println("\nTrying to open the file...");
			FileInputStream fstream = new FileInputStream(workdir);
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String line = null;
			while ((line = br.readLine()) != null)
			{
				// System.out.println(line);
				if (!line.startsWith("Comp"))
				{
					String[] lineArray = line.split("\t");
					String readAccession = lineArray[0];
					// System.out.println("\nPosition at lineArray[0] which would represent readAccession: "+readAccession);
					String seq = lineArray[6];
					// System.out.println("\nPosition at lineArray[6] which would represent seq: "+seq);
					// writing sequence name FASTA style
					bw.write(readAccession);
					// newline method putting the plataform specific new line char
					bw.newLine();
					readCount++;
					// writes the sequence itself
					bw.write(seq);
					// newline method putting the plataform specific new line char
					bw.newLine();
				}
			}
			in.close();
			bw.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		// finally
		// {
		// try
		// {
		// if (br != null)
		// br.close();
		// }
		// catch (Exception ex)
		// {
		// ex.printStackTrace();
		// }
		// }
		System.out.println("\ndone -- File "+shortTranscName+"_UsedReads.fasta created with "+readCount+" reads.");
		totalUsedReads = readCount;
		//PopulateTotalUsedReads(totalUsedReads);
		readCount = 0;
		return fastaFile;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static String tailoredFind(String input4FindCommand) throws IOException
	{
		Runtime rt = Runtime.getRuntime();
		Process pr = rt.exec(input4FindCommand);
		BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));
		String line = null;
		String workdir = null;
		while ((line = input.readLine()) != null)
		{
			System.out.println("\nFound and processing the file: "+line+".");
			workdir = line;
		}
		return workdir;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static void resultFlagstatFile(String input4FlagstatCommand, String transcName, SNPReport snpReport) throws IOException
	{
		File flagstatFile = new File(transcName+".flagstat.txt");
		// opening an output Stream to that file
		BufferedWriter bw = new BufferedWriter(new FileWriter(flagstatFile));
		Runtime rt = Runtime.getRuntime();
		Process pr = rt.exec(input4FlagstatCommand);
		BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));
		String str;
		int numLines = 0;
		int numTotalAvailReads = 0;
		int numUsedReads = 0;
		while ((str = input.readLine()) != null)
		{
			try
			{
				bw.write(str);
				bw.newLine();
				numLines++;
				if (numLines == 1)
				{
					//System.out.println("\nRetrieving the TOTAL NUMBER OF AVAILABLE READS for the mapping related to the "+transcName+".sorted.bam file...");
					//System.out.println("\n"+str);
					numTotalAvailReads = Integer.parseInt(str.substring(0, str.indexOf(" ")));
					//System.out.println("\n"+numTotalAvailReads);
				}
				if (numLines == 3)
				{
					//System.out.println("\nRetrieving the NUMBER OF MAPPED/USED READS in the mapping related to the "+transcName+".sorted.bam file...");
					//System.out.println("\n"+str+"\n");
					numUsedReads = Integer.parseInt(str.substring(0, str.indexOf(" ")));
					//System.out.println("\n"+numUsedReads);
				}
				// populate lookup of transcripts (keys) versus lists of TOTAL AVAILABLE READS IN THE MAPPING and TOTAL USED READS IN THE MAPPING (values)
				float percentageOfMappedUsedReads = ((numUsedReads/(float)numTotalAvailReads)*100);
				int roundedPercentOfMappedUsedReads = Math.round(percentageOfMappedUsedReads);
				populateAvailVsUsedReadsLookup(numTotalAvailReads, numUsedReads, roundedPercentOfMappedUsedReads, snpReport);
			}
			catch (IOException ex)
			{
			}
		}
		bw.close();
		numLines = 0;
		numTotalAvailReads = 0;
		numUsedReads = 0;
		totalUsedReads = 0;
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	public static void populateAvailVsUsedReadsLookup(int numTotalAvailReads, int numUsedReads, float roundedPercentOfMappedUsedReads, SNPReport snpReport)
	{
		snpReport.numTotalAvailReads = numTotalAvailReads;
		snpReport.numUsedReads = numUsedReads;
		snpReport.roundedPercentOfMappedUsedReads = roundedPercentOfMappedUsedReads;
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	private static void makeBitScoreLookup(File resultOfAllBlasts)
	{
		System.out.println("\nReading a given resultOfAllBlasts.txt file and populating a list of its unique entries...");
		// Manipulation of the proper resultOfAllBlasts.txt file to create a lookup table of lenghts (of reads) and perfect alignment bit scores amongst the lines of the file
		Integer lineCount = 0;
		try
		{
			FileInputStream fstream = new FileInputStream(resultOfAllBlasts);
			//HashMap<Integer, Float> bitScoreAndLengthLookup = new HashMap<Integer, Float>();
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String line = null;
			while ((line = br.readLine()) != null)
			{
				// System.out.println(line);
				lineCount++;
				String[] lineArray = line.split("\t");
				Integer alignmentLength = Integer.parseInt(lineArray[2]);
				Integer startOfAlignment = Integer.parseInt(lineArray[3]);
				Integer endOfAlignment = Integer.parseInt(lineArray[4]);
				float bitScore = Float.parseFloat(lineArray[5]);
				// System.out.println("\nPosition at lineArray[2] which would represent alignmentLength: "+alignmentLength);
				// System.out.println("\nPosition at lineArray[3] which would represent startOfAlignment: "+startOfAlignment);
				// System.out.println("\nPosition at lineArray[4] which would represent endOfAlignment: "+endOfAlignment);
				// System.out.println("\nPosition at lineArray[5] which would represent bitScore: "+bitScore);
				if ((startOfAlignment == 1) && (endOfAlignment == alignmentLength))
				{
					bitScoreLookup.put(alignmentLength, bitScore);
				}
			}
			in.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		System.out.println("\ndone -- "+lineCount+" lines of file "+resultOfAllBlasts+" processed and bitScoreLookup table created with "+bitScoreLookup.size()+" rows.");
		lineCount = 0;
		if(debug)
		{
		        System.out.println("\nbitScoreLookup entries are...");
		         for (Integer alignmentLength : bitScoreLookup.keySet())
		         {
		                String key = alignmentLength.toString();
		                String value = bitScoreLookup.get(alignmentLength).toString();
		                System.out.println(key+" "+value);  
		         } 
		}
	}
}
