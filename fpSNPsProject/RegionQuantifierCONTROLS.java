package fps;

import java.io.*;
import java.util.*;
//import utils.validation.transcripts.new RunJob_Antonio();
import fps.RunJob_Antonio;
import net.sf.picard.reference.*;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import utils.entities.SNPSequence;
import utils.fasta.FastaParser;

public class RegionQuantifierCONTROLS
{
	private static BufferedWriter writer;
	private static Timer timer;
	static final String USAGE = "java fps.RegionQuantifierCONTROLS <List of SNPs file> <FASTA file> <Truncate name at space character? true | false> <BAM file> <Read length of the simulated dataset> <outputFile> <Run BLAST against any specific database? true | false> <Target BLAST database file | none>";
	public static int totalUsedReads = 0;
	public static HashMap<String, ArrayList<Integer>> usedReadsLookup = new HashMap<String, ArrayList<Integer>>();
	// this keeps a tab of all our SNPReport objects
	static ArrayList<SNPReport> snpReports = new ArrayList<SNPReport>();
	// this keeps a tab of all our SNPReport2 objects (an alternative reporting structure in case we use the database file)
	static ArrayList<SNPReport2> snpReports2 = new ArrayList<SNPReport2>();
	// lookup of lookups where at the top level we are storing against the SNP name a lookup of (sampleName, SNPresults)
	HashMap<String, SNPSequence> snpSeqLookup;
	// structure which will be further used to collect the quantities of different alleles found at the SNP position
	public static HashMap<String, ArrayList<Integer>> alleleQuantitiesLookup = new HashMap<String, ArrayList<Integer>>();
	// structure which will keep a lookup of the "transcriptName + snpPos" (key) while storing the HashMap of the overlapping read and its respective allele at the snpPos
	public static HashMap<String, HashMap<String, Character>> infoAboutReadAndAllele = new HashMap<String, HashMap<String, Character>>();
	public static String databaseFile = null;
	public static HashMap<String, String> chrmLookup = new HashMap<String, String>();
	public static int accumMisassembly = 0;
	public static int accumNoMisassembly = 0;
	public static int accumIsRevComp = 0;
	public static int accumIsNotRevComp = 0;
	public static int accumNotPossible2DetermineSNPPosOnGenome = 0;
	public static int accumNonACTGEvents = 0;
	public static int accumNormalACTGEvents = 0;
	public static int accumSNPSiteVisitedByBlast = 0;
	public static int accumSNPSiteNotVisitedByBlast = 0;
	public static float accumAveragePercentageMismappedReadsPerSNPEvent = 0;
	public static int numberOfSNPEvents = 0;
	public static int accumPresenceOfBothAlleles = 0;
	public static int accumPresenceOnlyAlternateAlleles = 0;
	public static int accumEnteringAlleleExtractor = 0;
	// Created just for the CONTROLS scenario...
	static ArrayList<String> controlContigs = new ArrayList<String>();			
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// debugger switch block
	static boolean debug = true;
	// static boolean debug = false;
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	public static void main(String[] args)
	{
		if (args.length != 8)
		{
			System.out.println(USAGE);
			System.exit(0);
		}
		try
		{
			// New block to try to get the number of running processes for debugging purposes... -----------------------------------------------------------------------
			File logFile = new File(args[5] + ".log");
			writer = new BufferedWriter(new FileWriter(logFile));
			TimerTask run = new TimerTask()
			{
				@Override
				public void run()
				{
					try
					{
						String cmd = "../testScript.sh";
						//new RunJob_Antonio().runJob(cmd);
						String [] runJobClassOutputs = new RunJob_Antonio().runJob(cmd);
						String runJobClassStdOut = runJobClassOutputs[0];
						String runJobClassStdErr = runJobClassOutputs[1];
						System.out.println("runJobClassStdOut = "+runJobClassStdOut);
						System.out.println("runJobClassStdErr = "+runJobClassStdErr);
						File lsofLogFile = new File("lsof.log");
						String line;
						BufferedReader br = new BufferedReader(new FileReader(lsofLogFile));
						while ((line = br.readLine()) != null)
						{
							//String[] lineArray;
							//String numberOfRunningProcesses = lineArray[0];
							System.out.println();
							writer.write("\nNumber of opened files in the system: "+line+" / Current iteration of the AlleleExtractor method: "+accumEnteringAlleleExtractor);
						}
						writer.flush();
						br.close();
					}
					catch (IOException e)
					{
						e.printStackTrace();
					}
				}
			};
			timer = new Timer();
			timer.schedule(run, 60000, 60000);
			//----------------------------------------------------------------------------------------------------------------------------------------------------------
			boolean runBLAST = Boolean.parseBoolean(args[6]);
			if (runBLAST == true)
			{
				databaseFile = new String(args[7]);
			}
			File listOfSNPsFile = new File(args[0]);
			File fastaFile = new File(args[1]);
			boolean truncateNameAtSpaceChar = Boolean.parseBoolean(args[2]);
			File sortedBAMFile = new File(args[3]);
			// File samFile = new File(args[3]);
			int readLengthSimulation = Integer.parseInt(args[4]);
			quantify(listOfSNPsFile, fastaFile, truncateNameAtSpaceChar, sortedBAMFile, readLengthSimulation);
			outputResults(new File(args[5]));
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	private static void outputResults(File outputFile) throws IOException
	{
		if (databaseFile == null)
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
			// print a header
			SNPReport.printHeader(writer);
			// for each SNP, print a row with the SNP report
			for (SNPReport snpReport : snpReports)
				snpReport.outputToFile(writer);
			// don't forget to close the output stream
			writer.close();
		}
		else
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
			// print a header
			SNPReport2.printHeader(writer);
			// for each SNP, print a row with the SNP report
			for (SNPReport2 snpReport2 : snpReports2)
				snpReport2.outputToFile(writer);
			// don't forget to close the output stream
			writer.close();
		}
		timer.cancel();
		writer.close();
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// private static void quantify(File listOfSNPsFile, File fastaFile, boolean truncateNameAtSpaceChar, File sortedBAMFile, int readLengthSimulation) throws IOException
	private static void quantify(File listOfSNPsFile, File fastaFile, boolean truncateNameAtSpaceChar, File sortedBAMFile, int readLengthSimulation) throws IOException
	{
		System.out.println("Tool execution is starting... Analysing the files related to the " + readLengthSimulation + "bp simulated reads dataset... MAY THE FORCE BE WITH US!!!");
		// creates the list of entries of the kind NODE_10001_length_277_cov_91.635376_143 (contig name and SNP position) from a given text file in the format
		// 50bp_TAIR10_aligned_cf-0.12_sedded.sorted_Mcomb_mergedSNPs.txt
		if (databaseFile != null)
		{
			System.out.println("\nSince the BLAST option WAS CHOSEN, creating a data structure to store each CHROMOSOME and its respective SEQUENCE related to the given BLAST dbFile " + databaseFile + ".");
			File dbFile = new File(databaseFile);
			truncateNameAtSpaceChar = true;
			if (debug)
				System.out.println("\nParsing a given chromosomes FASTA file and storing each entry and its respective sequence for further usage...");
			// read all chromosomes into a sequence lookup (key = chromosome name, value = sequence)
			chrmLookup = FastaParser.parseFile(dbFile, truncateNameAtSpaceChar);
		}
		else
		{
			chrmLookup = null;
			System.out.println("\nSince the BLAST option WAS NOT CHOSEN, no dbFile is present to base the creation of a data structure to store each CHROMOSOME and its respective SEQUENCE for further usage...");
		}
		System.out.println("\nCreating data structures comprising contigs and their respective (single or multiple) SNP positions...");
		ArrayList<String> snpList = parseSNPList(listOfSNPsFile);
		// ----------------------------------------- New block of code created specifically for the "RegionQuantifier" class
		// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		numberOfSNPEvents = snpList.size();
		// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (snpList.isEmpty())
		{
			System.out.println("\nSorry, some of the reference sequence SNP positions in the " + listOfSNPsFile.getName() + " file are NOT CONSISTENT with the position denoted in the contig name for that entry. PLEASE, RECHECK THIS!");
			System.out.println("\nProgram is now aborting...");
			System.exit(1);
		}
		// make lookup of contig names (keys) versus lists of SNPs within these (values)
		HashMap<String, TreeSet<Integer>> snpLookup = makeSNPLookup(snpList);
		// read all contigs into a sequence lookup (key = contig name, value = sequence)
		System.out.println("\nParsing a given assembly file and storing each entry and its respective sequence...");
		LinkedHashMap<String, String> sequenceLookup = FastaParser.parseFile(fastaFile, truncateNameAtSpaceChar);
		// for each contig in keyset of lookup
		for (String contigName : snpLookup.keySet())
		{
			if (debug) System.out.println("\nEntering in the OUTER LOOP step for contig " + contigName + ".");
			// look up its sequence in the sequence lookup
			String sequence = sequenceLookup.get(contigName);
			// output it in FASTA format
			File contigFastaFile = outputSeq2Fasta(contigName, sequence);
			String com4Faidx = "samtools faidx " + contigFastaFile.getAbsolutePath();
			//new RunJob_Antonio().runJob(com4Faidx);
			String [] runJobClassOutputs = new RunJob_Antonio().runJob(com4Faidx);
			String runJobClassStdOut = runJobClassOutputs[0];
			String runJobClassStdErr = runJobClassOutputs[1];
			System.out.println("runJobClassStdOut = "+runJobClassStdOut);
			System.out.println("runJobClassStdErr = "+runJobClassStdErr);
			System.out.println();
			System.out.println("\nTrying to get the reads which were originally mapped to the contig " + contigName + "...");
			// trying to get the aligned reads for this specific contig by exploring the respective SAM file
			String com4ExtractContigFromBAMFile = "samtools view -h -F 4 -o " + contigName + ".sam " + sortedBAMFile.getAbsolutePath() + " " + contigName;
			//new RunJob_Antonio().runJob(com4ExtractContigFromBAMFile);
			runJobClassOutputs = new RunJob_Antonio().runJob(com4ExtractContigFromBAMFile);
			runJobClassStdOut = runJobClassOutputs[0];
			runJobClassStdErr = runJobClassOutputs[1];
			System.out.println("runJobClassStdOut = "+runJobClassStdOut);
			System.out.println("runJobClassStdErr = "+runJobClassStdErr);
			System.out.println();
			File greppedSAMFile = new File(contigName + ".sam");
			// reheading greppedSAMFile
			File samReHeadedFile = samReHeader(contigName, greppedSAMFile);
			String comSAM2BAM = "samtools view -S -b -o " + contigName + ".bam " + samReHeadedFile.getAbsolutePath();
			//new RunJob_Antonio().runJob(comSAM2BAM);
			runJobClassOutputs = new RunJob_Antonio().runJob(comSAM2BAM);
			runJobClassStdOut = runJobClassOutputs[0];
			runJobClassStdErr = runJobClassOutputs[1];
			System.out.println("runJobClassStdOut = "+runJobClassStdOut);
			System.out.println("runJobClassStdErr = "+runJobClassStdErr);
			System.out.println();
			// sorting BAM
			File convBAMFile = new File(contigName + ".bam");
			// indexing BAM
			System.out.println();
			if (debug) System.out.println("\nCreating the required file: " + contigName + ".bam.bai");
			String comIndexBAM = "samtools index " + convBAMFile.getAbsolutePath();
			//new RunJob_Antonio().runJob(comIndexBAM);
			runJobClassOutputs = new RunJob_Antonio().runJob(comIndexBAM);
			runJobClassStdOut = runJobClassOutputs[0];
			runJobClassStdErr = runJobClassOutputs[1];
			System.out.println("runJobClassStdOut = "+runJobClassStdOut);
			System.out.println("runJobClassStdErr = "+runJobClassStdErr);
			System.out.println();
			File baiFile = new File(contigName + ".bam.bai");
			for (Integer snpPos : snpLookup.get(contigName))
			{
				// a new report object for this SNP only
				if (debug)
					System.out.println();
				if (debug)
					System.out.println("\nEntering in the INNER LOOP step for snpPos " + snpPos + ".");
				SNPReport snpReport = new SNPReport();
				SNPReport2 snpReport2 = new SNPReport2();
				if (databaseFile == null)
				{
					snpReports.add(snpReport);
				}
				else
				{
					snpReports2.add(snpReport2);
				}
				// get the length of the contig
				int finalPosOfContig = sequence.length();
				if (debug)
					System.out.println("\nThe contig " + contigName + " has a length of " + finalPosOfContig + " bases.");
				char refAlleleAtPos = sequence.charAt((snpPos - 1));
				// extract allele counts and reads information
				IndexedFastaSequenceFile faidx = new IndexedFastaSequenceFile(new File(contigName + ".fasta"));
				if (debug)
					System.out.println("\nEntering in the alleleExtractorLight method step for contig " + contigName + " and SNP position " + snpPos + ".");
				// ----------------- Log variable to count the number of entrances in the alleleExtractorLight method -------------------------------------------------------------------------------------------------------------------------------------------------------
				accumEnteringAlleleExtractor++;
				// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				// The following was removed from the below method "alleleExtractorLight" call: contigFastaFile,
				alleleExtractorLight(new File(convBAMFile.getAbsolutePath()), snpPos.intValue(), faidx, baiFile, contigName, refAlleleAtPos, finalPosOfContig, snpReport, snpReport2, chrmLookup);
				// report/output this appropriately
				if (debug)
				{
					if (databaseFile == null)
						snpReport.outputToStdout();
					else
						snpReport2.outputToStdout();
				}
			}
		}
		System.out.println("\ndone - Processed and stored each entry of the given assembly file...");
		System.out.println("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =");
		if (debug)
			System.out.println("\nEnding the execution of the MAIN method... Seems promising...");
		if (databaseFile != null)
		{
			int sumAccumPresenceOfBothAllelesPLUSaccumPresenceOnlyAlternateAlleles = accumPresenceOfBothAlleles+accumPresenceOnlyAlternateAlleles;
			int sumAccumMisassemblyPLUSaccumNoMisassemblyPLUSaccumNotPossible2DetermineSNPPosOnGenome = accumMisassembly+accumNoMisassembly+accumNotPossible2DetermineSNPPosOnGenome;
			int sumAccumIsRevCompPLUSaccumIsNotRevCompPLUSaccumNotPossible2DetermineSNPPosOnGenome = accumIsRevComp+accumIsNotRevComp+accumNotPossible2DetermineSNPPosOnGenome;
			int sumAccumNonACTGEventsPLUSaccumNormalACTGEventsPLUSaccumNotPossible2DetermineSNPPosOnGenome = accumNonACTGEvents+accumNormalACTGEvents+accumNotPossible2DetermineSNPPosOnGenome;
			int sumAccumSNPSiteVisitedByBlastPLUSaccumSNPSiteNotVisitedByBlastPLUSaccumNotPossible2DetermineSNPPosOnGenome = accumSNPSiteVisitedByBlast+accumSNPSiteNotVisitedByBlast+accumNotPossible2DetermineSNPPosOnGenome;
			System.out.println("\nNumber of SNP events: "+ numberOfSNPEvents +" (at least two alleles present: "+ accumPresenceOfBothAlleles +" || Only one allele present: "+ accumPresenceOnlyAlternateAlleles + " || Total: " + sumAccumPresenceOfBothAllelesPLUSaccumPresenceOnlyAlternateAlleles +")");
			System.out.println("\nNumber of MISASSEMBLIES at the SNP site: "+ accumMisassembly +" || Number of NON-MISASSEMBLIES at the SNP site: "+ accumNoMisassembly + " || Number of NO BLAST results at all: "+ accumNotPossible2DetermineSNPPosOnGenome + " || Total: " + sumAccumMisassemblyPLUSaccumNoMisassemblyPLUSaccumNotPossible2DetermineSNPPosOnGenome);
			System.out.println("\nNumber of isRevComp events: "+ accumIsRevComp + " || Number of non-isRevComp events: "+accumIsNotRevComp + " || Number of NO BLAST results at all: "+ accumNotPossible2DetermineSNPPosOnGenome + " || Total: " + sumAccumIsRevCompPLUSaccumIsNotRevCompPLUSaccumNotPossible2DetermineSNPPosOnGenome); 
			System.out.println("\nNumber of non-ACTG events: "+ accumNonACTGEvents + " || Number of ACTG events: "+ accumNormalACTGEvents + " || Number of NO BLAST results at all: "+ accumNotPossible2DetermineSNPPosOnGenome + " || Total: " + sumAccumNonACTGEventsPLUSaccumNormalACTGEventsPLUSaccumNotPossible2DetermineSNPPosOnGenome);
			System.out.println("\nNumber of SNP positions COVERED by BLAST: "+ accumSNPSiteVisitedByBlast +" || Number of SNP positions NOT COVERED by BLAST: "+ accumSNPSiteNotVisitedByBlast+ " || Number of NO BLAST results at all: "+ accumNotPossible2DetermineSNPPosOnGenome + " || Total: " + sumAccumSNPSiteVisitedByBlastPLUSaccumSNPSiteNotVisitedByBlastPLUSaccumNotPossible2DetermineSNPPosOnGenome);				
			// Now, zeroing the last used variables....
			accumMisassembly = 0;
			accumNoMisassembly = 0;
			accumIsRevComp = 0;
			accumIsNotRevComp = 0;
			accumNotPossible2DetermineSNPPosOnGenome = 0;
			accumNonACTGEvents = 0;
			accumNormalACTGEvents = 0;
			accumSNPSiteVisitedByBlast = 0;
			accumSNPSiteNotVisitedByBlast = 0;
			sumAccumPresenceOfBothAllelesPLUSaccumPresenceOnlyAlternateAlleles = 0;
			sumAccumMisassemblyPLUSaccumNoMisassemblyPLUSaccumNotPossible2DetermineSNPPosOnGenome = 0;
			sumAccumIsRevCompPLUSaccumIsNotRevCompPLUSaccumNotPossible2DetermineSNPPosOnGenome = 0;
			sumAccumNonACTGEventsPLUSaccumNormalACTGEventsPLUSaccumNotPossible2DetermineSNPPosOnGenome = 0;
			sumAccumSNPSiteVisitedByBlastPLUSaccumSNPSiteNotVisitedByBlastPLUSaccumNotPossible2DetermineSNPPosOnGenome = 0;
			// accumAveragePercentageMismappedReadsPerSNPEvent = 0;		 
		}
		numberOfSNPEvents = 0;
		accumPresenceOfBothAlleles = 0;
		accumPresenceOnlyAlternateAlleles = 0;
		accumEnteringAlleleExtractor = 0;
		System.out.println("\nCleaning the unnecessary intermediate .sam files...");
		String cmd2DeleteSAMs = "../deleteSAMs.sh";
		//new RunJob_Antonio().runJob(cmd2DeleteSAMs);
		String [] runJobClassOutputs = new RunJob_Antonio().runJob(cmd2DeleteSAMs);
		String runJobClassStdOut = runJobClassOutputs[0];
		String runJobClassStdErr = runJobClassOutputs[1];
		System.out.println("runJobClassStdOut = "+runJobClassStdOut);
		System.out.println("runJobClassStdErr = "+runJobClassStdErr);
		System.out.println();
		System.out.println("\nGood job!!! Tool execution FINISHED smoothly! HAVE A NICE DAY!!!");
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	private static ArrayList<String> parseSNPList(File listOfSNPsFile)
	{
		System.out.println("\nReading a given list Of contigs and SNPs file and populating a list of its entries...");
		BufferedReader br = null;
		ArrayList<String> listTest = new ArrayList<String>();
		//String nameOfFileWOExtension = stripExtension(listOfSNPsFile.getName());
		//int counterUsefulLinesOfSNPsFile = 0;
		int counterRejectedLinesOfSNPsFile = 0;
		try
		{
			// reads each line of the file and populate the array list named listTest
			String line;
			br = new BufferedReader(new FileReader(listOfSNPsFile));
			while ((line = br.readLine()) != null)
			{
				String[] lineArray;
				lineArray = line.split("\t");
				String contigName = lineArray[0];
				//System.out.println("\ncontigName = "+contigName);
				int refSeqPos = Integer.parseInt(lineArray[1]);
				//System.out.println("\nrefSeqPos = "+refSeqPos);
				String refFieldContent = lineArray[3];
				//System.out.println("\nrefFieldContent = "+refFieldContent);
				String altFieldContent = lineArray[4];
				//System.out.println("\naltFieldContent = "+altFieldContent);
				//if (!(line.startsWith("SNP") || line.startsWith(nameOfFileWOExtension + ".vcf")))
				if (refFieldContent.length() == 1 && altFieldContent.length() == 1)
				{
					//counterUsefulLinesOfSNPsFile++;
					String contigNameFull = contigName+"_"+refSeqPos;
					//lineArray = line.split("\t");
					//String contigNameFull = lineArray[0];
					//String positionLastUnderscoreContigName = contigNameFull.substring(contigNameFull.lastIndexOf("_") + 1);
					//int posLastUnderscoreContigName = Integer.parseInt(positionLastUnderscoreContigName);
					// if (debug) System.out.println("\nposLastUnderscoreContigName: "+posLastUnderscoreContigName);
					//String contigNameWOPos = contigNameFull.substring(0, contigNameFull.lastIndexOf("_"));
					// if (debug) System.out.println("\ncontigNameWOPos: "+contigNameWOPos);
					//int refSeqPos = Integer.parseInt(lineArray[2]);
					// if (debug) System.out.println("\nrefSeqPos: "+refSeqPos);
					// String isHomozygous = lineArray[8];
					//if (posLastUnderscoreContigName == refSeqPos)
					//{
						// if (debug) System.out.println("\nCondition posLastUnderscoreContigName == refSeqPos");
						listTest.add(contigNameFull.trim());
					//}
					//else
					//{
						//if (debug)
							//System.out.println("\nTHIS IS NOT SO PROMISING!...For contig " + contigNameWOPos + ", the SNP site " + posLastUnderscoreContigName + " is NOT CONSISTENT with the reference sequence position " + refSeqPos + ". Please check this issue!!!");
						//counterRejectedLinesOfSNPsFile++;
						//if (debug)
							//System.out.println("counterRejectedLinesOfSNPsFile: " + counterRejectedLinesOfSNPsFile);
					//}
				}
				else 
				{
					counterRejectedLinesOfSNPsFile++;
				}
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
		//if (!(counterUsefulLinesOfSNPsFile == listTest.size()) && (counterRejectedLinesOfSNPsFile == 0))
			//listTest.clear();
		//System.out.println("\ndone -- The processed file had " + listTest.size() + " lines and was internally CONSISTENT!");
		System.out.println("\ndone -- The processed file had " + listTest.size() + " valid input lines and "+counterRejectedLinesOfSNPsFile+" which were rejected...");
		//counterUsefulLinesOfSNPsFile = 0;
		counterRejectedLinesOfSNPsFile = 0;
		return listTest;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// make lookup of contigs (keys) versus lists of SNPs within these (values) from a given snpList produced by the parseSNPList method
	private static HashMap<String, TreeSet<Integer>> makeSNPLookup(ArrayList<String> snpList)
	{
		System.out.println("\nMaking SNP Lookup 'data structure'...");
		// Stores contig names as keys and lists of SNPs positions as values
		HashMap<String, TreeSet<Integer>> snpLookup = new HashMap<String, TreeSet<Integer>>();
		// NODE_10001_length_277_cov_91.635376_143
		// or
		// NODE_10009_length_1777_cov_26.356218_1125
		// NODE_10009_length_1777_cov_26.356218_1210
		// iterates over list of SNPs
		for (String entryName : snpList)
		{
			// extracting contig names and SNPs positions list
			String contigName = entryName.substring(0, entryName.lastIndexOf("_"));
			String position = entryName.substring(entryName.lastIndexOf("_") + 1);
			// get the list of SNP positions for this contig and check whether it exists
			TreeSet<Integer> snpPositions = snpLookup.get(contigName);
			// if the contig doesn't exist, creates a new entry for it and, also, for its snpPositions list
			if (snpPositions == null)
			{
				// new instance of snpPositions object
				snpPositions = new TreeSet<Integer>();
				// adding the list of contigs and respective snpPositions
				snpLookup.put(contigName, snpPositions);
			}
			// adding the position to list of snpPositions for this contig
			snpPositions.add(Integer.parseInt(position));
		}
		System.out.println("\ndone -- Entries processed: " + snpLookup.size() + ".");
		return snpLookup;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	public static File outputSeq2Fasta(String sequenceName, String sequence) throws IOException
	{
		// naming the file with the sequence name and adding the .fasta extension
		File fastaFile = new File(sequenceName + ".fasta");
		// opening an output Stream to that file
		BufferedWriter bw = new BufferedWriter(new FileWriter(fastaFile));
		// writing sequence name FASTA style
		bw.write(">" + sequenceName);
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
	// private static void alleleExtractorLight(File splittedBAMFile, int snpPos, IndexedFastaSequenceFile faidx, File splittedBAMBAIFile, String contigName, char refAlleleAtPos, int
	// finalPosOfContig, SNPReport snpReport)
	// The following item was removed from the "alleleExtractorLight" method below: File contigFastaFile
	private static void alleleExtractorLight(File contigNameSortedBAMFile, int snpPos, IndexedFastaSequenceFile faidx, File baiFile, String contigName, char refAlleleAtPos, int finalPosOfContig, SNPReport snpReport, SNPReport2 snpReport2, HashMap<String, String> chrmLookup)
	{
		System.out.println("\nExtracting alleles and reads information STEP...");
		try
		{
			// ----------------------------------------- New block of code created specifically for the "RegionQuantifier" class ------------------------------------------------------------------------------------------------------------------------------------------
			if (debug)
				System.out.println("\nENTERING in the 'BLAST of assembled contig against the dbFile to retrieve the original genomic coordinates of the FPSNP site' substep of the 'Region Quantifier' module...");
			boolean isRevComp = false;
			int dist2SNPFromEndOfContig = 0;
			int snpPosOnGenome = 0;
			boolean misAssembly = false;
			String chrmInGenome = "";
			String chrmID = "";
			int lowerBound = 0;
			int upperBound = 0;
			int accumFragFromDiffChrm = 0;
			int accumFragFromSameChrm = 0;
			int accumFragFromWrongRegionSameChrm = 0;
			int accumMismappedRead = 0;
			int accumCorrectlyMappedRead = 0;
			int accumMismatchedRead = 0;
			int accumMatchedRead = 0;
			Character refAlleleOnGenome4Comparison = null;
			Character refAlleleAtPos4Comparison = null;
			char refAlleleOnGenome = ' ';
			char complementedRefAlleleOnGenome = ' ';
			int accumMisAssemblyMismatchedRead = 0;
			int accumMisAssemblyMatchedRead = 0;
			int accumMisAssemblyDiffChrmMismapMatchedRead = 0;
			int accumMisAssemblyDiffChrmMismapMismatchedRead = 0;
			int accumMisAssemblySameChrmCorrectlyMappedMatchedRead = 0;
			int accumMisAssemblySameChrmCorrectlyMappedMismatchedRead = 0;
			int accumMisAssemblySameChrmMismapMatchedRead = 0;
			int accumMisAssemblySameChrmMismapMismatchedRead = 0;
			int accumCorrectAssyMismatchedRead = 0;
			int accumCorrectAssyMatchedRead = 0;
			int accumCorrectAssyDiffChrmMismapMatchedRead = 0;
			int accumCorrectAssyDiffChrmMismapMismatchedRead = 0;
			int accumCorrectAssySameChrmCorrectlyMappedMatchedRead = 0;
			int accumCorrectAssySameChrmCorrectlyMappedMismatchedRead = 0;
			int accumCorrectAssySameChrmMismapMatchedRead = 0;
			int accumCorrectAssySameChrmMismapMismatchedRead = 0;
			int accumOddSituation4MisAssemblyScenario = 0;
			int accumOddSituation4CorrectAssyScenario = 0;
			int accumCasesOfCorrectAssyButMismapMismatch = 0;
			int accumCasesOfMisAssemblyAndMismapMatch = 0;
			int accumNetMismappedReads = 0;
			float percentageOfNetMismappedReads = 0;
			int sStart = 0;
			int sEnd = 0;
			// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			// Flags created for debugging purposes...
			boolean blastResultIsEmpty = false;
			boolean nonACTGEvent = false;
			boolean fpSNPSiteNotVisitedByBlast = false;
			// Helper variables created for the "nonACTGEvent" test block
			final char charA = 'A';
			final char charC = 'C';
			final char charT = 'T';
			final char charG = 'G';
			// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			if (databaseFile != null)
			{
				if (debug)
					System.out.println("\nENTERING the 'Retrieve the original genomic coordinates of the FPSNP site without BLAST...' substep aiming to fetch the hardwired genomic coordinates for the CONTROLS...");
				if (contigName.equalsIgnoreCase("Chr1"))
				{
					sEnd = 30427671;
				}
				else if (contigName.equalsIgnoreCase("Chr2"))
				{
					sEnd = 19698289;
				}
				else if (contigName.equalsIgnoreCase("Chr3"))
				{
					sEnd = 23459830;
				}
				else if (contigName.equalsIgnoreCase("Chr4"))
				{
					sEnd = 18585056;
				}
				else if (contigName.equalsIgnoreCase("Chr5"))
				{
					sEnd = 26975502;
				}
				sStart = 1;
				chrmInGenome = contigName;
				System.out.println("\nOriginal chromosome in genome is: " + chrmInGenome + ".");
				System.out.println("\nContig START position in genome is: " + sStart + ".");
				System.out.println("\nContig END position in genome is: " + sEnd + ".");
				if (sStart > sEnd)
				{
					isRevComp = true;
					accumIsRevComp++;
				}
				if (debug)
					System.out.println("\nContig is in REVCOMP test: " + isRevComp + ".");
				if (isRevComp)
				{
					dist2SNPFromEndOfContig = finalPosOfContig - snpPos;
					if (debug)
						System.out.println("\ndist2SNPFromEndOfContig with REVCOMP condition " + isRevComp + " is: " + dist2SNPFromEndOfContig + ".");
					snpPosOnGenome = sEnd + dist2SNPFromEndOfContig;
					if (debug)
						System.out.println("\nsnpPosOnGenome with REVCOMP condition " + isRevComp + " is: " + snpPosOnGenome + ".");
				}
				else
				{
					dist2SNPFromEndOfContig = finalPosOfContig - snpPos;
					snpPosOnGenome = sStart + snpPos - 1;
					if (debug)
						System.out.println("\nsnpPosOnGenome with REVCOMP condition " + isRevComp + " is: " + snpPosOnGenome + ".");
					accumIsNotRevComp++;
				}
				if (debug)
					System.out.println("\nEXITING the 'Retrieve the original genomic coordinates of the FPSNP site without BLAST...' substep aiming to fetch the hardwired genomic coordinates for the CONTROLS...");
				if (debug)
					System.out.println("\nENTERING the 'Retrieve chromosome sequence from the 'chrmLookup' data structure for the verification of the allele present in the FPSNP site' substep aiming to decide whether the scenario is of a MISASSEMBLY of the contig or not...");
				for (String chrm : chrmLookup.keySet())
				{
					if (chrmInGenome.equals(chrm))
					{
						// look up its sequence in the sequence lookup
						String tmpSequence = chrmLookup.get(chrm);
						if (debug)
						{
							char refAlleleOnGenomeMinus5 = tmpSequence.charAt(snpPosOnGenome - 5);
							if (isRevComp)
								refAlleleOnGenomeMinus5 = complementAllele(refAlleleOnGenomeMinus5);
							System.out.println("\nrefAlleleOnGenomeMinus5 = " + refAlleleOnGenomeMinus5);
							char refAlleleOnGenomeMinus4 = tmpSequence.charAt(snpPosOnGenome - 4);
							if (isRevComp)
								refAlleleOnGenomeMinus4 = complementAllele(refAlleleOnGenomeMinus4);
							System.out.println("refAlleleOnGenomeMinus4 = " + refAlleleOnGenomeMinus4);
							char refAlleleOnGenomeMinus3 = tmpSequence.charAt(snpPosOnGenome - 3);
							if (isRevComp)
								refAlleleOnGenomeMinus3 = complementAllele(refAlleleOnGenomeMinus3);
							System.out.println("refAlleleOnGenomeMinus3 = " + refAlleleOnGenomeMinus3);
							char refAlleleOnGenomeMinus2 = tmpSequence.charAt(snpPosOnGenome - 2);
							if (isRevComp)
								refAlleleOnGenomeMinus2 = complementAllele(refAlleleOnGenomeMinus2);
							System.out.println("refAlleleOnGenomeMinus2 = " + refAlleleOnGenomeMinus2);
							char refAlleleOnGenomeMinus1 = tmpSequence.charAt(snpPosOnGenome - 1);
							if (isRevComp)
								refAlleleOnGenomeMinus1 = complementAllele(refAlleleOnGenomeMinus1);
							System.out.println("refAlleleOnGenomeMinus1 = " + refAlleleOnGenomeMinus1);
						}
						refAlleleOnGenome = tmpSequence.charAt(snpPosOnGenome - 1);
						if (isRevComp)
						{
							System.out.println("\nisRevComp status is: " + isRevComp + ", then I'm complementing the genomic found base " + refAlleleOnGenome + " just for computational purposes here...");
							complementedRefAlleleOnGenome = complementAllele(refAlleleOnGenome);
							System.out.println("\ncomplementedRefAlleleOnGenome = " + complementedRefAlleleOnGenome);
						}
						else
							System.out.println("\nrefAlleleOnGenome (should be equal to the allele value of snpPosOnGenomeMinus1) = " + refAlleleOnGenome);
						if (debug)
						{
							char refAlleleOnGenomePlus0 = tmpSequence.charAt(snpPosOnGenome);
							if (isRevComp)
								refAlleleOnGenomePlus0 = complementAllele(refAlleleOnGenomePlus0);
							System.out.println("\nrefAlleleOnGenomePlus0 = " + refAlleleOnGenomePlus0);
							char refAlleleOnGenomePlus1 = tmpSequence.charAt(snpPosOnGenome + 1);
							if (isRevComp)
								refAlleleOnGenomePlus1 = complementAllele(refAlleleOnGenomePlus1);
							System.out.println("refAlleleOnGenomePlus1 = " + refAlleleOnGenomePlus1);
							char refAlleleOnGenomePlus2 = tmpSequence.charAt(snpPosOnGenome + 2);
							if (isRevComp)
								refAlleleOnGenomePlus2 = complementAllele(refAlleleOnGenomePlus2);
							System.out.println("refAlleleOnGenomePlus2 = " + refAlleleOnGenomePlus2);
							char refAlleleOnGenomePlus3 = tmpSequence.charAt(snpPosOnGenome + 3);
							if (isRevComp)
								refAlleleOnGenomePlus3 = complementAllele(refAlleleOnGenomePlus3);
							System.out.println("refAlleleOnGenomePlus3 = " + refAlleleOnGenomePlus3);
							char refAlleleOnGenomePlus4 = tmpSequence.charAt(snpPosOnGenome + 4);
							if (isRevComp)
								refAlleleOnGenomePlus4 = complementAllele(refAlleleOnGenomePlus4);
							System.out.println("refAlleleOnGenomePlus4 = " + refAlleleOnGenomePlus4);
							char refAlleleOnGenomePlus5 = tmpSequence.charAt(snpPosOnGenome + 5);
							if (isRevComp)
								refAlleleOnGenomePlus5 = complementAllele(refAlleleOnGenomePlus5);
							System.out.println("refAlleleOnGenomePlus5 = " + refAlleleOnGenomePlus5);
						}
						if (isRevComp)
							refAlleleOnGenome4Comparison = Character.toLowerCase(complementedRefAlleleOnGenome);
						else
							refAlleleOnGenome4Comparison = Character.toLowerCase(refAlleleOnGenome);
						refAlleleAtPos4Comparison = Character.toLowerCase(refAlleleAtPos);
						if (!refAlleleOnGenome4Comparison.equals(refAlleleAtPos4Comparison))
						{
							misAssembly = true;
							accumMisassembly++;
						}
						else
							accumNoMisassembly++;
                                                // ------------------------------------------------------- Block created to quantify occurrences of non ACTG events --------------------------------------------------------------------------------------------------------------------------------------------
						if (debug) System.out.println("refAlleleOnGenome4Comparison = "+refAlleleOnGenome4Comparison);
						if (!refAlleleOnGenome4Comparison.equals(Character.toLowerCase(charA)) && !refAlleleOnGenome4Comparison.equals(Character.toLowerCase(charC)) && !refAlleleOnGenome4Comparison.equals(Character.toLowerCase(charT)) && !refAlleleOnGenome4Comparison.equals(Character.toLowerCase(charG)))
						{
							nonACTGEvent = true;
							accumNonACTGEvents++;
						}
						else
							accumNormalACTGEvents++;
						// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
						// ------------------------------------------------------- Block created to quantify occurrences where the FPSNP site (on the genome) occurs in a region for which BLAST didn't encounter any similarities --------------------------------------------------------------------------------------------------------------------------------------------
						if (isRevComp) 
						{
					        	if ((snpPosOnGenome -1 > sStart) || (snpPosOnGenome -1 < sEnd))
					        	{
					        		fpSNPSiteNotVisitedByBlast = true;
					        		accumSNPSiteNotVisitedByBlast++;
					        	}
					        	else
					        		accumSNPSiteVisitedByBlast++;
						}
						else
						{
							if ((snpPosOnGenome -1 < sStart) || (snpPosOnGenome -1 > sEnd))
					        	{
					        		fpSNPSiteNotVisitedByBlast = true;
					        		accumSNPSiteNotVisitedByBlast++;
					        	}
					        	else
					        		accumSNPSiteVisitedByBlast++;
						}
						// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
					}
					// For some very few events, not all the information will be retrieved at this point... These "rows" will be further disconsidered, for the sake of mismapping calculations, when computing the results of the output file using spreadsheet proper filters.
				}
				if (debug)
					System.out.println("\nmisAssembly status of the contig is: " + misAssembly + " => refAlleleOnGenome4Comparison = " + refAlleleOnGenome4Comparison + " / refAlleleAtPos4Comparison (of contig) = " + refAlleleAtPos4Comparison + ".");
				if (debug)
					System.out.println("\nEXITING the 'Retrieve chromosome sequence from the 'chrmLookup' data structure for the verification of the allele present in the FPSNP site' substep...");
			}
			// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			char refAlleleAtLast = refAlleleAtPos;
			SAMFileReader samReader = new SAMFileReader(contigNameSortedBAMFile, baiFile, true);
			SAMFileHeader header = samReader.getFileHeader();
			List<SAMSequenceRecord> sequences = header.getSequenceDictionary().getSequences();
			HashMap<String, String> overlappingReadsLookup = new HashMap<String, String>();
			HashMap<String, Character> readsAndAllelesLookup = new HashMap<String, Character>();
			int accumBaseA = 0;
			int accumBaseC = 0;
			int accumBaseT = 0;
			int accumBaseG = 0;
			int accumUndetermined = 0;
			int acummRefAllele = 0;
			int totalDiffAlleles = 0;
			int numberOfOverlappingReads = 0;
			int counterOccurrenceAlleles = 0;
			boolean presenceOfBothAlleles = false;
			boolean presenceOnlyAlternateAlleles = false;
			boolean presenceOnlyReferenceAllele = false;
			boolean presenceOfOddOccurrence = false;
			boolean presenceOf3OrMoreAlleles = false;
			int[] testPresenceFor3OrMoreAlleles = new int[5];
			for (int i = 1; i == 5; i++)
			{
				testPresenceFor3OrMoreAlleles[i] = 0;
			}
			// iterate over all the contigs in the references
			for (SAMSequenceRecord rec : sequences)
			{
				String seqName = rec.getSequenceName();
				if (debug)
					System.out.println("\nseqName: " + seqName);
				// an iterator over all the reads that overlap the SNP position
				samReader.setValidationStringency(ValidationStringency.LENIENT);
				SAMRecordIterator samRecordIterator = samReader.queryOverlapping(seqName, snpPos, snpPos);
				System.out.println("\nThe following reads are overlapping the position " + snpPos + ":");
				System.out.println();
				while (samRecordIterator.hasNext())
				{
					SAMRecord record = samRecordIterator.next();
					String readName = record.getReadName();
					if (debug)
						System.out.println(readName);
					numberOfOverlappingReads++;
					// block to create a hashmap from the UNIQUE overlapping reads only
					String readSequence = record.getReadString();
					overlappingReadsLookup.put(readSequence, readName);
					// extract the allele at the SNP position
					int readStart = record.getAlignmentStart();
					int snpPosOnRead = snpPos - readStart + 1;
					if (debug)
						System.out.println("\noriginal read:");
					if (debug)
						System.out.println(record.getReadString());
					String cigarCorrectedRead = getCIGARCorrectedRead(record, faidx);
					if (debug)
						System.out.println("cigar corrected read:");
					if (debug)
						System.out.println(cigarCorrectedRead);
					if (debug)
						System.out.println("cigar string:");
					if (debug)
						System.out.println(record.getCigarString());
					if (debug)
						System.out.println("snpPosOnRead = " + snpPosOnRead);
					if (debug)
						System.out.println("readStart = " + readStart);
					char allele = cigarCorrectedRead.charAt(snpPosOnRead - 1);
					System.out.println("read " + readName + " has allele " + allele + ".");
					// block to create a hashmap of the reads and respective alleles at the SNP position
					readsAndAllelesLookup.put(readName, allele);
					// ----------------------------------------- New block of code created specifically for the "RegionQuantifier" class ------------------------------------------------------------------------------------------------------------------------------------------
					Character allele4Comparison = Character.toLowerCase(allele);
					Character refAlleleAtLast4Comparison = Character.toLowerCase(refAlleleAtLast);
					if (databaseFile != null)
					{
						if (debug)
							System.out.println("\nENTERING the 'Mismapped / Mismatched read occurrences computation' substep for the read...");
						System.out.println();
						// ####### Block changed when the read simulator was changed from SHERMAN to SimSeq ############################
						//SimSeq read label example: Chr1_7357887_7357972_R_14f38f6_R2
						//For old Sherman read label: chrmID = readName.substring(readName.indexOf("_") + 1, readName.lastIndexOf(":"));
						String token[] = readName.split("_");
						chrmID = token[0];
						if (debug)
							System.out.println("chrmID = " + chrmID);
						lowerBound = Integer.parseInt(token[1]);
						upperBound = Integer.parseInt(token[2]);
						// ##############################################################################################################
						//For old Sherman read label: lowerBound = Integer.parseInt(readName.substring(readName.indexOf(":") + 1, readName.lastIndexOf("-")));
						if (debug)
							System.out.println("lowerBound = " + lowerBound);
						//For old Sherman read label: upperBound = Integer.parseInt(readName.substring(readName.indexOf("-") + 1, readName.lastIndexOf("_")));
						if (debug)
							System.out.println("upperBound = " + upperBound);
						if (misAssembly)
						{
							if (!allele4Comparison.equals(refAlleleAtLast4Comparison) && allele4Comparison.equals(refAlleleOnGenome4Comparison))
							{
								accumMismatchedRead++;
								accumMisAssemblyMismatchedRead++;
								if (!chrmID.equalsIgnoreCase(chrmInGenome))
								{
									accumMismappedRead++;
									accumFragFromDiffChrm++;
									accumMisAssemblyDiffChrmMismapMismatchedRead++;
								}
								else
								{
									accumFragFromSameChrm++;
									if ((lowerBound <= snpPosOnGenome) && (upperBound >= snpPosOnGenome))
									{
										accumCorrectlyMappedRead++;
										accumMisAssemblySameChrmCorrectlyMappedMismatchedRead++;
									}
									else
									{
										accumMismappedRead++;
										accumMisAssemblySameChrmMismapMismatchedRead++;
										accumFragFromWrongRegionSameChrm++;
									}
								}
							}
							else if (allele4Comparison.equals(refAlleleAtLast4Comparison) && !allele4Comparison.equals(refAlleleOnGenome4Comparison))
							{
								accumMatchedRead++;
								accumMisAssemblyMatchedRead++;
								if (!chrmID.equalsIgnoreCase(chrmInGenome))
								{
									accumMismappedRead++;
									accumFragFromDiffChrm++;
									accumMisAssemblyDiffChrmMismapMatchedRead++;
								}
								else
								{
									accumFragFromSameChrm++;
									if ((lowerBound <= snpPosOnGenome) && (upperBound >= snpPosOnGenome))
									{
										accumCorrectlyMappedRead++;
										accumMisAssemblySameChrmCorrectlyMappedMatchedRead++;
									}
									else
									{
										accumMismappedRead++;
										accumMisAssemblySameChrmMismapMatchedRead++;
										accumFragFromWrongRegionSameChrm++;
									}
								}
							}
							else
							{
								accumOddSituation4MisAssemblyScenario++;
								accumMismatchedRead++;
								accumMisAssemblyMismatchedRead++;
								if (!chrmID.equalsIgnoreCase(chrmInGenome))
								{
									accumMismappedRead++;
									accumFragFromDiffChrm++;
									accumMisAssemblyDiffChrmMismapMismatchedRead++;
								}
								else
								{
									accumFragFromSameChrm++;
									if ((lowerBound <= snpPosOnGenome) && (upperBound >= snpPosOnGenome))
									{
										accumCorrectlyMappedRead++;
										accumMisAssemblySameChrmCorrectlyMappedMismatchedRead++;
									}
									else
									{
										accumMismappedRead++;
										accumMisAssemblySameChrmMismapMismatchedRead++;
										accumFragFromWrongRegionSameChrm++;
									}
								}							
								if (debug)
									System.out.println("\nOddSituation4MisAssemblyScenario...");
							}
						}
						else
						{
							if (!allele4Comparison.equals(refAlleleAtLast4Comparison) && !allele4Comparison.equals(refAlleleOnGenome4Comparison))
							{
								accumMismatchedRead++;
								accumCorrectAssyMismatchedRead++;
								if (!chrmID.equalsIgnoreCase(chrmInGenome))
								{
									accumMismappedRead++;
									accumFragFromDiffChrm++;
									accumCorrectAssyDiffChrmMismapMismatchedRead++;
								}
								else
								{
									accumFragFromSameChrm++;
									if ((lowerBound <= snpPosOnGenome) && (upperBound >= snpPosOnGenome))
									{
										accumCorrectlyMappedRead++;
										accumCorrectAssySameChrmCorrectlyMappedMismatchedRead++;
									}
									else
									{
										accumMismappedRead++;
										accumCorrectAssySameChrmMismapMismatchedRead++;
										accumFragFromWrongRegionSameChrm++;
									}
								}
							}
							else if (allele4Comparison.equals(refAlleleAtLast4Comparison) && allele4Comparison.equals(refAlleleOnGenome4Comparison))
							{
								accumMatchedRead++;
								accumCorrectAssyMatchedRead++;
								if (!chrmID.equalsIgnoreCase(chrmInGenome))
								{
									accumMismappedRead++;
									accumFragFromDiffChrm++;
									accumCorrectAssyDiffChrmMismapMatchedRead++;
								}
								else
								{
									accumFragFromSameChrm++;
									if ((lowerBound <= snpPosOnGenome) && (upperBound >= snpPosOnGenome))
									{
										accumCorrectlyMappedRead++;
										accumCorrectAssySameChrmCorrectlyMappedMatchedRead++;
									}
									else
									{
										accumMismappedRead++;
										accumCorrectAssySameChrmMismapMatchedRead++;
										accumFragFromWrongRegionSameChrm++;
									}
								}
							}
							else
							{
								accumOddSituation4CorrectAssyScenario++;
								if (debug)
									System.out.println("\nOddSituation4CorrectAssyScenario...");
							}
						}
						if (debug)
							System.out.println("\nFinally, EXITING the 'Mismapped / Mismatched read occurrences computation' substep for the read...");
						System.out.println();
					}
					// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
					infoAboutReadAndAllele.put(contigName + "_" + snpPos, readsAndAllelesLookup);
					// Trying to quantify the different alleles found in each read overlapping the SNP position
					if ((allele == 'A') || (allele == 'a'))
					{
						accumBaseA++;
					}
					else if ((allele == 'C') || (allele == 'c'))
					{
						accumBaseC++;
					}
					else if ((allele == 'T') || (allele == 't'))
					{
						accumBaseT++;
					}
					else if ((allele == 'G') || (allele == 'g'))
					{
						accumBaseG++;
					}
					else
					{
						accumUndetermined++;
					}
				}
				samReader.close();
				samRecordIterator.close();
			}
			testPresenceFor3OrMoreAlleles[0] = accumBaseA;
			testPresenceFor3OrMoreAlleles[1] = accumBaseC;
			testPresenceFor3OrMoreAlleles[2] = accumBaseT;
			testPresenceFor3OrMoreAlleles[3] = accumBaseG;
			testPresenceFor3OrMoreAlleles[4] = accumUndetermined;
			for (int i = 0; i < testPresenceFor3OrMoreAlleles.length; i++)
			{
				if (testPresenceFor3OrMoreAlleles[i] > 0)
				{
					counterOccurrenceAlleles++;
				}
			}
			if (counterOccurrenceAlleles >= 3)
			{
				presenceOf3OrMoreAlleles = true;
			}
			// System.out.println();
			if ((refAlleleAtLast == 'A') || (refAlleleAtLast == 'a'))
			{
				acummRefAllele = accumBaseA;
				totalDiffAlleles = accumBaseC + accumBaseT + accumBaseG + accumUndetermined;
			}
			else if ((refAlleleAtLast == 'C') || (refAlleleAtLast == 'c'))
			{
				acummRefAllele = accumBaseC;
				totalDiffAlleles = accumBaseA + accumBaseT + accumBaseG + accumUndetermined;
			}
			else if ((refAlleleAtLast == 'T') || (refAlleleAtLast == 't'))
			{
				acummRefAllele = accumBaseT;
				totalDiffAlleles = accumBaseA + accumBaseC + accumBaseG + accumUndetermined;
			}
			else if ((refAlleleAtLast == 'G') || (refAlleleAtLast == 'g'))
			{
				acummRefAllele = accumBaseG;
				totalDiffAlleles = accumBaseA + accumBaseC + accumBaseT + accumUndetermined;
			}
			else
			{
				acummRefAllele = accumUndetermined;
				totalDiffAlleles = accumBaseA + accumBaseC + accumBaseT + accumBaseG;
			}
			// populate lookup of transcripts/contigs (keys) versus positions and respective lists of quantities of each found allele
			if ((acummRefAllele > 0) && (totalDiffAlleles > 0))
			{
				presenceOfBothAlleles = true;
				accumPresenceOfBothAlleles++;
			}
			else if ((acummRefAllele == 0) && (totalDiffAlleles > 0))
			{
				presenceOnlyAlternateAlleles = true;
				accumPresenceOnlyAlternateAlleles++;
			}
			else if ((acummRefAllele > 0) && (totalDiffAlleles == 0))
			{
				presenceOnlyReferenceAllele = true;
			}
			else
			{
				presenceOfOddOccurrence = true;
			}
			int contigLength = finalPosOfContig;
			// ----------------------------------------- New block of code created specifically for the "RegionQuantifier" class
			// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			if (databaseFile != null)
			{
				accumCasesOfCorrectAssyButMismapMismatch = accumCorrectAssySameChrmMismapMismatchedRead + accumCorrectAssyDiffChrmMismapMismatchedRead;
				accumCasesOfMisAssemblyAndMismapMatch = accumMisAssemblySameChrmMismapMatchedRead + accumMisAssemblyDiffChrmMismapMatchedRead;
				if (misAssembly)
				{
					accumNetMismappedReads = accumCasesOfMisAssemblyAndMismapMatch;
					percentageOfNetMismappedReads = Math.round(100 * ((float) accumNetMismappedReads / (float) accumMatchedRead));
				}
				else
				{
					accumNetMismappedReads = accumCasesOfCorrectAssyButMismapMismatch;
					percentageOfNetMismappedReads = Math.round(100 * ((float) accumNetMismappedReads / (float) accumMismatchedRead));
				}
				if (debug)
					System.out.println("For debugging purposes, printing the temporary result of percentageOfNetMismappedReads: " + percentageOfNetMismappedReads + ".");
				populateAlleleQuantitiesLookup2(contigName, contigLength, snpPos, refAlleleAtLast, accumBaseA, accumBaseC, accumBaseT, accumBaseG, accumUndetermined, acummRefAllele, totalDiffAlleles, numberOfOverlappingReads, presenceOf3OrMoreAlleles, presenceOfBothAlleles, presenceOnlyAlternateAlleles, presenceOnlyReferenceAllele, presenceOfOddOccurrence, sStart, sEnd, blastResultIsEmpty, nonACTGEvent, fpSNPSiteNotVisitedByBlast, isRevComp, dist2SNPFromEndOfContig, snpPosOnGenome, misAssembly, chrmInGenome, accumFragFromDiffChrm, accumFragFromSameChrm, accumFragFromWrongRegionSameChrm, accumMismappedRead, accumCorrectlyMappedRead, accumMismatchedRead, accumMatchedRead, refAlleleOnGenome4Comparison, accumMisAssemblyMismatchedRead, accumMisAssemblyMatchedRead, accumMisAssemblyDiffChrmMismapMatchedRead, accumMisAssemblyDiffChrmMismapMismatchedRead, accumMisAssemblySameChrmCorrectlyMappedMatchedRead, accumMisAssemblySameChrmCorrectlyMappedMismatchedRead, accumMisAssemblySameChrmMismapMatchedRead, accumMisAssemblySameChrmMismapMismatchedRead, accumCorrectAssyMismatchedRead, accumCorrectAssyMatchedRead, accumCorrectAssyDiffChrmMismapMatchedRead, accumCorrectAssyDiffChrmMismapMismatchedRead, accumCorrectAssySameChrmCorrectlyMappedMatchedRead, accumCorrectAssySameChrmCorrectlyMappedMismatchedRead, accumCorrectAssySameChrmMismapMatchedRead, accumCorrectAssySameChrmMismapMismatchedRead, accumOddSituation4MisAssemblyScenario, accumOddSituation4CorrectAssyScenario, accumNetMismappedReads, percentageOfNetMismappedReads, snpReport2);
			}
			else
			{
				populateAlleleQuantitiesLookup(contigName, contigLength, snpPos, refAlleleAtLast, accumBaseA, accumBaseC, accumBaseT, accumBaseG, accumUndetermined, acummRefAllele, totalDiffAlleles, numberOfOverlappingReads, presenceOf3OrMoreAlleles, presenceOfBothAlleles, presenceOnlyAlternateAlleles, presenceOnlyReferenceAllele, presenceOfOddOccurrence, snpReport);
			}
			// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			accumBaseA = 0;
			accumBaseC = 0;
			accumBaseT = 0;
			accumBaseG = 0;
			accumUndetermined = 0;
			acummRefAllele = 0;
			totalDiffAlleles = 0;
			numberOfOverlappingReads = 0;
			counterOccurrenceAlleles = 0;
			presenceOfBothAlleles = false;
			presenceOnlyAlternateAlleles = false;
			presenceOnlyReferenceAllele = false;
			presenceOfOddOccurrence = false;
			presenceOf3OrMoreAlleles = false;
			for (int i = 1; i == 5; i++)
			{
				testPresenceFor3OrMoreAlleles[i] = 0;
			}
			// ----------------------------------------- New block of code created specifically for the "RegionQuantifier" class ------------------------------------------------------------------------------------------------------------------------------------------
			if (databaseFile != null)
			{
				isRevComp = false;
				dist2SNPFromEndOfContig = 0;
				snpPosOnGenome = 0;
				misAssembly = false;
				chrmInGenome = "";
				chrmID = "";
				lowerBound = 0;
				upperBound = 0;
				accumFragFromDiffChrm = 0;
				accumFragFromSameChrm = 0;
				accumFragFromWrongRegionSameChrm = 0;
				accumMismappedRead = 0;
				accumCorrectlyMappedRead = 0;
				accumMismatchedRead = 0;
				accumMatchedRead = 0;
				refAlleleOnGenome4Comparison = null;
				refAlleleAtPos4Comparison = null;
				refAlleleOnGenome = ' ';
				complementedRefAlleleOnGenome = ' ';
				accumMisAssemblyMismatchedRead = 0;
				accumMisAssemblyMatchedRead = 0;
				accumMisAssemblyDiffChrmMismapMatchedRead = 0;
				accumMisAssemblyDiffChrmMismapMismatchedRead = 0;
				accumMisAssemblySameChrmCorrectlyMappedMatchedRead = 0;
				accumMisAssemblySameChrmCorrectlyMappedMismatchedRead = 0;
				accumMisAssemblySameChrmMismapMatchedRead = 0;
				accumMisAssemblySameChrmMismapMismatchedRead = 0;
				accumCorrectAssyMismatchedRead = 0;
				accumCorrectAssyMatchedRead = 0;
				accumCorrectAssyDiffChrmMismapMatchedRead = 0;
				accumCorrectAssyDiffChrmMismapMismatchedRead = 0;
				accumCorrectAssySameChrmCorrectlyMappedMatchedRead = 0;
				accumCorrectAssySameChrmCorrectlyMappedMismatchedRead = 0;
				accumCorrectAssySameChrmMismapMatchedRead = 0;
				accumCorrectAssySameChrmMismapMismatchedRead = 0;
				accumOddSituation4MisAssemblyScenario = 0;
				accumOddSituation4CorrectAssyScenario = 0;
				accumCasesOfCorrectAssyButMismapMismatch = 0;
				accumCasesOfMisAssemblyAndMismapMatch = 0;
				accumNetMismappedReads = 0;
				percentageOfNetMismappedReads = 0;
				sStart = 0;
				sEnd = 0;
				blastResultIsEmpty = false;
				nonACTGEvent = false;
				fpSNPSiteNotVisitedByBlast = false;
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// Tip extracted from kindly available code at: https://www.biostars.org/p/17842/
	private static String getCIGARCorrectedRead(SAMRecord rec, IndexedFastaSequenceFile faidx)
	{
		// System.out.println("getCIGARCorrectedRead");
		byte ref_sequence[] = faidx.getSubsequenceAt(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd()).getBases();
		byte read_sequence[] = rec.getReadBases();
		int ref_pos = 0;
		int read_pos = 0;
		Cigar cigar = rec.getCigar();
		StringBuilder refseq = new StringBuilder();
		StringBuilder readseq = new StringBuilder();
		// System.out.println("CigarElements:");
		for (CigarElement c : cigar.getCigarElements())
		{
			// System.out.println(c.getOperator()+""+c.getLength());
			switch (c.getOperator())
			{
				case EQ:
				case M:
				{
					for (int i = 0; i < c.getLength(); ++i)
					{
						refseq.append((char) ref_sequence[ref_pos++]);
						readseq.append((char) read_sequence[read_pos++]);
					}
					break;
				}
				case I:
				case S:
				{
					for (int i = 0; i < c.getLength(); ++i)
					{
						refseq.append("-");
						char ch = (char) read_sequence[read_pos++];
						readseq.append(ch);
						// System.out.println("appending to read " + ch );
					}
					break;
				}
				case P:
					break;
				case H:
					break;
				case D:
				{
					for (int i = 0; i < c.getLength(); ++i)
					{
						refseq.append((char) ref_sequence[ref_pos++]);
						readseq.append("-");
					}
					break;
				}
				case N:
				{
					for (int i = 0; i < c.getLength(); ++i)
					{
						readseq.append("-");
					}
					break;
				}
				default:
				{
					throw new IllegalArgumentException(c.getOperator() + "\n" + refseq + "\n" + readseq + "\n" + new String(ref_sequence) + "\n" + new String(read_sequence) + "\n" + cigar);
				}
			}
		}
		// convert the string buffer to a string
		String newReadSeq = readseq.toString();
		// check whether the first CIGAR element is a soft clip
		CigarElement el = cigar.getCigarElements().get(0);
		// if yes, trim a sequence of the corresponding length from the start of the read so that the coordinate system is in keeping with what we see in Tablet and for our SNP discovery
		if (el.getOperator().toString().equals("S"))
		{
			// System.out.println("soft clip found at read start");
			int length = el.getLength();
			newReadSeq = newReadSeq.substring(length, newReadSeq.length());
		}
		return newReadSeq;
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	public static void populateAlleleQuantitiesLookup(String transcriptName, int transcriptLength, int snpPos, char refAlleleAtLast, int accumBaseA, int accumBaseC, int accumBaseT, int accumBaseG, int accumUndetermined, int acummRefAllele, int totalDiffAlleles, int numberOfOverlappingReads, boolean presenceOf3OrMoreAlleles, boolean presenceOfBothAlleles, boolean presenceOnlyAlternateAlleles, boolean presenceOnlyReferenceAllele, boolean presenceOfOddOccurrence, SNPReport snpReport)
	{
		snpReport.transcriptName = transcriptName;
		snpReport.transcriptLength = transcriptLength;
		snpReport.snpPos = snpPos;
		snpReport.refAlleleAtLast = refAlleleAtLast;
		snpReport.accumBaseA = accumBaseA;
		snpReport.accumBaseC = accumBaseC;
		snpReport.accumBaseT = accumBaseT;
		snpReport.accumBaseG = accumBaseG;
		snpReport.accumUndetermined = accumUndetermined;
		snpReport.acummRefAllele = acummRefAllele;
		snpReport.totalDiffAlleles = totalDiffAlleles;
		snpReport.numberOfOverlappingReads = numberOfOverlappingReads;
		snpReport.presenceOf3OrMoreAlleles = presenceOf3OrMoreAlleles;
		snpReport.presenceOfBothAlleles = presenceOfBothAlleles;
		snpReport.presenceOnlyAlternateAlleles = presenceOnlyAlternateAlleles;
		snpReport.presenceOnlyReferenceAllele = presenceOnlyReferenceAllele;
		snpReport.presenceOfOddOccurrence = presenceOfOddOccurrence;
		// String transcriptAndpositionLiteral = transcriptName+"_"+String.valueOf(snpPos);
		// alleleQuantitiesLookup.put(transcriptName, contentsRelated2AlleleQuantities);
		// infoAboutRefAllele.put(transcriptAndpositionLiteral, refAlleleAtLast);
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	public static void populateAlleleQuantitiesLookup2(String contigName, int contigLength, int snpPos, char refAlleleAtLast, int accumBaseA, int accumBaseC, int accumBaseT, int accumBaseG, int accumUndetermined, int acummRefAllele, int totalDiffAlleles, int numberOfOverlappingReads, boolean presenceOf3OrMoreAlleles, boolean presenceOfBothAlleles, boolean presenceOnlyAlternateAlleles, boolean presenceOnlyReferenceAllele, boolean presenceOfOddOccurrence, int sStart, int sEnd, boolean blastResultIsEmpty, boolean nonACTGEvent, boolean fpSNPSiteNotVisitedByBlast, boolean isRevComp, int dist2SNPFromEndOfContig, int snpPosOnGenome, boolean misAssembly, String chrmInGenome, int accumFragFromDiffChrm, int accumFragFromSameChrm, int accumFragFromWrongRegionSameChrm, int accumMismappedRead, int accumCorrectlyMappedRead, int accumMismatchedRead, int accumMatchedRead, Character refAlleleOnGenome4Comparison, int accumMisAssemblyMismatchedRead, int accumMisAssemblyMatchedRead, int accumMisAssemblyDiffChrmMismapMatchedRead, int accumMisAssemblyDiffChrmMismapMismatchedRead, int accumMisAssemblySameChrmCorrectlyMappedMatchedRead, int accumMisAssemblySameChrmCorrectlyMappedMismatchedRead, int accumMisAssemblySameChrmMismapMatchedRead, int accumMisAssemblySameChrmMismapMismatchedRead, int accumCorrectAssyMismatchedRead, int accumCorrectAssyMatchedRead, int accumCorrectAssyDiffChrmMismapMatchedRead, int accumCorrectAssyDiffChrmMismapMismatchedRead, int accumCorrectAssySameChrmCorrectlyMappedMatchedRead, int accumCorrectAssySameChrmCorrectlyMappedMismatchedRead, int accumCorrectAssySameChrmMismapMatchedRead, int accumCorrectAssySameChrmMismapMismatchedRead, int accumOddSituation4MisAssemblyScenario, int accumOddSituation4CorrectAssyScenario, int accumNetMismappedReads, float percentageOfNetMismappedReads, SNPReport2 snpReport2)
	{
		snpReport2.contigName = contigName;
		snpReport2.contigLength = contigLength;
		snpReport2.snpPos = snpPos;
		snpReport2.refAlleleAtLast = refAlleleAtLast;
		snpReport2.accumBaseA = accumBaseA;
		snpReport2.accumBaseC = accumBaseC;
		snpReport2.accumBaseT = accumBaseT;
		snpReport2.accumBaseG = accumBaseG;
		snpReport2.accumUndetermined = accumUndetermined;
		snpReport2.acummRefAllele = acummRefAllele;
		snpReport2.totalDiffAlleles = totalDiffAlleles;
		snpReport2.numberOfOverlappingReads = numberOfOverlappingReads;
		snpReport2.presenceOf3OrMoreAlleles = presenceOf3OrMoreAlleles;
		snpReport2.presenceOfBothAlleles = presenceOfBothAlleles;
		snpReport2.presenceOnlyAlternateAlleles = presenceOnlyAlternateAlleles;
		snpReport2.presenceOnlyReferenceAllele = presenceOnlyReferenceAllele;
		snpReport2.presenceOfOddOccurrence = presenceOfOddOccurrence;
		// String transcriptAndpositionLiteral = transcriptName+"_"+String.valueOf(snpPos);
		// alleleQuantitiesLookup.put(transcriptName, contentsRelated2AlleleQuantities);
		// infoAboutRefAllele.put(transcriptAndpositionLiteral, refAlleleAtLast);
		// ----------------------------------------- New block of code created specifically for the "RegionQuantifier" class
		// ------------------------------------------------------------------------------------------------------------------------------------------
		// snpReport2.chrmID = chrmID;
		snpReport2.sStart = sStart;
		snpReport2.sEnd = sEnd;
		snpReport2.blastResultIsEmpty = blastResultIsEmpty;
		snpReport2.nonACTGEvent = nonACTGEvent;
		snpReport2.fpSNPSiteNotVisitedByBlast = fpSNPSiteNotVisitedByBlast;
		snpReport2.isRevComp = isRevComp;
		snpReport2.dist2SNPFromEndOfContig = dist2SNPFromEndOfContig;
		snpReport2.snpPosOnGenome = snpPosOnGenome;
		snpReport2.misAssembly = misAssembly;
		snpReport2.chrmInGenome = chrmInGenome;
		snpReport2.accumFragFromDiffChrm = accumFragFromDiffChrm;
		snpReport2.accumFragFromSameChrm = accumFragFromSameChrm;
		snpReport2.accumFragFromWrongRegionSameChrm = accumFragFromWrongRegionSameChrm;
		snpReport2.accumMismappedRead = accumMismappedRead;
		snpReport2.accumCorrectlyMappedRead = accumCorrectlyMappedRead;
		snpReport2.accumMismatchedRead = accumMismatchedRead;
		snpReport2.accumMatchedRead = accumMatchedRead;
		snpReport2.refAlleleOnGenome4Comparison = refAlleleOnGenome4Comparison;
		snpReport2.accumMisAssemblyMismatchedRead = accumMisAssemblyMismatchedRead;
		snpReport2.accumMisAssemblyMatchedRead = accumMisAssemblyMatchedRead;
		snpReport2.accumMisAssemblyDiffChrmMismapMatchedRead = accumMisAssemblyDiffChrmMismapMatchedRead;
		snpReport2.accumMisAssemblyDiffChrmMismapMismatchedRead = accumMisAssemblyDiffChrmMismapMismatchedRead;
		snpReport2.accumMisAssemblySameChrmCorrectlyMappedMatchedRead = accumMisAssemblySameChrmCorrectlyMappedMatchedRead;
		snpReport2.accumMisAssemblySameChrmCorrectlyMappedMismatchedRead = accumMisAssemblySameChrmCorrectlyMappedMismatchedRead;
		snpReport2.accumMisAssemblySameChrmMismapMatchedRead = accumMisAssemblySameChrmMismapMatchedRead;
		snpReport2.accumMisAssemblySameChrmMismapMismatchedRead = accumMisAssemblySameChrmMismapMismatchedRead;
		snpReport2.accumCorrectAssyMismatchedRead = accumCorrectAssyMismatchedRead;
		snpReport2.accumCorrectAssyMatchedRead = accumCorrectAssyMatchedRead;
		snpReport2.accumCorrectAssyDiffChrmMismapMatchedRead = accumCorrectAssyDiffChrmMismapMatchedRead;
		snpReport2.accumCorrectAssyDiffChrmMismapMismatchedRead = accumCorrectAssyDiffChrmMismapMismatchedRead;
		snpReport2.accumCorrectAssySameChrmCorrectlyMappedMatchedRead = accumCorrectAssySameChrmCorrectlyMappedMatchedRead;
		snpReport2.accumCorrectAssySameChrmCorrectlyMappedMismatchedRead = accumCorrectAssySameChrmCorrectlyMappedMismatchedRead;
		snpReport2.accumCorrectAssySameChrmMismapMatchedRead = accumCorrectAssySameChrmMismapMatchedRead;
		snpReport2.accumCorrectAssySameChrmMismapMismatchedRead = accumCorrectAssySameChrmMismapMismatchedRead;
		snpReport2.accumOddSituation4MisAssemblyScenario = accumOddSituation4MisAssemblyScenario;
		snpReport2.accumOddSituation4CorrectAssyScenario = accumOddSituation4CorrectAssyScenario;
		snpReport2.accumNetMismappedReads = accumNetMismappedReads;
		snpReport2.percentageOfNetMismappedReads = percentageOfNetMismappedReads;
		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	private static File samReHeader(String contigName, File greppedSAMFile) throws IOException
	{
		// Manipulation of the sam file resulted from "samtools view -h -F 4 <sortedBAMFile> > contigName.sam to keep only the @SQ line for that contigName
		// naming the file with the sequence name and adding the .fasta extension
		File samReHeadedFile = new File(contigName + "_reHeaded.sam");
		// opening an output Stream to that file
		BufferedWriter bw = new BufferedWriter(new FileWriter(samReHeadedFile));
		Integer lineCount = 0;
		System.out.println();
		System.out.println("\nReading the produced SAM file of the mapping for the given contig and fixing its header...");
		try
		{
			// reading the contigName SAM file
			// System.out.println("\nTrying to open the file...");
			FileInputStream fstream = new FileInputStream(greppedSAMFile);
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String line = null;
			while ((line = br.readLine()) != null)
			{
				// System.out.println(line);
				if (!line.startsWith("@SQ"))
				{
					bw.write(line);
					lineCount++;
					// newline method putting the plataform specific new line char
					bw.newLine();
				}
				if (line.startsWith("@SQ"))
				{
					String[] lineArray = line.split("\t");
					String str = lineArray[1];
					// if (debug) System.out.println("str: "+str);
					String str2Compare = "SN:" + contigName;
					if (str.equalsIgnoreCase(str2Compare))
					{
						bw.write(line);
						lineCount++;
						// newline method putting the plataform specific new line char
						bw.newLine();
					}
				}
			}
			in.close();
			bw.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		System.out.println("\ndone -- File " + contigName + "_reHeaded.sam rebuilt with " + lineCount + " lines.");
		System.out.println();
		lineCount = 0;
		return samReHeadedFile;
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	private static String stripExtension(String str)
	{
		// Tip extracted from kindly available code at: http://stackoverflow.com/questions/924394/how-to-get-file-name-without-the-extension
		// Handling null case...
		if (str == null)
			return null;
		// Get position of the last dot ...
		int pos = str.lastIndexOf(".");
		// if none dot is present, just return the string as it is...
		if (pos == -1)
			return str;
		// otherwise, return the string up to the dot...
		return str.substring(0, pos);
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	private static char complementAllele(char c)
	{
		char complementedAllele = ' ';
		if ((c == 'A') || (c == 'a'))
		{
			complementedAllele = 'T';
		}
		else if ((c == 'C') || (c == 'c'))
		{
			complementedAllele = 'G';
		}
		else if ((c == 'T') || (c == 't'))
		{
			complementedAllele = 'A';
		}
		else if ((c == 'U') || (c == 'u'))
		{
			complementedAllele = 'A';
		}
		else if ((c == 'G') || (c == 'g'))
		{
			complementedAllele = 'C';
		}
		else if ((c == 'Y') || (c == 'y'))
		{
			complementedAllele = 'R';
		}
		else if ((c == 'R') || (c == 'r'))
		{
			complementedAllele = 'Y';
		}
		else if ((c == 'S') || (c == 's'))
		{
			complementedAllele = 'S';
		}
		else if ((c == 'W') || (c == 'w'))
		{
			complementedAllele = 'W';
		}
		else if ((c == 'K') || (c == 'k'))
		{
			complementedAllele = 'M';
		}
		else if ((c == 'M') || (c == 'm'))
		{
			complementedAllele = 'K';
		}
		else if ((c == 'B') || (c == 'b'))
		{
			complementedAllele = 'V';
		}
		else if ((c == 'V') || (c == 'v'))
		{
			complementedAllele = 'B';
		}
		else if ((c == 'D') || (c == 'd'))
		{
			complementedAllele = 'H';
		}
		else if ((c == 'H') || (c == 'h'))
		{
			complementedAllele = 'D';
		}
		else if ((c == 'N') || (c == 'n'))
		{
			complementedAllele = 'N';
		}
		return complementedAllele;
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
}
