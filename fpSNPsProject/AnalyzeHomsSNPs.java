package fps;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
//import java.nio.file.*;
import java.util.*;
import utils.entities.Read;
import utils.fasta.*;
import utils.validation.transcripts.*;
import utilsAntonio.entities.*;
import net.sf.picard.reference.*;
//import net.sf.samtools.*;
//import java.math.*;

public class AnalyzeHomsSNPs
{
	static final String USAGE = "java fps.AnalyzeHomsSNPs <FASTA file> <Truncate name at space character? true | false> <List of Homozygous SNPs file> <Path to the assembly directory> <AFG file produced by Velvet> <Sequences file produced by Velvet> <Number of mismatches for Bowtie2 mapping> <Number of Bowtie2 running threads> <outputFile> <Run BLAST against any specific database? true | false> <Target BLAST database file | none>";
	public static int totalUsedReads = 0;
	public static HashMap<String, ArrayList<Integer>> usedReadsLookup = new HashMap<String, ArrayList<Integer>>();
	//this keeps a tab of all our SNPReport objects
	static ArrayList<SNPReport> snpReports = new ArrayList<SNPReport>();
	public static HashMap<Integer, Float> bitScoreLookup = new HashMap<Integer, Float>();
	public static String databaseFile = null;
	//Debugging switch block ----------------------------------------------
	static boolean debug = true;
	//static boolean debug = false;
	//---------------------------------------------------------------------
	public static void main(String[] args)
	{
		if (args.length != 11) {
			System.out.println(USAGE);
			System.exit(0);
		}
		try
		{
			boolean runBLAST = Boolean.parseBoolean(args[9]);
			if (runBLAST == true)
			{
				databaseFile = new String(args[10]);
			}
			File fastaFile = new File(args[0]);
			boolean truncateNameAtSpaceChar = Boolean.parseBoolean(args[1]);
			File listHomSNPsFile = new File(args[2]);
			String path2AssyDir = new String(args[3]);
			File velvetAFGFile = new File(args[4]);
			File velvetSequencesFile = new File(args[5]);
			int numberOfMismatches = Integer.parseInt(args[6]);
			int numberOfThreads = Integer.parseInt(args[7]);
			analyze(fastaFile, truncateNameAtSpaceChar, listHomSNPsFile, path2AssyDir, velvetAFGFile, velvetSequencesFile, numberOfMismatches, numberOfThreads, databaseFile);
			outputResults(new File(args[8]));
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
	
	private static void analyze(File fastaFile, boolean truncateNameAtSpaceChar, File listHomSNPsFile, String path2AssyDir, File velvetAFGFile, File velvetSequencesFile, int numberOfMismatches, int numberOfThreads, String databaseFile) throws IOException
	{
		float coefficient4Bowtie2 = calculateBowtie2Coefficient(numberOfMismatches);
		System.out.println("\ncoefficient4Bowtie2 = "+coefficient4Bowtie2);
		String currentDir = System.getProperty("user.dir");
		System.out.println("\ncurrentDir = "+currentDir);
		System.out.println("\nParsing a given assembly file and storing each entry...");
		// read all contigs into sequence lookup (key = contig name, value = sequence)
		LinkedHashMap<String,String> sequenceLookup = FastaParser.parseFile(fastaFile, truncateNameAtSpaceChar);
		// creates the list of entries NODE_10001_length_277_cov_91.635376_143 (contig name and SNP position) from a given text file in the format
		ArrayList<String> snpList = parseSNPList(listHomSNPsFile);
		// make lookup of contigs (keys) versus lists of SNPs within these (values)
		// Method used for old work with transcripts: HashMap<String, TreeSet<Integer>> snpLookup = makeSNPLookup(snpList);
		// New method call included here to deal with Velvet assembled contigs...
		// ----------------------------------------------------------------------
		HashMap<String, TreeSet<Integer>> snpLookup = makeSNPLookup2(snpList);
		// ----------------------------------------------------------------------
		// File velvetSequencesFileConverted2Fasta = convertVelvetSequencesFile2Fasta(velvetSequencesFile);
		// read all Velvet Sequences into velvetSequencesLookup (key = readHeader, value = sequence)
		// LinkedHashMap<String,String> velvetSequencesLookup = FastaParser.parseFile(velvetSequencesFileConverted2Fasta, truncateNameAtSpaceChar);
		HashMap<Integer,VelvetSeqsFileEntry> velvetSeqsLookup = makeVelvetSeqsFileLookup(velvetSequencesFile);
		// ----------------------------------------------------------------------
		// Temporary validation method for makeVelvetSeqsFileLookup
		testVelvetSeqsLookup(velvetSeqsLookup);
		// ----------------------------------------------------------------------
		// for each contig in keyset of lookup
		for (String contigName : snpLookup.keySet())
		{
			// look up its sequence in the sequence lookup
			String sequence = sequenceLookup.get(contigName);
			// output it in FASTA format
			File contigFastaFile = outputSeq2Fasta(contigName,sequence);
			// index this file for Bowtie2 and by SAMTOOLS
			System.out.println("\nProcessing the file: "+contigFastaFile.getAbsolutePath());
			String com4AlignerIndex = "bowtie2-build "+contigFastaFile.getAbsolutePath()+" "+contigName;
			// If BWA: String com4BWAIndex = "bwa index "+contigFastaFile.getAbsolutePath();
			RunJob.runJob(com4AlignerIndex);
			// index this file by SAMTOOLS
			String com4Faidx = "samtools faidx "+contigFastaFile.getAbsolutePath();
			RunJob.runJob(com4Faidx);
			// get the proper reads file from the assembly output (preferably already in FASTA format...)
			System.out.println("\n\nTrying to retrieve the used reads in the original contig assembly...");
			// retrieving the contig number
			String contigNumber = contigName.split("_")[1];
			System.out.println("\ncontigNumber = "+contigNumber);
			// Firstly, trying to get an splitted assembly file related to the contig... 
			String com4runningVelvetContribAsblySplitter = "/home/ar41690/antonio_bio_software/velvet_1.2.10/contrib/afg_handling/asmbly_splitter.pl "+contigNumber+" "+velvetAFGFile;
			System.out.println("\ncom4runningVelvetContribAsblySplitter = "+com4runningVelvetContribAsblySplitter);
			RunJob.runJob(com4runningVelvetContribAsblySplitter);
			String velvetAFGFilenameWOExtension = stripExtension(velvetAFGFile.getName());
			System.out.println("\nvelvetAFGFilenameWOExtension = "+velvetAFGFilenameWOExtension);
			String com4movingProducedAFGSplittedFile = "mv "+path2AssyDir+"/"+velvetAFGFilenameWOExtension+"_"+contigNumber+".afg .";
			System.out.println("\ncom4movingProducedAFGSplittedFile = "+com4movingProducedAFGSplittedFile);
			RunJob.runJob(com4movingProducedAFGSplittedFile);
			File splittedAFGFile = new File(velvetAFGFilenameWOExtension+"_"+contigNumber+".afg");
			// Trying to get the reads that took part in the contig assembly based on the retrieved portion of the original AFG file produced by Velvet...
			// Firstly, trying to parse the given AFG file in order to retrieve the RED events...
			ArrayList<Integer> REDentries4contig = parseAFGFile4REDentries(splittedAFGFile);
			// Once there is a list of RED entries, try to find (map) the corresponding reads in the Velvet Sequences file and produce a new file with the subset of reads for this contig...
			extractRawReadsFromSequencesFile(contigName,velvetSeqsLookup,REDentries4contig);
			// -------------------------------------------------------------------------------
			// NOT USED/NEEDED IN THIS NEW CODE UNTIL FURTHER NOTICE: File handledRawCompReadsFile = new File(contigName+"_rawReads.fa");
			// map all reads with Bowtie2 to the contig... All this block will be eventually refactored later in order to call a specific helper method for each task.........
			System.out.println("\nExecuting the mappings of each set of reads related to the respective contig "+contigName+" assembly event...");
			// aligning the reads
			// if BWA: System.out.println("\nGenerating the required file: "+contigName+".sai");
			// if BWA: String com4BWAAlignment = "bwa aln -n "+numberOfMismatches+" -l 100 -t "+numberOfThreads+" -o 0 -f "+contigName+".sai "+contigFastaFile.getAbsolutePath() +" "+handledRawCompReadsFile.getAbsolutePath();
			// if BWA: RunJob.runJob(com4BWAAlignment);
			// if BWA: File resultFromBWAAlignment = new File(contigName+".sai");
			// if BWA: String com4BWASAM = "bwa samse -f "+contigName+".sam "+contigFastaFile.getAbsolutePath()+" "+resultFromBWAAlignment.getAbsolutePath()+" "+handledRawCompReadsFile.getAbsolutePath();
			// if BWA: RunJob.runJob(com4BWASAM);
			// if BWA: File bwaSAMFile = new File(contigName+".sam");
			// Solution to get the basename of the file, but not used here...
			// ------------------------------------------------------------------------------
			//String[] filenameTokens = contigFastaFile.getName().split("\\.(?=[^\\.]+$)");
			//String contigFastaFileBasename = filenameTokens[0];
			// if (debug) System.out.println("\ncontigFastaFileBasename = "+contigFastaFileBasename);
			// ------------------------------------------------------------------------------
			String com4Bowtie2Alignment = "bowtie2 -f -x "+contigName+" -U "+contigName+"_rawReads.fa -S "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sam --score-min L,0,"+coefficient4Bowtie2+" -p "+numberOfThreads+" --un-gz ./"+contigName+"_rawReadsAligned_mis"+numberOfMismatches+"_unal.sam.gz --al-gz ./"+contigName+"_rawReadsAligned_mis"+numberOfMismatches+"_once.sam.gz  --un-conc-gz ./"+contigName+"_rawReadsAligned_mis"+numberOfMismatches+"_unconc.sam.gz --al-conc-gz ./"+contigName+"_rawReadsAligned_mis"+numberOfMismatches+"_conc.sam.gz --rg-id "+contigName+"_mis_"+numberOfMismatches+" --rg SM:Pool1";
			System.out.println("\ncom4Bowtie2Alignment = "+com4Bowtie2Alignment);
			RunJob.runJob(com4Bowtie2Alignment);
			File bowtie2SAMFile = new File(contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sam");
			// converting SAM2BAM
			// if BWA: if (debug) System.out.println("\nCreating the required file: "+contigName+".unsorted.bam");
			// if BWA: String comSAM2BAM = "samtools view -S -b -o "+contigName+".unsorted.bam "+bwaSAMFile.getAbsolutePath();
			System.out.println("\nCreating the required file: "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".unsorted.bam");
			String comSAM2BAM = "samtools view -S -b -o "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".unsorted.bam "+bowtie2SAMFile.getAbsolutePath();
			RunJob.runJob(comSAM2BAM);
			// sorting BAM
			// if BWA: File convBAMFile = new File(contigName+".unsorted.bam");
			// if BWA: if (debug) System.out.println("\nCreating the required file: "+contigName+".sorted.bam");
			// if BWA: String comSortBAM = "samtools sort "+convBAMFile.getAbsolutePath()+" "+contigName+".sorted";
			File convBAMFile = new File(contigName+"_rawReadsAligned_mis"+numberOfMismatches+".unsorted.bam");
			System.out.println("\nCreating the required file: "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam");
			String comSortBAM = "samtools sort "+convBAMFile.getAbsolutePath()+" "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted";
			RunJob.runJob(comSortBAM);
			// indexing BAM
			// if BWA: if (debug) System.out.println("\nCreating the required file: "+contigName+".sorted.bam.bai");
			// if BWA: File sortedBAMFile = new File(contigName+".sorted.bam");
			System.out.println("\nCreating the required file: "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam.bai");
			File sortedBAMFile = new File(contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam");
			String comIndexBAM = "samtools index "+sortedBAMFile.getAbsolutePath();
			RunJob.runJob(comIndexBAM);
			// if BWA: File indexedBAMFile = new File(contigName+".sorted.bam.bai");
			// File indexedBAMFile = new File(contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam.bai");
			//System.out.println("\n -------- CHECKPOINT -------- : Confirming that the totalUsedReads variable has value = 0 before being used for another contig event: "+ totalUsedReads);
			//System.out.println("\n -------- CHECKPOINT -------- : Confirming the values of the usedReadsLookup HASHMAP....");
			//System.out.println("\n -------- CHECKPOINT -------- : Getting the values for the given contigName "+contigName+" in usedReadsLookup HASHMAP: "+usedReadsLookup.get(contigName));
			// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			// Cleaning up stage (to be coded after the proper inventory of all created files....)
			// contigFastaFile.delete();
			// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			
			// for each SNP in the transcript
			for (Integer snpPos : snpLookup.get(contigName))
			{
				//a new report object for this SNP only
				System.out.println("\nEntering at the loop step for SNP position "+snpPos+".");
				SNPReport snpReport = new SNPReport();
				snpReports.add(snpReport);
				// run samtools flagstat to get the percentage/count of reads mapped
				// if BWA: if (debug) System.out.println("\nGetting the proper statistics related to the file: "+contigName+".sorted.bam");
				// if BWA: String comFlagstat = "samtools flagstat "+sortedBAMFile.getAbsolutePath()+" > "+contigName+".sorted.bam.flagstat.txt";
				System.out.println("\nGetting the proper statistics related to the file: "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam");
				String comFlagstat = "samtools flagstat "+sortedBAMFile.getAbsolutePath()+" > "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam.flagstat.txt";
				// if BWA: resultFlagstatFile(comFlagstat, contigName, snpReport);
				resultFlagstatFile(comFlagstat, contigName, snpReport, numberOfMismatches);
				// if BWA: ATTENTION BECAUSE IT WOULD NEED A MODIFICATION IN THE METHOD REGARDING THE NAMING OF THE FILE OR A NEW METHOD CALL
				// get the reference allele at the SNP position
				int finalPosOfContig = sequence.length();
				//System.out.println("\nThe contig "+contigName+" has a length of "+finalPosOfContig+" bases.");
				char refAlleleAtPos = sequence.charAt((snpPos-1));
				// ------------------------------------------------------------------------------- NOT USED ANYMORE BLOCK - kept here for reference -----------------------------------------------------------------------------------------------------------------------------
				//System.out.println("\nThe found reference allele for the contig "+contigName+" at position "+position+" is: "+refAlleleAtPos+".");
				// get all reads that intersect with the SNP position (it will be handled by AlleleExtractor class)
				// extract all the alleles at the SNP position (it will be handled by AlleleExtractor class)
				// count how many alternate alleles we have (it will be handled by AlleleExtractor class)
				// for each alternate allele (it will be handled by AlleleExtractor class)
				// record its count (it will be handled by AlleleExtractor class)
				// ------------------------------------------------------------------------------- END OF NOT USED ANYMORE BLOCK ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				IndexedFastaSequenceFile faidx = new IndexedFastaSequenceFile(new File(contigFastaFile.getAbsolutePath()));
				//extract allele counts
				AlleleExtractor.extractAlleles(new File(sortedBAMFile.getAbsolutePath()), snpPos.intValue(), faidx, contigName, refAlleleAtPos, finalPosOfContig, snpReport, databaseFile);
				// report/output this appropriately
				//if(debug)
					snpReport.outputToStdout();
			}
		}
		/*System.out.println("\ndone - Processed and stored each entry of the given contigs assembly file...");
		System.out.println("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =");
		System.out.println("\nNow entering in the stage 2 of the analysis...");
		System.out.println("= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =");
		System.out.println();
		System.out.println("\nMaking BLAST operations to compute the best alignments and associated bit scores from all the previous BLAST operations executed for each set of reads for each contig...");
		String com4OverallBLAST = "/home/ar41690/phd_libraries/scan.sh _OverlappingReads.fasta";
		RunJob.runJob(com4OverallBLAST);
		String com4appendingBLASTResults = "/home/ar41690/phd_libraries/appendingBlasts.sh _OverReads_vs_BLASTDB.txt";
		RunJob.runJob(com4appendingBLASTResults);
		System.out.println("\ndone - Processed and BLASTed all the necessary files for the best alignments calculation stage...");
		System.out.println("\nNow populating the best 'bit scores' lookup table...");
		File resultOfAllBlasts = new File("resultOfAllBLASTs.txt");
		// make lookup of whole lenghts of reads (keys) versus better BLAST bit scores obtained (values)
		makeBitScoreLookup(resultOfAllBlasts);*/
		System.out.println("\nGood job!!! Tool execution FINISHED!");
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static void testVelvetSeqsLookup(HashMap<Integer, VelvetSeqsFileEntry> velvetSeqsLookup)
	{
		System.out.println("\nTesting velvetSeqsLookup accuracy...");
		VelvetSeqsFileEntry tempVelvetSeqsFileEntry = new VelvetSeqsFileEntry();
		tempVelvetSeqsFileEntry = velvetSeqsLookup.get(1);
		System.out.println("tempVelvetSeqsFileEntry.t1 = "+tempVelvetSeqsFileEntry.t1);
		System.out.println("tempVelvetSeqsFileEntry.t2 = "+tempVelvetSeqsFileEntry.t2);
		System.out.println("tempVelvetSeqsFileEntry.readNameWithGtSign = "+tempVelvetSeqsFileEntry.readNameWithGtSign);
		System.out.println("tempVelvetSeqsFileEntry.sequence = "+tempVelvetSeqsFileEntry.sequence);
		tempVelvetSeqsFileEntry = velvetSeqsLookup.get(2);
		System.out.println("tempVelvetSeqsFileEntry.t1 = "+tempVelvetSeqsFileEntry.t1);
		System.out.println("tempVelvetSeqsFileEntry.t2 = "+tempVelvetSeqsFileEntry.t2);
		System.out.println("tempVelvetSeqsFileEntry.readNameWithGtSign = "+tempVelvetSeqsFileEntry.readNameWithGtSign);
		System.out.println("tempVelvetSeqsFileEntry.sequence = "+tempVelvetSeqsFileEntry.sequence);
		tempVelvetSeqsFileEntry = velvetSeqsLookup.get(3);
		System.out.println("tempVelvetSeqsFileEntry.t1 = "+tempVelvetSeqsFileEntry.t1);
		System.out.println("tempVelvetSeqsFileEntry.t2 = "+tempVelvetSeqsFileEntry.t2);
		System.out.println("tempVelvetSeqsFileEntry.readNameWithGtSign = "+tempVelvetSeqsFileEntry.readNameWithGtSign);
		System.out.println("tempVelvetSeqsFileEntry.sequence = "+tempVelvetSeqsFileEntry.sequence);
		tempVelvetSeqsFileEntry = velvetSeqsLookup.get(4);
		System.out.println("tempVelvetSeqsFileEntry.t1 = "+tempVelvetSeqsFileEntry.t1);
		System.out.println("tempVelvetSeqsFileEntry.t2 = "+tempVelvetSeqsFileEntry.t2);
		System.out.println("tempVelvetSeqsFileEntry.readNameWithGtSign = "+tempVelvetSeqsFileEntry.readNameWithGtSign);
		System.out.println("tempVelvetSeqsFileEntry.sequence = "+tempVelvetSeqsFileEntry.sequence);
		tempVelvetSeqsFileEntry = velvetSeqsLookup.get(5);
		System.out.println("tempVelvetSeqsFileEntry.t1 = "+tempVelvetSeqsFileEntry.t1);
		System.out.println("tempVelvetSeqsFileEntry.t2 = "+tempVelvetSeqsFileEntry.t2);
		System.out.println("tempVelvetSeqsFileEntry.readNameWithGtSign = "+tempVelvetSeqsFileEntry.readNameWithGtSign);
		System.out.println("tempVelvetSeqsFileEntry.sequence = "+tempVelvetSeqsFileEntry.sequence);
		System.out.println("\ndone testing VelvetSeqsLookup accuracy...");
	}

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static HashMap<Integer,VelvetSeqsFileEntry> makeVelvetSeqsFileLookup(File velvetSequencesFile)
	{
		System.out.println("\nConverting the original Velvet Sequences file into a Lookup data structure...");
		// Manipulation of the Velvet Sequences file to convert it to the Lookup data structure using the iids as keys and the whole combination of a Sequence entry as values
		HashMap<Integer, VelvetSeqsFileEntry> tempVelvetSeqsLookup = new HashMap<Integer, VelvetSeqsFileEntry>();
		// System.out.println("\nReading the Velvet Sequences file...");
		VelvetSeqsFileEntry tempVelvetSeqsFileEntry = new VelvetSeqsFileEntry();
		try
		{
			// System.out.println("\nTrying to open the file...");
			BufferedReader br = new BufferedReader(new FileReader(velvetSequencesFile));
			String line = null;
			String firstPortion = null;
			String secondPortion = null;
			String thirdPortion = null;
			String sequenceAtLast = null;
			while ((line = br.readLine()) != null)
			{
				// System.out.println(line);
				if (line.startsWith(">"))
				{
					String[] lineArray = line.split("\t");
					String gtSign_plus_readID = lineArray[0];
					int t1 = Integer.parseInt(lineArray[1]);
					int t2 = Integer.parseInt(lineArray[2]);
					//if (debug) System.out.println("\nPosition at lineArray[0] which would represent gtSign_plus_readID: "+gtSign_plus_readID);
					//if (debug) System.out.println("\nPosition at lineArray[1] which would represent t1: "+t1);
					//if (debug) System.out.println("\nPosition at lineArray[2] which would represent t2: "+t2);
					firstPortion = br.readLine();
					//if (debug) System.out.println("firstPortion = "+firstPortion);
					secondPortion = br.readLine();
					//if (debug) System.out.println("secondPortion = "+secondPortion);
					thirdPortion = br.readLine();
					//if (debug) System.out.println("thirdPortion = "+thirdPortion);
					sequenceAtLast = firstPortion+secondPortion+thirdPortion;
					//if (debug) System.out.println("sequenceAtLast = "+sequenceAtLast);
					// ------------- inserting the tokens into the data structure
					//if the t1 index doesn't exist, creates a new entry for it and, also, for its respective tempVelvetSeqsFileEntry...
					tempVelvetSeqsFileEntry = tempVelvetSeqsLookup.get(t1);
					if (tempVelvetSeqsFileEntry == null)
					{
						// System.out.println("I'm entering the branch where, supposedly, tempVelvetSeqsFileEntry == null and I should include a new entry in the tempVelvetSeqsLookup data structure...");
						// new instance of snpPositions object
						tempVelvetSeqsFileEntry = new VelvetSeqsFileEntry();
						tempVelvetSeqsFileEntry.t1 = t1;
						//if (debug) System.out.println("\nWhat I'm assigning to tempVelvetSeqsFileEntry.t1 = "+t1); 
						tempVelvetSeqsFileEntry.t2 = t2;
						//if (debug) System.out.println("\nWhat I'm assigning to tempVelvetSeqsFileEntry.t2 = "+t2);
						tempVelvetSeqsFileEntry.readNameWithGtSign = gtSign_plus_readID.trim();
						//if (debug) System.out.println("\nWhat I'm assigning to tempVelvetSeqsFileEntry.readNameWithGtSign = "+gtSign_plus_readID);
						tempVelvetSeqsFileEntry.sequence = sequenceAtLast.trim();
						//if (debug) System.out.println("\nWhat I'm assigning to tempVelvetSeqsFileEntry.sequence = "+sequenceAtLast);
						// putting the index t1 and respective tempVelvetSeqsFileEntry
						tempVelvetSeqsLookup.put(t1, tempVelvetSeqsFileEntry);
					}
				}
			}
			br.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		System.out.println("\ndone -- tempVelvetSeqsLookup structure created with "+tempVelvetSeqsLookup.size()+" entries.");
		//PopulateTotalUsedReads(totalUsedReads);
		return tempVelvetSeqsLookup;
	}
		
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static float calculateBowtie2Coefficient(int numberOfMismatches)
	{
		final int penalty = -6; // obtained from Bowtie2 manual
		final float readLength = 150f; // I'm just Hard-wiring it for now for the sake of this specific study...
		int recalculatedPenalty = numberOfMismatches * penalty;
		float coeff = recalculatedPenalty/readLength;
		return coeff;
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
	
	// make lookup of contigs (keys) versus lists of SNPs within these (values) from a given snpList produced by the parseSNPList method
	private static HashMap<String,TreeSet<Integer>> makeSNPLookup(ArrayList<String> snpList)
	{
		System.out.println("\nMaking SNP Lookup 'data structure'...");
		// Stores contigs names as keys and lists of SNPs positions as values
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
		
	// make lookup of contigs (keys) versus lists of SNPs within these (values) from a given snpList produced by the parseSNPList method
	private static HashMap<String,TreeSet<Integer>> makeSNPLookup2(ArrayList<String> snpList)
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
	
	private static File outputSeq2Fasta(String sequenceName,String sequence) throws IOException
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
	
	private static File convertVelvetSequencesFile2Fasta(File velvetSequencesFile) throws IOException
	{
		System.out.println("\nConverting the original Velvet Sequences file into a FASTA type one...");
		// Manipulation of the Velvet Sequences file to convert it to the FASTA format
		// renaming the file and adding the .fasta extension
		File velvetSequencesFileConverted2Fasta = new File("velvetSequencesConverted2.fasta");
		// opening an output Stream to that file
		BufferedWriter bw = new BufferedWriter(new FileWriter(velvetSequencesFileConverted2Fasta));
		Integer readCount = 0;
		// System.out.println("\nReading the Velvet Sequences file...");
		try
		{
			// System.out.println("\nTrying to open the file...");
			BufferedReader br = new BufferedReader(new FileReader(velvetSequencesFile));
			String line = null;
			String firstPortion = null;
			String secondPortion = null;
			String thirdPortion = null;
			String sequenceAtLast = null;
			while ((line = br.readLine()) != null)
			{
				// System.out.println(line);
				if (line.startsWith(">"))
				{
					String[] lineArray = line.split("\t");
					String gtSign_plus_readID = lineArray[0];
					String t1 = lineArray[1];
					String t2 = lineArray[2];
					// if (debug) System.out.println("\nPosition at lineArray[0] which would represent gtSign_plus_readID: "+gtSign_plus_readID);
					// if (debug) System.out.println("\nPosition at lineArray[1] which would represent t1: "+t1);
					// if (debug) System.out.println("\nPosition at lineArray[2] which would represent t2: "+t2);
					firstPortion = br.readLine();
					// if (debug) System.out.println("firstPortion = "+firstPortion);
					secondPortion = br.readLine();
					// if (debug) System.out.println("secondPortion = "+secondPortion);
					thirdPortion = br.readLine();
					// if (debug) System.out.println("thirdPortion = "+thirdPortion);
					sequenceAtLast = firstPortion+secondPortion+thirdPortion;
					// if (debug) System.out.println("sequenceAtLast = "+sequenceAtLast);
					// ------------- writing sequence name FASTA style
					bw.write(gtSign_plus_readID+"$"+t1+"%"+t2);
					// newline method putting the plataform specific new line char
					bw.newLine();
					// writes the the sequence
					bw.write(sequenceAtLast);
					// newline method putting the plataform specific new line char
					bw.newLine();
					readCount++;
				}
			}
			br.close();
			bw.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		System.out.println("\ndone -- File velvetSequencesConverted2.fasta created with "+readCount+" reads.");
		//PopulateTotalUsedReads(totalUsedReads);
		readCount = 0;
		return velvetSequencesFileConverted2Fasta;
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
	
	private static void resultFlagstatFile(String input4FlagstatCommand,String contigName,SNPReport snpReport,int numberOfMismatches) throws IOException
	{
		// if BWA: File flagstatFile = new File(contigName+".sorted.bam.flagstat.txt");
		File flagstatFile = new File(contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam.flagstat.txt");
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
					//System.out.println("\nRetrieving the TOTAL NUMBER OF AVAILABLE READS for the mapping related to the "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam" file...");
					// if BWA: the above line would have to be changed to use the proper .sorted.bam file name previously coded...
					//System.out.println("\n"+str);
					numTotalAvailReads = Integer.parseInt(str.substring(0, str.indexOf(" ")));
					//System.out.println("\n"+numTotalAvailReads);
				}
				if (numLines == 3)
				{
					//System.out.println("\nRetrieving the NUMBER OF MAPPED/USED READS in the mapping related to the "+contigName+"_rawReadsAligned_mis"+numberOfMismatches+".sorted.bam" file...");
					// if BWA: the above line would have to be changed to use the proper .sorted.bam file name previously coded...
					//System.out.println("\n"+str+"\n");
					numUsedReads = Integer.parseInt(str.substring(0, str.indexOf(" ")));
					//System.out.println("\n"+numUsedReads);
				}
				// populate lookup of contigs (keys) versus lists of TOTAL AVAILABLE READS IN THE MAPPING and TOTAL USED READS IN THE MAPPING (values)
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
	
	private static ArrayList<Integer> parseAFGFile4REDentries(File inputAFGFile)
	{
		System.out.println("\nReading a given input AFG file and populating a list of its RED entries...");
		// reads the given input AFG file to create a TreeSet (of Integers) for each RED entry of the file...
		// But, firstly, getting some information regarding the file...
		String str = null;
		int numReads = 0;
		int numContigs = 0;
		//int totalNumTiles = 0;
		ArrayList<Integer> listOfREDiids = new ArrayList<Integer>();
		try
		{
			System.out.println("\nparsing AFG file...");
			BufferedReader reader = new BufferedReader(new FileReader(inputAFGFile));
			//NumberFormat formatter = new DecimalFormat("###,###");
			//int count = 0;
			while ((str = reader.readLine()) != null)
			{
				//count++;
				//if(count % 10000 == 0)
					//System.out.print("\r# lines parsed = " + formatter.format(count));
				
				if (str.startsWith("{RED"))
				{
					numReads++;
					if(debug) System.out.println(" new read");
					String read_eid = null;
					int read_iid = 0;
					while ((str = reader.readLine()) != null && str.length() > 0 && !str.startsWith("}"))
					{
						// external ID AKA read name
						if (str.startsWith("eid"))
						{
							read_eid = str.substring(str.indexOf(":")+1).trim();
							if(debug) System.out.println("read_eid = "+read_eid);
						}
						// internal ID -- for cross referencing within AFG file
						else if (str.startsWith("iid"))
						{
							read_iid = Integer.parseInt(str.substring(str.indexOf(":")+1));
							if(debug) System.out.println("read_iid = "+read_iid);
							//now add this read to the TreeSet
							if(debug) System.out.println("adding to list the read "+read_iid);
							listOfREDiids.add(read_iid);
						}
					}
				}
				else
					if (str.startsWith("{CTG"))
				{
					numContigs++;
					//processContig();
				}			
			}
			reader.close();
			System.out.println("\n\ndone");
			System.out.println("total num reads found  = " + numReads);
			System.out.println("total num contigs found  = " + numContigs);
			//System.out.println("total num tiles found  = " + totalNumTiles);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return listOfREDiids;
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
	
	private static void extractRawReadsFromSequencesFile(String contigName,HashMap<Integer,VelvetSeqsFileEntry> velvetSeqsLookup,ArrayList<Integer> REDentries4contig) throws IOException
	{
		System.out.println("\nTrying to build the "+contigName+"_rawReads.fa file");
		int numOfCreatedReads = 0;
		File fastaFile = new File(contigName+"_rawReads.fa");
		// opening an output Stream to that file
		BufferedWriter bw = new BufferedWriter(new FileWriter(fastaFile));
		VelvetSeqsFileEntry tempVelvetSeqsFileEntry = new VelvetSeqsFileEntry();
		for (Integer iid : REDentries4contig)
		{
			if (debug) System.out.println("\niid inside REDentries4contig structure = "+iid);
			// gets velvetSeqsLookup entry to compose the FASTA file to be written:
			// Too slow block removed... Trying another approach... --------------------------------------------------------------
			// Too slow! for (Integer key : velvetSeqsLookup.keySet())
			// Too slow!{
				// Too slow! if (iid.equals(key)) //look up its sequence in the sequence lookup
				// Too slow!{
			// Too slow block removed... Trying another approach... --------------------------------------------------------------
			// if (debug) System.out.println("I've entered in the branch where, supposedly, iid "+iid+" = key "+key+"..."); 
			tempVelvetSeqsFileEntry = velvetSeqsLookup.get(iid);
			//output it in FASTA format
			// -------------------------------------------------------------------
			if (debug) System.out.println("tempVelvetSeqsFileEntry.readNameWithGtSign = "+tempVelvetSeqsFileEntry.readNameWithGtSign);
			if (debug) System.out.println("tempVelvetSeqsFileEntry.t1 = "+tempVelvetSeqsFileEntry.t1);
			if (debug) System.out.println("tempVelvetSeqsFileEntry.t2 = "+tempVelvetSeqsFileEntry.t2);
			if (debug) System.out.println("tempVelvetSeqsFileEntry.sequence = "+tempVelvetSeqsFileEntry.sequence);
			// ---------------------------------------------------------------------
			// extracting readName and iid...
			bw.write(tempVelvetSeqsFileEntry.readNameWithGtSign);
			// newline method putting the plataform specific new line char
			bw.newLine();
			// writes the sequence
			bw.write(tempVelvetSeqsFileEntry.sequence);
			// newline method putting the plataform specific new line char
			bw.newLine();
			numOfCreatedReads++;
			// Too slow!{	}
			// Too slow!{}
		}
		// closing the output Stream
		bw.close();
		System.out.println("\nThe number of reads present in the "+fastaFile.getName()+" is "+numOfCreatedReads+".");	
	}
	
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
	
}


