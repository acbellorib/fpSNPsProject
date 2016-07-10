package fps;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;
import utils.validation.transcripts.*;
import net.sf.picard.reference.*;
import net.sf.samtools.*;
import utils.entities.*;

public class AlleleExtractor
{
	// this nasty data structure contains the result
	// essentially this is a lookup of lookups where at the top level we are storing against the SNP name a lookup of (sampleName, SNPresults)
	HashMap<String, SNPSequence> snpSeqLookup;
	// structure which will be further used to collect the quantities of different alleles found at the SNP position
	public static HashMap<String, ArrayList<Integer>> alleleQuantitiesLookup = new HashMap<String, ArrayList<Integer>>();
	// NOT USED so far: public static HashMap<String, Character> infoAboutRefAllele = new HashMap<String, Character>();
	// structure which will keep a lookup of the "transcriptName + snpPos" (key) while storing the HashMap of the overlapping read and its respective allele at the snpPos
	public static HashMap<String, HashMap<String, Character>> infoAboutReadAndAllele = new HashMap<String, HashMap<String, Character>>();
	//public static HashMap<String, boolean[]> betweenBLASTStagesRecord = new HashMap<String, boolean[]>();
	//-------------------------------------------------- Debugging option switch BLOCK ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	static boolean debug = true;
	//static boolean debug = false;
	//-------------------------------------------------- END OF debugging option switch BLOCK ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
	public static void extractAlleles(File bamFile, int snpPos, IndexedFastaSequenceFile faidx, String transcriptName, char refAlleleAtLast, int transcriptLength, SNPReport snpReport, String databaseFile)
	{
		System.out.println("\nExtracting alleles STEP...");
		try
		{
			File baiFile = new File(bamFile.getAbsolutePath() + ".bai");
			SAMFileReader samReader = new SAMFileReader(bamFile, baiFile, true);
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
				//an iterator over all the reads that overlap the SNP position 
				SAMRecordIterator samRecordIterator = samReader.queryOverlapping(seqName, snpPos, snpPos);
				//System.out.println("exploring reads at position " + snpPos);
				System.out.println("\nThe following reads are overlapping the position "+snpPos+":");
				//System.out.println();
				while (samRecordIterator.hasNext())
				{
					SAMRecord record = samRecordIterator.next();
					String readName = record.getReadName();
					//System.out.println(readName);
					numberOfOverlappingReads++;
					// block to create a hashmap from the UNIQUE overlapping reads only
					String readSequence = record.getReadString();
					overlappingReadsLookup.put(readSequence, readName);
					// extract the allele at the SNP position
					int readStart = record.getAlignmentStart();
					int snpPosOnRead = snpPos - readStart + 1;
					//System.out.println("\noriginal read:");
					//System.out.println(record.getReadString());
					String cigarCorrectedRead = getCIGARCorrectedRead(record, faidx);		
					//System.out.println("cigar corrected read:");
					//System.out.println(cigarCorrectedRead);		
					//System.out.println("cigar string:");
					//System.out.println(record.getCigarString());	
					//System.out.println("snpPosOnRead = " + snpPosOnRead);
					//System.out.println("readStart = " + readStart);
					char allele = cigarCorrectedRead.charAt(snpPosOnRead - 1);
					System.out.println("read " + readName + " has allele " + allele);
					// block to create a hashmap of the reads and respective alleles at the SNP position
					readsAndAllelesLookup.put(readName, allele);
					infoAboutReadAndAllele.put(transcriptName+"_"+snpPos, readsAndAllelesLookup);
					// Trying to quantify the different alleles found in each read overlapping the SNP position
					if ((allele == 'A') || (allele == 'a'))
					{
						accumBaseA++;
					} else if ((allele == 'C') || (allele == 'c'))
					{
						accumBaseC++;
					} else if ((allele == 'T') || (allele == 't'))
					{
						accumBaseT++;
					} else if ((allele == 'G') || (allele == 'g'))
					{
						accumBaseG++;
					} else
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
			for (int i =0; i < testPresenceFor3OrMoreAlleles.length; i++)
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
			//System.out.println();
			if ((refAlleleAtLast == 'A') || (refAlleleAtLast == 'a'))  
			{
				acummRefAllele = accumBaseA;
				totalDiffAlleles = accumBaseC + accumBaseT + accumBaseG + accumUndetermined;
			} 
			else if  ((refAlleleAtLast == 'C') || (refAlleleAtLast == 'c'))  
			{
				acummRefAllele = accumBaseC;
				totalDiffAlleles = accumBaseA + accumBaseT + accumBaseG + accumUndetermined;
			}
			else if  ((refAlleleAtLast == 'T') || (refAlleleAtLast == 't'))  
			{
				acummRefAllele = accumBaseT;
				totalDiffAlleles = accumBaseA + accumBaseC + accumBaseG + accumUndetermined;
			}
			else if  ((refAlleleAtLast == 'G') || (refAlleleAtLast == 'g'))  
			{
				acummRefAllele = accumBaseG;
				totalDiffAlleles = accumBaseA + accumBaseC + accumBaseT + accumUndetermined;
			}
			else 
			{
				acummRefAllele = accumUndetermined;
				totalDiffAlleles = accumBaseA + accumBaseC + accumBaseT + accumBaseG;
			}
			// populate lookup of transcripts (keys) versus positions and respective lists of quantities of each found allele
			if ((acummRefAllele > 0) && (totalDiffAlleles > 0)) 
			{
				presenceOfBothAlleles = true;
			}
			else if ((acummRefAllele == 0) && (totalDiffAlleles > 0))
			{
				presenceOnlyAlternateAlleles = true;
			}
			else if ((acummRefAllele > 0) && (totalDiffAlleles == 0))
			{
				presenceOnlyReferenceAllele = true;
			}
			else
			{
				presenceOfOddOccurrence = true;
			}
			populateAlleleQuantitiesLookup(transcriptName, transcriptLength, snpPos,refAlleleAtLast, accumBaseA, accumBaseC, accumBaseT, accumBaseG, accumUndetermined, acummRefAllele, totalDiffAlleles, numberOfOverlappingReads, presenceOf3OrMoreAlleles, presenceOfBothAlleles, presenceOnlyAlternateAlleles, presenceOnlyReferenceAllele, presenceOfOddOccurrence, snpReport);
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
			if (databaseFile != null) {
				File dbFile = new File(databaseFile);
				File overlappingReadsFastaFile = createFastaFileFromOverlappingReads(transcriptName, snpPos, overlappingReadsLookup);
				System.out.println("\nBLASTing the reads associated with the "+transcriptName+"_"+snpPos+" STEP...");
				String com4IntermediateBLAST = "blastn -query "+overlappingReadsFastaFile.getAbsolutePath()+" -db "+dbFile.getAbsolutePath()+" -outfmt 6 -out "+transcriptName+"_"+snpPos+"_OverReads_vs_BLASTDB.txt";
				// original path to the AK file just for the sake of backup: /mnt/shared/projects/barley/201109_fullLengthCDNAs/HvuFLcDNA_rep.fa
				RunJob.runJob(com4IntermediateBLAST);
				File blastResults4TrascriptAndSNPposition = new File(transcriptName+"_"+snpPos+"_OverReads_vs_BLASTDB.txt");
				File bestBLASTResultsFilteredFile = createBestBLASTResultsFilteredFile(blastResults4TrascriptAndSNPposition, transcriptName, snpPos);
				populateBetweenBLASTStagesRecord(blastResults4TrascriptAndSNPposition, bestBLASTResultsFilteredFile, snpReport);
				paralogSearchReport(transcriptName, snpPos, refAlleleAtLast, bestBLASTResultsFilteredFile, snpReport);
				// **********************************************************************************************************************************************************************************************************************************************************************************************************************************
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	private static String getCIGARCorrectedRead(SAMRecord rec, IndexedFastaSequenceFile faidx)
	{
		//System.out.println("getCIGARCorrectedRead");
		
		byte ref_sequence[] = faidx.getSubsequenceAt(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd()).getBases();
		byte read_sequence[] = rec.getReadBases();
		int ref_pos = 0;
		int read_pos = 0;
		Cigar cigar = rec.getCigar();
		StringBuilder refseq = new StringBuilder();
		StringBuilder readseq = new StringBuilder();
		//System.out.println("CigarElements:");
		for (CigarElement c : cigar.getCigarElements())
		{
			//System.out.println(c.getOperator()+""+c.getLength());
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
						//System.out.println("appending to read " + ch );
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
					throw new IllegalArgumentException(c.getOperator() + "\n" +

					refseq + "\n" + readseq + "\n" + new String(ref_sequence) + "\n" + new String(read_sequence) + "\n" + cigar);
				}
			}
		}
		//convert the string buffer to a string
		String newReadSeq = readseq.toString();
		//check whether the first CIGAR element is a soft clip
		CigarElement el = cigar.getCigarElements().get(0);
		//if yes, trim a sequence of the corresponding length from the start of the read so that the coordinate system is in keeping with what we see in Tablet and for our SNP discovery
		if(el.getOperator().toString().equals("S"))
		{
			//System.out.println("soft clip found at read start");
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
		//String transcriptAndpositionLiteral = transcriptName+"_"+String.valueOf(snpPos);
		//alleleQuantitiesLookup.put(transcriptName, contentsRelated2AlleleQuantities);
		//infoAboutRefAllele.put(transcriptAndpositionLiteral, refAlleleAtLast);
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
		public static File createFastaFileFromOverlappingReads(String transcriptName, int snpPos, HashMap <String,String> overlappingReadsLookup) throws IOException
		{
			// naming the file with the sequence name + _OverlappingReads.fasta name/extension
			File overlappingReadsFile = new File(transcriptName+"_"+snpPos+"_OverlappingReads.fasta");
			// opening an output Stream to that file
			BufferedWriter bw = new BufferedWriter(new FileWriter(overlappingReadsFile));
			for (String readSequence : overlappingReadsLookup.keySet())
			{
				// writing sequence name FASTA style
				bw.write(">"+overlappingReadsLookup.get(readSequence));
				// newline method putting the plataform specific new line char
				bw.newLine();
				bw.write(readSequence);
				// newline method putting the plataform specific new line char
				bw.newLine();
				// closing the output Stream
			}
			bw.close();
			System.out.println("\ndone -- Fasta file processed with: "+overlappingReadsLookup.size()+" reads.");
			return overlappingReadsFile;
		}
		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
		private static File createBestBLASTResultsFilteredFile(File blastResults4TrascriptAndSNPposition, String transcriptName, int snpPos)
		{
			System.out.println("\nReading the "+blastResults4TrascriptAndSNPposition.getName()+" file and retrieving only the best BLAST results...");
			Integer lineCount = 0;
			Integer capturedLineCount = 0;
			File bestBLASTResultsFilteredFile = new File(transcriptName+"_"+snpPos+"_bestBLAST.txt");
			try
			{
				// opening an output Stream to the file to be written
				BufferedWriter bw = new BufferedWriter(new FileWriter(bestBLASTResultsFilteredFile));
				// reading the input file
				FileInputStream fstream = new FileInputStream(blastResults4TrascriptAndSNPposition);
				DataInputStream in = new DataInputStream(fstream);
				BufferedReader br = new BufferedReader(new InputStreamReader(in));
				String line = null;
				while ((line = br.readLine()) != null)
				{
					//Float similarity = new Float("0.0f");
					//Float bitScore = new Float("0.0f");
					// System.out.println(line);
					lineCount++;
					String[] lineArray = line.split("\t");
					String readName = lineArray[0];
					String ak = lineArray[1];
					//similarity = Float.parseFloat(lineArray[2]);
					String similarity_Brute = lineArray[2].trim();
					System.out.println("\nsimilarity_Brute = "+similarity_Brute);
					Integer similarity = Integer.parseInt(stripExtension(similarity_Brute));
					System.out.println("\nsimilarity = "+similarity);
					Integer alignmentLength = Integer.parseInt(lineArray[3]);
					Integer startOfAlignment = Integer.parseInt(lineArray[6]);
					Integer endOfAlignment = Integer.parseInt(lineArray[7]);
					//bitScore = Float.parseFloat(lineArray[11]);
					String bitScore_Brute = stripExtension(lineArray[11].trim());
					System.out.println("\nbitScore_Brute = "+bitScore_Brute);
					Integer bitScore = Integer.parseInt(bitScore_Brute);
					System.out.println("\nbitScore = "+bitScore);
					System.out.println("\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
					System.out.println("\nPosition at lineArray[0] which would represent readName: "+readName);
					System.out.println("\nPosition at lineArray[1] which would represent ak: "+ak);
					System.out.println("\nPosition at lineArray[2] which would represent similarity: "+similarity);
					System.out.println("\nPosition at lineArray[3] which would represent alignmentLength: "+alignmentLength);
					System.out.println("\nPosition at lineArray[6] which would represent startOfAlignment: "+startOfAlignment);
					System.out.println("\nPosition at lineArray[7] which would represent endOfAlignment: "+endOfAlignment);
					System.out.println("\nPosition at lineArray[11] which would represent bitScore: "+bitScore);
					System.out.println("\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
					if (startOfAlignment.equals(1)) System.out.println("\nYes, startOfAlignment = 1");
					if (endOfAlignment.equals(alignmentLength)) System.out.println("\nYes, endOfAlignment = alignmentLength");
					if (similarity.equals(100.0)) System.out.println("\nYes, similarity = 100.0");
					if (similarity.equals(100.00)) System.out.println("\nYes, similarity = 100.00");
					if (similarity.equals(100)) System.out.println("\nYes, similarity = 100");
					System.out.println("\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
					if ((startOfAlignment.equals(1)) && (endOfAlignment.equals(alignmentLength)) && (similarity.equals(100.0) || similarity.equals(100.00) || similarity.equals(100)))
					{
						if (debug) System.out.println("\nI've entered in the branch where I should capture a valid line...");
						// naming the file with the blastResults4TrascriptAndSNPposition.getName() name + _bestBLAST.txt name/extension
						capturedLineCount++;
						// writing the line with the readName and the respective perfectly matched AK
						bw.write(readName+"\t"+ak+"\t"+similarity+"\t"+bitScore);
						// newline method putting the plataform specific new line char
						bw.newLine();
						// closing the output Stream
					}
				}
				in.close();
				bw.close();
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
			System.out.println("\ndone -- File "+transcriptName+"_"+snpPos+"_bestBLAST.txt created with "+capturedLineCount+" filtered lines out of the original "+lineCount+" lines.");
			lineCount = 0;
			capturedLineCount = 0;
			return bestBLASTResultsFilteredFile;
		}
		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		private static void populateBetweenBLASTStagesRecord(File blastResults4TrascriptAndSNPposition, File bestBLASTResultsFilteredFile, SNPReport snpReport)
		{
			boolean hasBLASTHitBeforeFilter = true;
			boolean hasBLASTHitAfterFilter = true;
			long length1 = blastResults4TrascriptAndSNPposition.length();
			long length2 = bestBLASTResultsFilteredFile.length();
			if (length1 == 0) 
			{
				hasBLASTHitBeforeFilter = false;
			}
			if (length2 == 0) 
			{
				hasBLASTHitAfterFilter = false;
			}
			boolean[] betweenBLASTStagesRecord = {hasBLASTHitBeforeFilter, hasBLASTHitAfterFilter};
			snpReport.hasBLASTResultsBeforeFilter = betweenBLASTStagesRecord[0];
			snpReport.hasBLASTResultsAfterFilter = betweenBLASTStagesRecord[1];
			System.out.println("Checking the value of the 'snpReport.hasBLASTResultsBeforeFilter' bin... "+snpReport.hasBLASTResultsBeforeFilter);
			System.out.println("Checking the value of the 'snpReport.hasBLASTResultsAfterFilter' bin... "+snpReport.hasBLASTResultsAfterFilter);
		}
		
		// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		private static void paralogSearchReport(String transcriptName, int snpPos, char refAlleleAtLast, File bestBLASTResultsFilteredFile, SNPReport snpReport)
		{
			//make new hashmap of AKs as keys and hashsets of reads as values
			HashMap<String, HashSet<String>> akVersusReadsLookup = new HashMap<String, HashSet<String>>();
			//iterate over file line by line
			try
			{
				if (debug)
				{
					System.out.println("\nMaking akVersusReadsLookup 'data structure'...");
				}
				// reading the input file to populate the desired map
				FileInputStream fstream = new FileInputStream(bestBLASTResultsFilteredFile);
				DataInputStream in = new DataInputStream(fstream);
				BufferedReader br = new BufferedReader(new InputStreamReader(in));
				String line = null;
				int lineCount = 0;
				if (debug)
				{
					System.out.println("\n-------------------------------------------- Beginning of the file analysis -----------------------------------------------");
				}
				while ((line = br.readLine()) != null)
				{
					System.out.println(line);
					lineCount++;
					String[] lineArray = line.split("\t");
					String readName = lineArray[0];
					String ak = lineArray[1];
					//lookup AK in map
					HashSet<String> uniqueListOfReads = akVersusReadsLookup.get(ak);
					//if this is null, make a new one
					if(uniqueListOfReads == null)
					{
						//make new one if not found in the akVersusReadsLookup map
						uniqueListOfReads = new HashSet<String>();
						//add it to the akVersusReadsLookup map
						akVersusReadsLookup.put(ak, uniqueListOfReads);
					}
					//add the read to the existing hashset in the map for that specific ak
					uniqueListOfReads.add(readName);
				}
				if (debug)
				{
					System.out.println("\n-------------------------------------------------- End of the file analysis --------------------------------------------------");
					System.out.println("\nFile had "+lineCount+" lines.");
					System.out.println("\nNow, revealing the contents of the akVersusReadsLookup built 'data structure' for debugging purposes...");
					// Get a set of the entries 
					Set<Entry<String, HashSet<String>>> set = akVersusReadsLookup.entrySet(); 
					// Get an iterator 
					Iterator<Entry<String, HashSet<String>>> i = set.iterator(); 
					// Display elements 
					while(i.hasNext()) 
					{ 
						System.out.println();
						Entry<String, HashSet<String>> me = i.next(); 
						System.out.print(me.getKey()+": "); 
						System.out.println(me.getValue());
						//System.out.print(akVersusReadsLookup.keySet());
						//System.out.print(akVersusReadsLookup.values());
					}
				}
				in.close();
				//now do something with the populated map
				lineCount = 0;
				System.out.println("\nFound "+akVersusReadsLookup.size()+" UNIQUE BLAST hits in the target database for the "+transcriptName+"_"+snpPos+" entry...");
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
			System.out.println("\n--------------------------------------------------------------------------------------------------------------");
			System.out.println("Now, trying to find the alleles associated with the previously found target database BLAST hits...");
			System.out.println("--------------------------------------------------------------------------------------------------------------");
			if (akVersusReadsLookup.isEmpty() == false) 
			{
				HashSet<String> tempReadsHS = new HashSet<String>();
				// Get a set of the entries in the akVersusReadsLookup map
				Set<Entry<String, HashSet<String>>> set2 = akVersusReadsLookup.entrySet(); 
				// Get an iterator 
				Iterator<Entry<String, HashSet<String>>> i = set2.iterator(); 
				// Process the elements of akVersusReadsLookup map
				ArrayList<String> tempListOfAKReps = new ArrayList<String>();
				while(i.hasNext()) 
				{ 
					// Auxiliary variables of this processing (and "zeroing" them before each AK iteration)
					int accumBaseA2 = 0;
					int accumBaseC2 = 0;
					int accumBaseT2 = 0;
					int accumBaseG2 = 0;
					int accumUndetermined2 = 0;
					int accumRefAlleleEqualsA = 0;
					int accumRefAlleleEqualsC = 0;
					int accumRefAlleleEqualsT = 0;
					int accumRefAlleleEqualsG = 0;
					int accumRefAlleleEqualsUndetermined = 0;
					char flagRefAllele = 'N';
					String tempAssembledString4ParalogReport = null;
					// For refreshing the mind: the objects (and a copy of there declarations) which have some sort of persistence to be searched through: infoAboutReadAndAllele.put(transcriptName+"_"+snpPos, readsAndAllelesLookup);
					// readsAndAllelesLookup.put(readName, allele);
					HashMap<String, Character> tempReadsAndAllelesLookup = infoAboutReadAndAllele.get(transcriptName+"_"+snpPos);
					// iterate over the reads for each AK
					Entry<String, HashSet<String>> me = i.next(); 
					tempReadsHS = me.getValue();
					Iterator<String> it = tempReadsHS.iterator();
					while(it.hasNext())
					{
						String entry = it.next();
						char alleleOfRead = tempReadsAndAllelesLookup.get(entry);
						// Trying to quantify the different alleles found in each read overlapping the SNP position
						if ((alleleOfRead == 'A') || (alleleOfRead == 'a'))
						{
							accumBaseA2++;
							if (alleleOfRead == refAlleleAtLast) 
							{
								accumRefAlleleEqualsA++;
								flagRefAllele = 'Y';
							}
							
						} else if ((alleleOfRead == 'C') || (alleleOfRead == 'c'))
						{
							accumBaseC2++;
							if (alleleOfRead == refAlleleAtLast) 
							{
								accumRefAlleleEqualsC++;
								flagRefAllele = 'Y';
							}
						} else if ((alleleOfRead == 'T') || (alleleOfRead == 't'))
						{
							accumBaseT2++;
							if (alleleOfRead == refAlleleAtLast) 
							{
								accumRefAlleleEqualsT++;
								flagRefAllele = 'Y';
							}
						} else if ((alleleOfRead == 'G') || (alleleOfRead == 'g'))
						{
							accumBaseG2++;
							if (alleleOfRead == refAlleleAtLast) 
							{
								accumRefAlleleEqualsG++;
								flagRefAllele = 'Y';
							}
						} else
						{
							accumUndetermined2++;
							if (alleleOfRead == refAlleleAtLast) 
							{
								accumRefAlleleEqualsUndetermined++;
								flagRefAllele = 'Y';
							}
						}
					// last line of inner read while
					}
					if (debug)
					{
						System.out.println("\n--------------------------------------------------------------------------------------------------------------");
						System.out.println("Printing a very basic report for the "+me.getKey()+" target database BLAST hit which was found...");
						System.out.println("\nDoes "+me.getKey()+" have the reference allele "+refAlleleAtLast+" present at the SNP position? "+flagRefAllele);
						System.out.println("\nQuantity of hits for the 'A' allele: "+accumBaseA2+" and, if this allele is the same of the reference allele one, its quantity as well for debugging purposes: "+accumRefAlleleEqualsA);
						System.out.println("\nQuantity of hits for the 'C' allele: "+accumBaseC2+" and, if this allele is the same of the reference allele one, its quantity as well for debugging purposes: "+accumRefAlleleEqualsC);
						System.out.println("\nQuantity of hits for the 'T' allele: "+accumBaseT2+" and, if this allele is the same of the reference allele one, its quantity as well for debugging purposes: "+accumRefAlleleEqualsT);
					        System.out.println("\nQuantity of hits for the 'G' allele: "+accumBaseG2+" and, if this allele is the same of the reference allele one, its quantity as well for debugging purposes: "+accumRefAlleleEqualsG);
						System.out.println("\nQuantity of hits for any 'undetermined' allele: "+accumUndetermined2+" and, if this allele is the same of the reference allele one, its quantity as well for debugging purposes: "+accumRefAlleleEqualsUndetermined);
						System.out.println("--------------------------------------------------------------------------------------------------------------");
					}
					tempAssembledString4ParalogReport = me.getKey()+":"+flagRefAllele+":"+accumBaseA2+":"+accumBaseC2+":"+accumBaseT2+":"+accumBaseG2+":"+accumUndetermined2+";";
					tempListOfAKReps.add(tempAssembledString4ParalogReport);
				//last line of AK while
				}
				if (debug)
				{
					System.out.println("\nRevealing the contents of tempListOfAKReps 'data structure' for debugging purposes...");
					System.out.println();
					Iterator<String> itArrayList = tempListOfAKReps.iterator();
					String listString = "";
					while (itArrayList.hasNext()) 
					{
						listString += itArrayList.next();
					}
					System.out.println(listString);
					listString = "";
				}
				Iterator<String> itArrayList = tempListOfAKReps.iterator();
				String listString = "";
				while (itArrayList.hasNext()) 
				{
					listString += itArrayList.next();
				}
				// snpReport member taken from this method
				snpReport.paralogReport = listString;
				listString = "";
				tempListOfAKReps.clear();
			// last line if akVersusReadsLookup.isEmpty() == false
			} 
			else 
			{
				if (debug) 
				{
					System.out.println("\n--------------------------------------------------------------------------------------------------------------");
					System.out.println("NO HIT FOUND in the target BLAST database for the "+transcriptName+"_"+snpPos+" entry...");
					System.out.println("--------------------------------------------------------------------------------------------------------------");
				}
			}
		// last line of the method paralogSearchReport
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
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// "Garbage" during development time:
	/*tempAssembledString4ParalogReport = me.getKey()+":"+flagRefAllele+":"+accumBaseA2+":"+accumBaseC2+":"+accumBaseT2+":"+accumBaseG2+":"+accumUndetermined2+";";
	tempListOfAKReps.add(tempAssembledString4ParalogReport);
	//last line of AK while
	}
	if (debug)
	{
		System.out.println("\nRevealing the contents of tempListOfAKReps 'data structure' for debugging purposes...");
		System.out.println();
		Iterator<String> itArrayList = tempListOfAKReps.iterator();
		String listString = "";
		while (itArrayList.hasNext()) 
		{
			listString += itArrayList.next();
		}
		System.out.println(listString);
		listString = "";
	}
	Iterator<String> itArrayList = tempListOfAKReps.iterator();
	String listString = "";
	while (itArrayList.hasNext()) 
	{
		listString += itArrayList.next();
	}
		
		// snpReport member taken from this method
		//snpReport.paralogReport = listString;
		//listString = "";
		//tempListOfAKReps.clear();
		 */
}
