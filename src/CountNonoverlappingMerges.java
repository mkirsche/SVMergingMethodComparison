import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class CountNonoverlappingMerges
{
	static String tableFn = "";
	static String vcfFilelist = "";
	static String discSuppVec = "";
	static String ofn = "";
	static String software = ".";
	static boolean append = false;
	
	static int specificReads = -1;
	static int specificLength = -1;
		
	/*
	 * Parse command line arguments
	 */
	static void parseArgs(String[] args)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf('=');
			if(equalsIdx == -1)
			{
				if(arg.toLowerCase().endsWith("append"))
				{
					append = true;
				}
			}
			else
			{
				String key = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("table_file"))
				{
					tableFn = val;
				}
				else if(key.equalsIgnoreCase("vcf_filelist"))
				{
					vcfFilelist = val;
				}
				else if(key.equalsIgnoreCase("disc_supp_vec"))
				{
					discSuppVec = val;
				}
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
				else if(key.equalsIgnoreCase("software"))
				{
					software = val;
				}
			}
			
		}
		
		if(tableFn.length() == 0 || vcfFilelist.length() == 0)
		{
			usage();
			System.exit(0);
		}
	}
	
	/*
	 * Print a usage menu
	 */
	static void usage()
	{
		System.out.println();
		System.out.println("Usage: java -cp src CountMergingErrors [args]");
		System.out.println("Required args:");
		System.out.println("  table_file    (String) - the TSV built from the BuildMergingTable script");
		System.out.println("  vcf_filelist  (String) - a txt file with each line containing the filename of a merged vcf");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  disc_supp_vec (String) - the support vector corresponding to discordant reads");
		System.out.println("  out_file      (String) - where to write results in table format");
		System.out.println("  software  [.] (String) - which merging software to list in the table (if table is produced)");
		System.out.println("  --append               - append to an existing table instead of making a new one");
		System.out.println();
	}
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		ArrayList<String> filenames = PipelineManager.getFilesFromList(vcfFilelist);
		
		Scanner input = new Scanner(new FileInputStream(new File(tableFn)));
				
		// A list of all merged variants
		ArrayList<DetailedMergedVariant> variantList = new ArrayList<DetailedMergedVariant>();
		
		// For each input file, a map from variant ID to the index in the merged variant list to which it belongs
		HashMap<String, Integer>[] componentMap = new HashMap[filenames.size()];
		
		for(int i = 0; i<filenames.size(); i++)
		{
			componentMap[i] = new HashMap<String, Integer>();
		}
		
		input.nextLine();
				
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			DetailedMergedVariant v = new DetailedMergedVariant(line, filenames);
			for(int i = 0; i<v.ids.size(); i++)
			{
				if(v.ids.get(i).length() > 0)
				{
					String[] ids = v.ids.get(i).split(";");
					for(String id : ids)
					{
						componentMap[i].put(id, variantList.size());
					}
				}
			}
			
			variantList.add(v);
		}
		
		input.close();
		
		for(int i = 0; i<filenames.size(); i++)
		{
			input = new Scanner(new FileInputStream(new File(filenames.get(i))));
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.length() == 0 || line.startsWith("#"))
				{
					continue;
				}
				VcfEntry entry = new VcfEntry(line);
				String id = entry.getId();
				if(componentMap[i].containsKey(id))
				{
					int mergedIndex = componentMap[i].get(id);
					variantList.get(mergedIndex).integrateVariant(entry, i);
					if(discSuppVec.length() > 0 && i == discSuppVec.indexOf('1'))
					{
						variantList.get(mergedIndex).updateDiscordanceReduction(entry);
					}
				}
				else if(id.startsWith("0_") && componentMap[i].containsKey(id.substring(2)))
				{
					int mergedIndex = componentMap[i].get(id.substring(2));
					variantList.get(mergedIndex).integrateVariant(entry, i);
					if(discSuppVec.length() > 0 && i == discSuppVec.indexOf('1'))
					{
						variantList.get(mergedIndex).updateDiscordanceReduction(entry);
					}
				}
				else
				{
					//System.out.println(i+" "+id);
				}
			}
		}
		
		int badStrand = 0, badType = 0, badBoth = 0;
		int badStrandDiscordant = 0, badTypeDiscordant = 0, badBothDiscordant = 0;
		int nonoverlap = 0, nonoverlapDiscordant = 0;
		int ism = 0, ismDiscordant = 0, ismCombinedDisc = 0;
		
		// Output info for each merged variant
		for(DetailedMergedVariant v : variantList)
		{
			v.finalizeVariant();
			int totalStrand = v.numPlusMinus + v.numPlusPlus + v.numMinusMinus + v.numMinusPlus;
			int maxStrand = Math.max(Math.max(v.numPlusMinus, v.numPlusPlus), Math.max(v.numMinusMinus, v.numMinusPlus));
			
			int totalType = v.numIns + v.numDel + v.numInv + v.numDup + v.numTra;
			int maxType = Math.max(Math.max(Math.max(v.numIns, v.numDel), Math.max(v.numInv, v.numDup)), v.numTra);
			
			if(totalType != maxType && totalStrand != maxStrand)
			{
				badBoth++;
				if(v.reducesDiscordance)
				{
					badBothDiscordant++;
				}
			}
			else if(totalStrand != maxStrand)
			{
				badStrand++;
				if(v.reducesDiscordance)
				{
					badStrandDiscordant++;
				}
			}
			else if(totalType != maxType)
			{
				badType++;
				if(v.reducesDiscordance)
				{
					badTypeDiscordant++;
				}
			}
			else if(v.numVars > v.supp)
			{
				ism++;
				if(v.reducesDiscordance)
				{
					ismDiscordant++;
				}
				if(v.suppVec.equals(discSuppVec))
				{
					ismCombinedDisc+= v.numVars - 1;
				}
			}
			else if(v.numVars > 1 && !v.overlap)
			{
				nonoverlap++;
				if(v.reducesDiscordance)
				{
					nonoverlapDiscordant++;
				}
			}
			
		}
		
		System.out.println();
		System.out.println("Errors for " + software);
		
		System.out.println("Mixed strand and type: " + badBoth);
		System.out.println("Mixed strand only: " + badStrand);
		System.out.println("Mixed type only: " + badType);
		System.out.println("Non-overlapping merges: " + nonoverlap);
		System.out.println("Intrasample merges: " + ism);

		
		if(discSuppVec.length() > 0)
		{
			System.out.println("Mixed strand and type affecting discordance: " + badBothDiscordant);
			System.out.println("Mixed strand only affecting discordance: " + badStrandDiscordant);
			System.out.println("Mixed type only affecting discordance: " + badTypeDiscordant);
			System.out.println("Non-overlapping merges affecting discordance: " + nonoverlapDiscordant);
			System.out.println("Intrasample merges affecting discordance: " + ismDiscordant);
		}
		
		System.out.println();
		
		if(ofn.length() > 0)
		{
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(ofn, append)));
			
			// Print header if making a new table
			if(!append)
			{
				out.print("SOFTWARE\tMIXED_STRAND_AND_TYPE\tMIXED_STRAND_ONLY\tMIXED_TYPE_ONLY\tNONOVERLAP\tISM");
				if(discSuppVec.length() > 0)
				{
					out.print("\tDISC_MIXED_STRAND_AND_TYPE\tDISC_MIXED_STRAND_ONLY\tDISC_MIXED_TYPE_ONLY\tDISC_NONOVERLAP\tDISC_ISM\tDISC_COMBINED_ISM");
				}
				out.println();
			}
			
			// Print statistics for this software
			out.print(software + "\t" + badBoth + "\t" + badStrand + "\t" + badType + "\t" + nonoverlap + "\t" + ism);
			if(discSuppVec.length() > 0)
			{
				out.print("\t" + badBothDiscordant + "\t" + badStrandDiscordant + "\t" + badTypeDiscordant + "\t" + nonoverlapDiscordant + 
						"\t" + ismDiscordant + "\t" + ismCombinedDisc);
			}
			out.println();
			out.close();
		}
		
	}
	
	static class DetailedMergedVariant
	{
		String index;
		ArrayList<String> ids;
		HashSet<Integer> specificSampleSet;
		String suppVec;
		int supp;
		
		int numIns = 0, numDel = 0, numInv = 0, numTra = 0, numDup = 0;
		String chr = "";
		String type = "";
		
		int maxStart = -1, minStart = -1, maxEnd = -1, minEnd = -1;
		int start = -1, end = -1, len = -1;
		
		int numPlusPlus = 0, numMinusMinus = 0, numPlusMinus = 0, numMinusPlus = 0;
		
		int numVars = 0;
		
		boolean isSpecific = false;
		boolean isPrecise = false;

		int specificCount = 0;
		
		boolean reducesDiscordance = false;
		
		boolean overlap = false;
				
		void updateDiscordanceReduction(VcfEntry entry) throws Exception
		{
			boolean entrySpecific = false, entryPrecise = false;
			if(entry.hasInfoField("IS_SPECIFIC") && entry.getInfo("IS_SPECIFIC").equals("1"))
			{
				entrySpecific = true;
			}
			
			if(entry.tabTokens[7].contains(";PRECISE;") || entry.tabTokens[7].startsWith("PRECISE;"))
			{
				entryPrecise = true;
			}
			if(entrySpecific && entryPrecise)
			{
				reducesDiscordance = true;
			}
		}
		
		void integrateVariant(VcfEntry entry, int sampleNum) throws Exception
		{
			String chr = entry.getChromosome();
			int pos = (int)entry.getPos();
			int end = (int)entry.getEnd();
			String type = entry.getNormalizedType();
			int len = entry.getLength();
			
			if(type.equals("TRA") && this.chr.length() > 0 && !chr.equals(this.chr))
			{
				chr = entry.getChr2();
				int temp = pos;
				pos = end;
				end = temp;
			}
			
			if(pos <= maxEnd && end >= maxStart)
			{
				overlap = true;
			}
			
			if(pos > maxStart)
			{
				maxStart = pos;
			}
			if(minStart == -1 || pos < minStart)
			{
				minStart = pos;
			}
			
			if(end > maxEnd)
			{
				maxEnd = end;
			}
			if(minEnd == -1 || end < minEnd)
			{
				minEnd = end;
			}
			
			if(start == -1)
			{
				this.chr = chr;
				this.start = pos;
				this.end = end;
				this.len = len;
			}
			
			if(type.equals("INS"))
			{
				numIns++;
			}
			else if(type.equals("DEL"))
			{
				numDel++;
			}
			else if(type.equals("INV"))
			{
				numInv++;
			}
			else if(type.equals("TRA"))
			{
				numTra++;
			}
			else if(type.equals("DUP"))
			{
				numDup++;
			}
			
			String strand = entry.getStrand();
			if(strand.equals("++"))
			{
				numPlusPlus++;
			}
			else if(strand.equals("--"))
			{
				numMinusMinus++;
			}
			else if(strand.equals("+-"))
			{
				numPlusMinus++;
			}
			else if(strand.equals("-+"))
			{
				numMinusPlus++;
			}
			
			if(specificLength != -1)
			{
				boolean inSpecific = false;
				int readSupport = entry.getReadSupport();
				boolean longEnough = entry.getType().equals("TRA") || entry.getType().equals("BND") || 
						Math.abs(entry.getLength()) >= specificLength || entry.getLength() == 0;
				
				if(readSupport >= specificReads && longEnough)
				{
					inSpecific = true;
				}
				entry.setInfo("IS_SPECIFIC", inSpecific ? "1": "0");
			}
			if(entry.hasInfoField("IS_SPECIFIC") && entry.getInfo("IS_SPECIFIC").equals("1"))
			{
				
				specificCount++;
				specificSampleSet.add(sampleNum);
				isSpecific = true;
			}
			
			if(entry.tabTokens[7].contains(";PRECISE;") || entry.tabTokens[7].startsWith("PRECISE;"))
			{
				isPrecise = true;
			}
			
			numVars++;
			
		}
		
		void finalizeVariant()
		{
			// Get type basd on what shows up the most in merged variants
			type = "INS";
			int count = numIns;
			if(numDel > count)
			{
				type = "DEL";
				count = numDel;
			}
			if(numInv > count)
			{
				type = "INV";
				count = numInv;
			}
			if(numTra > count)
			{
				type = "TRA";
				count = numTra;
			}
			if(numDup > count)
			{
				type = "DUP";
				count = numDup;
			}
		}
		DetailedMergedVariant(String line, ArrayList<String> filenames)
		{
			String[] tokens = line.split("\t");
			index = tokens[0];
			ids = new ArrayList<String>();
			specificSampleSet = new HashSet<Integer>();
			StringBuilder suppVecBuilder = new StringBuilder("");
			for(int i = 1; i<tokens.length; i++)
			{
				boolean isPresent = !tokens[i].equals(".");
				ids.add(isPresent ? tokens[i] : "");
				if(isPresent)
				{
					supp++;
					suppVecBuilder.append("1"); 
				}
				else
				{
					suppVecBuilder.append("0");
				}
			}
			suppVec = suppVecBuilder.toString();
		}
		
		static String makeHeader(ArrayList<String> filenames)
		{
			StringBuilder res = new StringBuilder("");
			res.append("INDEX");
			for(String fn : filenames)
			{
				res.append("\t");
				res.append(fn);
			}
			res.append("\t" + "CHR");
			res.append("\t" + "START");
			res.append("\t" + "END");
			res.append("\t" + "LEN");
			res.append("\t" + "TYPE");
			res.append("\t" + "SUPP_VEC");
			res.append("\t" + "SUPP");
			res.append("\t" + "MIN_START");
			res.append("\t" + "MAX_START");
			res.append("\t" + "MIN_END");
			res.append("\t" + "MAX_END");
			res.append("\t" + "SPECIFIC_FLAG");
			res.append("\t" + "PRECISE_FLAG");
			res.append("\t" + "NUMVARS");
			res.append("\t" + "SPECIFIC_COUNT");
			res.append("\t" + "SPECIFIC_SAMPLES");
			return res.toString();
		}
		
		public String toString()
		{
			StringBuilder res = new StringBuilder("");
			res.append(index);
			for(String id : ids)
			{
				res.append("\t");
				res.append(id.length() > 0 ? id : ".");
			}
			res.append("\t" + chr);
			res.append("\t" + start);
			res.append("\t" + end);
			res.append("\t" + len);
			res.append("\t" + type);
			res.append("\t" + suppVec);
			res.append("\t" + supp);
			res.append("\t" + minStart);
			res.append("\t" + maxStart);
			res.append("\t" + minEnd);
			res.append("\t" + maxEnd);
			res.append("\t" + (isSpecific ? 1 : 0));
			res.append("\t" + (isPrecise ? 1 : 0));
			res.append("\t" + numVars);
			res.append("\t" + specificCount);
			res.append("\t" + specificSampleSet.size());
			return res.toString();
		}
	}
}
