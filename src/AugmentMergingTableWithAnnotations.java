import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class AugmentMergingTableWithAnnotations
{
	static String tableFn = "";
	static String vcfFilelist = "";
	static String ofn = "";
		
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
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
			}
			
		}
		
		if(tableFn.length() == 0 || ofn.length() == 0 || vcfFilelist.length() == 0)
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
		System.out.println("Usage: java -cp src AugmentMergingTable [args]");
		System.out.println("Required args:");
		System.out.println("  table_file   (String) - the TSV built from the BuildMergingTable script");
		System.out.println("  out_file     (String) - the name of the file to output the merging table to");
		System.out.println("  vcf_filelist (String) - a txt file with each line containing the filename of a merged vcf");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println();
	}
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		ArrayList<String> filenames = PipelineManager.getFilesFromList(vcfFilelist);
		
		Scanner input = new Scanner(new FileInputStream(new File(tableFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
				
		// A list of all merged variants
		ArrayList<DetailedMergedVariant> variantList = new ArrayList<DetailedMergedVariant>();
		
		// For each input file, a map from variant ID to the index in the merged variant list to which it belongs
		HashMap<String, Integer>[] componentMap = new HashMap[filenames.size()];
		
		for(int i = 0; i<filenames.size(); i++)
		{
			componentMap[i] = new HashMap<String, Integer>();
		}
		
		input.nextLine(); // TODO get the filenames from here instead of filelist
				
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
				}
				else if(id.startsWith("0_") && componentMap[i].containsKey(id.substring(2)))
				{
					int mergedIndex = componentMap[i].get(id.substring(2));
					variantList.get(mergedIndex).integrateVariant(entry, i);
				}
				else
				{
					//System.out.println(i+" "+id);
				}
			}
		}
		
		// Output header
		out.println(DetailedMergedVariant.makeHeader(filenames));
		
		// Output info for each merged variant
		for(DetailedMergedVariant v : variantList)
		{
			v.finalizeVariant();
			out.println(v);
		}
		
		
		out.close();
	}
	
	/*
	 * A representation of which input variants make up a merged variant, with additional details and stats 
	 */
	static class DetailedMergedVariant
	{
		String index;
		ArrayList<String> ids;
		HashSet<Integer> specificSampleSet;
		HashSet<String> annotations;
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
		
		boolean inCentromere = false;
		boolean inRepeat = false;
		boolean inGene = false;
		boolean inExon = false;
		
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
			for(int i = 1; i<tokens.length - 1; i++)
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
			
			annotations = new HashSet<String>();
			String[] annotationArray = tokens[tokens.length - 1].split(",");
			for(String s : annotationArray)
			{
				annotations.add(s);
			}
			
			inGene = annotations.contains("GENE_gene");
			inExon = annotations.contains("GENE_exon");
			inRepeat = false;
			for(String s : annotations)
			{
				if(s.startsWith("REPEAT_"))
				{
					inRepeat = true;
				}
			}
			inCentromere = annotations.contains("CENTROMERE");
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
			res.append("\t" + "GENE_FLAG");
			res.append("\t" + "EXON_FLAG");
			res.append("\t" + "REPEAT_FLAG");
			res.append("\t" + "CENTROMERE_FLAG");
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
			res.append("\t" + (inGene ? 1 : 0));
			res.append("\t" + (inExon ? 1 : 0));
			res.append("\t" + (inRepeat ? 1 : 0));
			res.append("\t" + (inCentromere ? 1 : 0));
			return res.toString();
		}
	}
}
