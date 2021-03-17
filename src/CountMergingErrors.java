import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class CountMergingErrors
{
	static String tableFn = "";
	static String vcfFilelist = "";
	static String discSuppVec = "";
	static String ofn = "";
	static String software = ".";
	static boolean append = false;
		
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
		ArrayList<AugmentMergingTable.DetailedMergedVariant> variantList = new ArrayList<AugmentMergingTable.DetailedMergedVariant>();
		
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
			AugmentMergingTable.DetailedMergedVariant v = new AugmentMergingTable.DetailedMergedVariant(line, filenames);
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
		
		// Output info for each merged variant
		for(AugmentMergingTable.DetailedMergedVariant v : variantList)
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
			
		}
		
		System.out.println();
		System.out.println("Errors for " + software);
		
		System.out.println("Mixed strand and type: " + badBoth);
		System.out.println("Mixed strand only: " + badStrand);
		System.out.println("Mixed type only: " + badType);
		
		if(discSuppVec.length() > 0)
		{
			System.out.println("Mixed strand and type affecting discordance: " + badBothDiscordant);
			System.out.println("Mixed strand only affecting discordance: " + badStrandDiscordant);
			System.out.println("Mixed type only affecting discordance: " + badTypeDiscordant);
		}
		
		System.out.println();
		
		if(ofn.length() > 0)
		{
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(ofn, append)));
			
			// Print header if making a new table
			if(!append)
			{
				out.print("SOFTWARE\tMIXED_STRAND_AND_TYPE\tMIXED_STRAND_ONLY\tMIXED_TYPE_ONLY");
				if(discSuppVec.length() > 0)
				{
					out.print("\tDISC_MIXED_STRAND_AND_TYPE\tDISC_MIXED_STRAND_ONLY\tDISC_MIXED_TYPE_ONLY");
				}
				out.println();
			}
			
			// Print statistics for this software
			out.print(software + "\t" + badBoth + "\t" + badStrand + "\t" + badType);
			if(discSuppVec.length() > 0)
			{
				out.print("\t" + badBothDiscordant + "\t" + badStrandDiscordant + "\t" + badTypeDiscordant);
			}
			out.println();
			out.close();
		}
		
	}
}
