/*
 * Builds a TSV giving the IDs making up each merged variant in a VCF
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;

public class BuildMergingTable
{
	static String vcfFn = "";
	static String vcfFilelist = "";
	static String ofn = "";
	static String mode = "jasmine";
	static String sampleList = "";
	
	static HashMap<String, Integer> sampleToIndex;
	static HashMap<String, String>[] coordsToId;
	
	static void buildSampleMap() throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(vcfFilelist)));
		int index = 0;
		sampleToIndex = new HashMap<String, Integer>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			String command ="cat " + line + " | grep '#' | tail -n 1 | cut -f 10";
			String[] fullCommand = new String[] {"/bin/sh", "-c", command};
			Process child = Runtime.getRuntime().exec(fullCommand);
			InputStream seqStream = child.getInputStream();
			Scanner seqInput = new Scanner(seqStream);
			String sample = seqInput.next();
			sampleToIndex.put(sample, index);
			index++;
			seqInput.close();
		}
	}
	
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
				if(key.equalsIgnoreCase("vcf_file"))
				{
					vcfFn = val;
				}
				else if(key.equalsIgnoreCase("vcf_filelist"))
				{
					vcfFilelist = val;
				}
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
				else if(key.equalsIgnoreCase("sample_list"))
				{
					sampleList = val;
				}
				else if(key.equalsIgnoreCase("mode"))
				{
					mode = val;
				}
			}
			
		}
		
		if(vcfFn.length() == 0 || ofn.length() == 0 || vcfFilelist.length() == 0)
		{
			usage();
			System.exit(0);
		}
		
		if(!mode.equals("jasmine") && !mode.equals("survivor") && !mode.equals("svtools") && !mode.equals("svimmer") && !mode.equals("jasmine_intra")
				&& !mode.equals("svpop") && !mode.equals("dbsvmerge"))
		{
			usage();
			System.exit(0);
		}	
		if(mode.equals("svpop") && sampleList.length() == 0)
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
		System.out.println("Usage: java -cp src BuildMergingTable [args]");
		System.out.println("Required args:");
		System.out.println("  vcf_file     (String) - the path to the merged VCF");
		System.out.println("  out_file     (String) - the name of the file to output the merging table to");
		System.out.println("  vcf_filelist (String) - a txt file with each line containing the filename of a merged vcf");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  mode       (String) [jasmine] - the merging software used: jasmine, survivor, svtools, jasmine_intra, svpop, dbsvmerge, or svimmer");
		System.out.println("  sample_list (String) []        - comma-separated list of sample names in same order as in table - required when using svpop");
		System.out.println();
	}
	
	/*
	 * Builds a map from coordinates string (in format chr1_pos-chr2_end) to ID for each sample
	 */
	@SuppressWarnings("unchecked")
	static void buildCoordsMap() throws Exception
	{
		Scanner filelistInput = new Scanner(new FileInputStream(new File(vcfFilelist)));
		ArrayList<String> filelist = PipelineManager.getFilesFromList(vcfFilelist);
		coordsToId = new HashMap[filelist.size()];
		for(int i = 0; i<filelist.size(); i++)
		{
			coordsToId[i] = new HashMap<String, String>();
			String vcf = filelist.get(i);
			Scanner input = new Scanner(new FileInputStream(new File(vcf)));
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.length() == 0 || line.startsWith("#"))
				{
					continue;
				}
				VcfEntry entry = new VcfEntry(line);
				String coords = entry.getChromosome() + "_" + entry.getPos() + "-" +
						entry.getInfo("CHR2") + "_" + entry.getEnd();
				coordsToId[i].put(coords, entry.getId());
			}
			input.close();
		}
		filelistInput.close();
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		if(mode.equals("svtools") || mode.equals("svimmer") || mode.equals("dbsvmerge"))
		{		
			buildSampleMap();
		}
		
		if(mode.equals("survivor"))
		{
			buildCoordsMap();
		}
		
		ArrayList<String> filenames = PipelineManager.getFilesFromList(vcfFilelist);
		
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		ArrayList<SimpleMergedVariant> variants = new ArrayList<SimpleMergedVariant>();
		
		String headerLine = "";
		
		String lastId = "";
		ArrayList<String> dbSamples = new ArrayList<String>(), dbIds = new ArrayList<String>();
		
		while(input.hasNext())
		{
			String line = input.nextLine();
			
			if(mode.equals("svpop"))
			{
				if(headerLine.length() == 0)
				{
					headerLine = line;
				}
				else
				{
					variants.add(new SimpleMergedVariant(headerLine, line));
				}
			}
			else if(mode.equals("dbsvmerge"))
			{
				String[] tokens = line.split("\t");
				String id = tokens[0];
				String sample = tokens[7];
				String varId = tokens[8];
				if(lastId.equals(id))
				{
					dbSamples.add(sample);
					dbIds.add(varId);
				}
				if(!lastId.equals(id) || !input.hasNext())
				{
					if(dbSamples.size() > 0)
					{
						variants.add(new SimpleMergedVariant(dbSamples, dbIds));
					}
					
					if(input.hasNext())
					{
						dbSamples = new ArrayList<String>();
						dbSamples.add(sample);
						dbIds = new ArrayList<String>();
						dbIds.add(varId);
					}
				}
				lastId = id;
			}
			else
			{
				if(line.length() == 0 || line.startsWith("#"))
				{
					continue;
				}
				variants.add(new SimpleMergedVariant(new VcfEntry(line)));
			}
		}
		
		Collections.sort(variants);
		
		out.println(SimpleMergedVariant.makeHeader(filenames));
		for(SimpleMergedVariant v : variants)
		{
			out.println(v);
		}
		
		input.close();
		out.close();
	}
	
	/*
	 * Get the list of IDs merged into a given variant.  There are different cases depending on the caller.
	 * There is an implemenation for the following modes (choices of merging software):
	 *   jasmine
	 */
	static ArrayList<String> getIds(VcfEntry entry) throws Exception
	{
		/*
		 * Get list of Jasmine IDs by parsing the support vector
		 */
		if(mode.equals("jasmine"))
		{
			String suppVec = entry.getInfo("SUPP_VEC");
			String idList = entry.getInfo("IDLIST");
			String[] ids = idList.split(",");
			int n = suppVec.length();
			
			ArrayList<String> res = new ArrayList<String>();
			
			int idIndex = 0;
			
			for(int i = 0; i<n; i++)
			{
				if(suppVec.charAt(i) == '0')
				{
					res.add("");
				}
				else
				{
					res.add(ids[idIndex]);
					idIndex++;
				}
			}
			
			return res;
		}
		
		if(mode.equals("jasmine_intra"))
		{
			String suppVec = entry.getInfo("SUPP_VEC");
			String idList = entry.getInfo("INTRASAMPLE_IDLIST");
			String[] ids = idList.split("\\.");
			int n = suppVec.length();
			
			ArrayList<String> res = new ArrayList<String>();
			
			int idIndex = 0;
			
			for(int i = 0; i<n; i++)
			{
				if(suppVec.charAt(i) == '0')
				{
					res.add("");
				}
				else
				{
					res.add(ids[idIndex].replaceAll(",", ";"));
					idIndex++;
				}
			}
			
			return res;
		}
		
		else if(mode.equals("survivor"))
		{
			String suppVec = entry.getInfo("SUPP_VEC");
			int n = suppVec.length();
			
			ArrayList<String> res = new ArrayList<String>();
						
			for(int i = 0; i<n; i++)
			{
				if(suppVec.charAt(i) == '0')
				{
					res.add("");
				}
				else
				{
					String[] coordsList = entry.tabTokens[9+i].split(":")[10].split(",");
					StringBuilder idList = new StringBuilder("");
					for(int j = 0; j<coordsList.length; j++)
					{
						if(j > 0) idList.append(';');
						idList.append(coordsToId[i].get(coordsList[j]));
					}
					res.add(idList.toString());
				}
			}
			
			return res;
		}
		
		else if(mode.equals("svtools"))
		{
			String sname = entry.getInfo("SNAME");
			String[] tokens = sname.split(",");
			StringBuilder[] answers = new StringBuilder[sampleToIndex.size()];
			for(int i = 0; i<answers.length; i++)
			{
				answers[i] = new StringBuilder("");
			}
			for(String token : tokens)
			{
				int colonIdx = token.indexOf(':');
				String sample = token.substring(0, colonIdx);
				String id = token.substring(1 + colonIdx);
				int sampleIndex = sampleToIndex.get(sample);
				if(answers[sampleIndex].length() > 0)
				{
					answers[sampleIndex].append(';');
				}
				answers[sampleIndex].append(id);
			}
			
			ArrayList<String> res = new ArrayList<String>();
			for(StringBuilder answer : answers) res.add(answer.toString());
			return res;
		}
		
		else if(mode.equals("svimmer"))
		{
			String sname = entry.getInfo("MERGED_IDS");
			String[] tokens = sname.split(",");
			StringBuilder[] answers = new StringBuilder[sampleToIndex.size()];
			for(int i = 0; i<answers.length; i++)
			{
				answers[i] = new StringBuilder("");
			}
			for(String token : tokens)
			{
				int colonIdx = token.indexOf(':');
				String sample = token.substring(0, colonIdx);
				String id = token.substring(1 + colonIdx);
				int sampleIndex = sampleToIndex.get(sample);
				if(answers[sampleIndex].length() > 0)
				{
					answers[sampleIndex].append(';');
				}
				answers[sampleIndex].append(id);
			}
			
			ArrayList<String> res = new ArrayList<String>();
			for(StringBuilder answer : answers) res.add(answer.toString());
			return res;
		}
		
		return new ArrayList<String>();
	}
	
	static ArrayList<String> getSvPopIds(String headerLine, String line)
	{
		HashMap<String, Integer> sampleToIndexMap = new HashMap<String, Integer>();
		String[] samples = sampleList.split(",");
		int numSamples = samples.length;
		for(int i = 0; i<numSamples; i++)
		{
			sampleToIndexMap.put(samples[i], i);
		}
		
		ArrayList<String> res = new ArrayList<String>();
		String[] headerTokens = headerLine.split("\t");
		int sampleIndex = -1, idIndex = -1;
		for(int i = 0; i<headerTokens.length; i++)
		{
			if(headerTokens[i].equals("MERGE_SAMPLES"))
			{
				sampleIndex = i;
			}
			else if(headerTokens[i].equals("MERGE_IDLIST"))
			{
				idIndex = i;
			}
		}
		if(sampleIndex == -1 || idIndex == -1)
		{
			return res;
		}
		String[] tokens = line.split("\t");
		samples = tokens[sampleIndex].split(",");
		String[] ids = tokens[idIndex].split(",");
		String[] outputTokens = new String[numSamples];
		for(int i = 0; i<outputTokens.length; i++)
		{
			outputTokens[i] = "";
		}
		//Arrays.fill(outputTokens, ".");
		for(int i = 0; i<ids.length; i++)
		{
			int sample = sampleToIndexMap.get(samples[i]);
			outputTokens[sample] = ids[i].substring(1 + ids[i].lastIndexOf('-'));
		}
		for(int i = 0; i<numSamples; i++)
		{
			res.add(outputTokens[i]);
		}
		return res;
	}
	
	/*
	 * Representation of which input IDs make up a merged variant
	 */
	static class SimpleMergedVariant implements Comparable<SimpleMergedVariant>
	{
		static int index;
		ArrayList<String> ids;
		SimpleMergedVariant(ArrayList<String> dbSamples, ArrayList<String> dbIds)
		{
			String[] outputTokens = new String[sampleToIndex.size()];
			for(int i = 0; i<outputTokens.length; i++)
			{
				outputTokens[i] = "";
			}
			for(int i = 0; i < dbSamples.size(); i++)
			{
				int index = sampleToIndex.get(dbSamples.get(i));
				outputTokens[index] = dbIds.get(i);
			}
			ids = new ArrayList<String>();
			for(int i = 0; i<outputTokens.length; i++)
			{
				ids.add(outputTokens[i]);
			}
		}
		SimpleMergedVariant(String headerLine, String line)
		{
			ids = getSvPopIds(headerLine, line);
		}
		SimpleMergedVariant(VcfEntry entry) throws Exception
		{
			ids = getIds(entry);
		}
		static String makeHeader(ArrayList<String> filenames)
		{
			index = 0;
			StringBuilder res = new StringBuilder("");
			res.append("INDEX");
			for(String fn : filenames)
			{
				res.append("\t");
				res.append(fn);
			}
			return res.toString();
		}
		public String toString()
		{
			StringBuilder res = new StringBuilder("");
			res.append(index);
			index++;
			for(String id : ids)
			{
				res.append("\t");
				res.append(id.length() > 0 ? id : ".");
			}
			return res.toString();
		}
		public int compareTo(SimpleMergedVariant o)
		{
			for(int i = 0; i<ids.size(); i++)
			{
				String myId = ids.get(i), theirId = o.ids.get(i);
				if(myId.equals("") && !theirId.equals(""))
				{
					return 1;
				}
				if(!myId.equals("") && theirId.equals(""))
				{
					return -1;
				}
				if(!myId.equals(theirId))
				{
					return myId.compareTo(theirId);
				}
			}
			return 0;
		}
	}

}
