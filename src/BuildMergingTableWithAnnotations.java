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
import java.util.TreeMap;
import java.util.TreeSet;

public class BuildMergingTableWithAnnotations
{
	static String vcfFn = "";
	static String vcfFilelist = "";
	static String ofn = "";
	static String mode = "jasmine";
	static String sampleList = "";
	static String bedtoolsFile = "";
	
	static HashMap<String, Integer> sampleToIndex;
	static HashMap<String, String>[] coordsToId;
	
	static TreeSet<String> annotationOptions;
	static TreeMap<String, TreeSet<String>> idToAnnotations;
	
	static String annotationsOutput = ""; 
	
	static void extractAnnotations() throws Exception
	{
		annotationOptions = new TreeSet<String>();
		idToAnnotations = new TreeMap<String, TreeSet<String>>();
		
		if(annotationsOutput.length() > 0 && new File(annotationsOutput).exists())
		{
			Scanner input = new Scanner(new FileInputStream(new File(annotationsOutput)));
			String[] options = input.nextLine().split(" ");
			for(String s : options) annotationOptions.add(s);
			while(input.hasNext())
			{
				String line = input.nextLine();
				String[] tokens = line.split(" ");
				String id = tokens[0];
				TreeSet<String> vals = new TreeSet<String>();
				for(int i = 1; i<tokens.length; i++)
				{
					vals.add(tokens[i]);
				}
				idToAnnotations.put(id, vals);
			}
			input.close();
		}
		else
		{
			Scanner input = new Scanner(new FileInputStream(new File(bedtoolsFile)));
			while(input.hasNextLine())
			{
				String line = input.nextLine();
				String[] tokens = line.split("\t");
				String id = tokens[0];
				int idx = 1;
				
				String annotation = "";
				while(idx < tokens.length && !tokens[idx].equals("1") && !tokens[idx].equals("2") && !tokens[idx].equals("3"))
				{
					idx++;
				}
				if(idx == tokens.length)
				{
					continue;
				}
				
				if(tokens[idx].equals("1"))
				{
					annotation = "CENTROMERE";
				}
				else if(tokens[idx].equals("2"))
				{
					annotation = "REPEAT_" + tokens[idx+7];
				}
				else if(tokens[idx].equals("3"))
				{
					annotation = "GENE_" + tokens[idx+3];
				}
				
				if(annotation.length() == 0)
				{
					continue;
				}
				annotationOptions.add(annotation);
				if(!idToAnnotations.containsKey(id))
				{
					idToAnnotations.put(id, new TreeSet<String>());
				}
				idToAnnotations.get(id).add(annotation);
			}
			input.close();
			writeAnnotations();
		}
	}
	
	static void writeAnnotations() throws Exception
	{
		if(annotationsOutput.length() == 0)
		{
			return;
		}
		PrintWriter out = new PrintWriter(new File(annotationsOutput));
		StringBuilder firstLine = new StringBuilder("");
		for(String s : annotationOptions)
		{
			firstLine.append(s+" ");
		}
		out.println(firstLine.toString().trim());
		
		for(String s : idToAnnotations.keySet())
		{
			out.print(s);
			for(String v : idToAnnotations.get(s))
			{
				out.print(" " + v);
			}
			out.println();
		}
		out.close();
	}
	
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
				else if(key.equalsIgnoreCase("bed_file"))
				{
					bedtoolsFile = val;
				}
				else if(key.equalsIgnoreCase("annotation_file"))
				{
					annotationsOutput = val;
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
		
		if(vcfFn.length() == 0 || ofn.length() == 0 || vcfFilelist.length() == 0 || bedtoolsFile.length() == 0)
		{
			usage();
			System.exit(0);
		}
		
		if(!mode.equals("jasmine") && !mode.equals("survivor") && !mode.equals("svtools") && !mode.equals("svimmer") && !mode.equals("jasmine_intra")
				&& !mode.equals("svpop"))
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
		System.out.println("  bed_file     (String) - the bed file with the VCF intersected with centromeres, repeats, and genes");
		System.out.println("  out_file     (String) - the name of the file to output the merging table to");
		System.out.println("  vcf_filelist (String) - a txt file with each line containing the filename of a merged vcf");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  mode            (String) [jasmine] - the merging software used: jasmine, survivor, svtools, jasmine_intra, svpop, or svimmer");
		System.out.println("  annotation_file (String) []        - file to output annotations to avoid reading the large bed file every time");
		System.out.println("  sample_list     (String) []        - comma-separated list of sample names in same order as in table - required when using svpop");
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
		
		if(mode.equals("svtools") || mode.equals("svimmer"))
		{		
			buildSampleMap();
		}
		
		if(mode.equals("survivor"))
		{
			buildCoordsMap();
		}
		
		extractAnnotations();
		
		ArrayList<String> filenames = PipelineManager.getFilesFromList(vcfFilelist);
		
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		ArrayList<SimpleMergedVariant> variants = new ArrayList<SimpleMergedVariant>();
		
		String headerLine = "";
		
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
		TreeSet<String> annotations;
		SimpleMergedVariant(String headerLine, String line)
		{
			ids = getSvPopIds(headerLine, line);
			annotations = new TreeSet<String>();
		}
		SimpleMergedVariant(VcfEntry entry) throws Exception
		{
			ids = getIds(entry);
			if(idToAnnotations.containsKey(entry.getId()))
			{
				annotations = idToAnnotations.get(entry.getId());
			}
			else
			{
				annotations = new TreeSet<String>();
			}
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
			res.append("\t" + "ANNOTATIONS");
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
			res.append("\t");
			StringBuilder annotationString = new StringBuilder("");
			if(annotations.size() == 0)
			{
				annotationString.append(".");
			}
			for(String s : annotations)
			{
				if(annotationString.length() > 0)
				{
					annotationString.append(",");
				}
				annotationString.append(s);
			}
			res.append(annotationString);
			
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
