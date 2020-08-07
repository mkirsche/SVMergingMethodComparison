/*
 * Builds a TSV giving the IDs making up each merged variant in a VCF
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;

public class BuildMergingTable
{
	static String vcfFn = "";
	static String vcfFilelist = "";
	static String ofn = "";
	static String mode = "jasmine";
	
	static HashMap<String, Integer> sampleToIndex;
	
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
		
		if(!mode.equals("jasmine") && !mode.equals("survivor") && !mode.equals("svtools") && !mode.equals("svimmer"))
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
		System.out.println("  mode (String) [jasmine] - the merging software used: jasmine, survivor, svtools, or svimmer");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		if(mode.equals("svtools") || mode.equals("svimmer"))
		{		
			buildSampleMap();
		}
		
		ArrayList<String> filenames = PipelineManager.getFilesFromList(vcfFilelist);
		
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		ArrayList<SimpleMergedVariant> variants = new ArrayList<SimpleMergedVariant>();
		
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			variants.add(new SimpleMergedVariant(new VcfEntry(line)));
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
					res.add(entry.tabTokens[9+i].split(":")[7]);
				}
			}
			
			return res;
		}
		
		else if(mode.equals("svtools"))
		{
			String sname = entry.getInfo("SNAME");
			String[] tokens = sname.split(",");
			String[] answers = new String[sampleToIndex.size()];
			Arrays.fill(answers, "");
			for(String token : tokens)
			{
				int colonIdx = token.indexOf(':');
				String sample = token.substring(0, colonIdx);
				String id = token.substring(1 + colonIdx);
				
				answers[sampleToIndex.get(sample)] = id;
			}
			
			ArrayList<String> res = new ArrayList<String>();
			for(String answer : answers) res.add(answer);
			return res;
		}
		
		else if(mode.equals("svimmer"))
		{
			String sname = entry.getInfo("MERGED_IDS");
			String[] tokens = sname.split(",");
			String[] answers = new String[sampleToIndex.size()];
			Arrays.fill(answers, "");
			for(String token : tokens)
			{
				int colonIdx = token.indexOf(':');
				String sample = token.substring(0, colonIdx);
				String id = token.substring(1 + colonIdx);
				
				answers[sampleToIndex.get(sample)] = id;
			}
			
			ArrayList<String> res = new ArrayList<String>();
			for(String answer : answers) res.add(answer);
			return res;
		}
		
		return new ArrayList<String>();
	}
	
	/*
	 * Representation of which input IDs make up a merged variant
	 */
	static class SimpleMergedVariant implements Comparable<SimpleMergedVariant>
	{
		static int index;
		ArrayList<String> ids;
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
