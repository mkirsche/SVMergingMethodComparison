import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class AddReadSupport
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
		System.out.println("Usage: java -cp src AddReadSupport [args]");
		System.out.println("Required args:");
		System.out.println("  table_file   (String) - the TSV built from the BuildMergingTable script");
		System.out.println("  out_file     (String) - the file to output the table of read support vs. allele frequency to");
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
						
		// For each input file, a map from variant ID to the index in the merged variant list to which it belongs
		HashMap<String, Integer>[] idToAFPerSample = new HashMap[filenames.size()];
		
		for(int i = 0; i<filenames.size(); i++)
		{
			idToAFPerSample[i] = new HashMap<String, Integer>();
		}
		
		input.nextLine(); // TODO get the filenames from here instead of filelist
				
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			AugmentMergingTable.DetailedMergedVariant v = new AugmentMergingTable.DetailedMergedVariant(line, filenames);
			
			int af = 0;
			for(int i = 0; i<v.ids.size(); i++)
			{
				if(v.ids.get(i).length() > 0)
				{
					af++;
				}
			}
			for(int i = 0; i<v.ids.size(); i++)
			{
				if(v.ids.get(i).length() > 0)
				{
					String[] ids = v.ids.get(i).split(";");
					for(String id : ids)
					{
						idToAFPerSample[i].put(id, af);
					}
				}
			}
		}
		
		input.close();
		
		out.println("SAMPLE\tTECH\tID\tAF\tRE");
		
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
				if(idToAFPerSample[i].containsKey(id))
				{
					String name = filenames.get(i);
					name = name.substring(1 + name.lastIndexOf('/'));
					if(name.contains("vGRCh38"))
					{
						name = name.substring(0, name.indexOf("vGRCh38"));
					}
					String tech = "Hifi";
					if(filenames.get(i).contains("_PB_")) tech = "CLR";
					if(filenames.get(i).contains("_ONT_")) tech = "ONT";
					int support = entry.getReadSupport();
					out.println(name + "\t" + tech + "\t" + id + "\t" + idToAFPerSample[i].get(id) + "\t" + support);
				}
			}
		}
		
		
		out.close();
	}
}
