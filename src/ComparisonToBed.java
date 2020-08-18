import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class ComparisonToBed {
	
	static String idFile = "";
	static String ofn = "";
	static String vcfFilelist = "";
	static int PADDING = 100;
	static boolean PRECISE = false;
	static boolean SPECIFIC = false;
	static boolean FIRST_PRECISE = false;
	static boolean FIRST_SPECIFIC = false;
	
	static void parseArgs(String[] args)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf('=');
			if(equalsIdx == -1)
			{
				if(arg.toLowerCase().endsWith("precise"))
				{
					PRECISE = true;
				}
				else if(arg.toLowerCase().endsWith("specific"))
				{
					SPECIFIC = true;
				}
				if(arg.toLowerCase().endsWith("fp"))
				{
					FIRST_PRECISE = true;
				}
				else if(arg.toLowerCase().endsWith("fs"))
				{
					FIRST_SPECIFIC = true;
				}
			}
			else
			{
				String key = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("id_file"))
				{
					idFile = val;
				}
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
				else if(key.equalsIgnoreCase("vcf_filelist"))
				{
					vcfFilelist = val;
				}
			}
			
		}
		
		if(idFile.length() == 0 || ofn.length() == 0 || vcfFilelist.length() == 0)
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
		System.out.println("Usage: java -cp src ComparisonToBed [args]");
		System.out.println("Required args:");
		System.out.println("  id_file      (String) - table listing samples and IDs of variant pairs of interest");
		System.out.println("  out_file     (String) - where to write BED file of regions of interest");
		System.out.println("  vcf_filelist (String) - txt file listing the merged VCFs");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  --specific - only include pairs where at least one variant has IS_SPECIFIC INFO field set to 1");
		System.out.println("  --precise  - only include pairs where at least one variant has PRECISE as an INFO field");
		System.out.println("  --fs       - only include pairs where the first variant has IS_SPECIFIC INFO field set to 1");
		System.out.println("  --fp       - only include pairs where the first variant has PRECISE as an INFO field");

		System.out.println();
	}
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		HashMap<String, ArrayList<VariantRegion>> idMap = new HashMap<String, ArrayList<VariantRegion>>();
		ArrayList<String> vcfFiles = PipelineManager.getFilesFromList(vcfFilelist);
		for(int i = 0; i<vcfFiles.size(); i++)
		{
			Scanner input = new Scanner(new FileInputStream(new File(vcfFiles.get(i))));
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.length() == 0 || line.startsWith("#"))
				{
					continue;
				}
				ArrayList<VariantRegion> fromLine = VariantRegion.fromLine(line);
				if(fromLine.size() > 0)
				{
					String id = (i+1) + "_" + fromLine.get(0).id;
					idMap.put(id, fromLine);
				}
			}
			input.close();
		}
		Scanner input = new Scanner(new FileInputStream(new File(idFile)));
		PrintWriter out = new PrintWriter(new File(ofn));
		input.nextLine();
		while(input.hasNext())
		{
			String line = input.nextLine();
			String[] tokens = line.split("\t");
			String id1 = tokens[0] + "_" + tokens[1];
			String id2 = tokens[2] + "_" + tokens[3];
			if(idMap.containsKey(id1) && idMap.containsKey(id2))
			{
				ArrayList<VariantRegion> regions1 = idMap.get(id1);
				ArrayList<VariantRegion> regions2 = idMap.get(id2);
				for(VariantRegion r1 : regions1)
				{
					for(VariantRegion r2 : regions2)
					{
						if(SPECIFIC && !r1.specific && !r2.specific)
						{
							continue;
						}
						if(PRECISE && !r1.precise && !r2.precise)
						{
							continue;
						}
						if(FIRST_SPECIFIC && !r1.specific)
						{
							continue;
						}
						if(FIRST_PRECISE && !r1.precise)
						{
							continue;
						}
						if(r1.chr.equals(r2.chr))
						{
							long minStart = Math.min(r1.start, r2.start);
							long maxEnd = Math.max(r1.end, r2.end);
							if(maxEnd <= 100000 + minStart)
							{
								out.println(r1.chr + "\t" + minStart + "\t" + maxEnd + "\t" + id1 + ";" + id2);
							}
						}
					}
				}
			}
			
		}
		input.close();
		out.close();
	}
	
	static class VariantRegion
	{
		String chr;
		long start, end;
		String id;
		boolean specific;
		boolean precise;
		
		VariantRegion(String chr, long start, long end, String id, boolean specific, boolean precise)
		{
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.id = id;
			this.specific = specific;
			this.precise = precise;
		}
		static ArrayList<VariantRegion> fromLine(String line) throws Exception
		{
			
			ArrayList<VariantRegion> res = new ArrayList<VariantRegion>();
			VcfEntry entry = new VcfEntry(line);
			String id = entry.getId();
			String type = entry.getNormalizedType();
			
			boolean specific = entry.hasInfoField("IS_SPECIFIC") && entry.getInfo("IS_SPECIFIC").equals("1");
			boolean precise = entry.tabTokens[7].contains(";PRECISE;") || entry.tabTokens[7].startsWith("PRECISE;");

			if(type.equalsIgnoreCase("TRA"))
			{
				String chr = entry.getChromosome();
				String chr2 = entry.getChr2();
				res.add(new VariantRegion(chr, Math.max(1, entry.getPos() - PADDING), entry.getPos() + PADDING, id, specific, precise));
				res.add(new VariantRegion(chr2, Math.max(1, entry.getEnd() - PADDING), entry.getEnd() + PADDING, id, specific, precise));
				return res;
			}
			
			long start = entry.getPos() - PADDING;				
			long end = entry.getEnd() + PADDING;
			
			// Avoid giving non-positive coords
			start = Math.max(start, 1);
			end = Math.max(end, 1);
			
			// Make sure entire insertion is covered
			if(entry.getNormalizedType().equals("INS"))
			{
				end = start + entry.getLength() + PADDING;
			}
			
			res.add(new VariantRegion(entry.getChromosome(), start, end, id, specific, precise));
			
			return res;
		}
	}
}
