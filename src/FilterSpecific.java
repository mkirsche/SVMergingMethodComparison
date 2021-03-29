import java.util.*;
import java.io.*;
public class FilterSpecific {
public static void main(String[] args) throws Exception
{
	String fn = "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasmine_support.txt";
	String augmentedtable = "/home/mkirsche/jasmine_data/figures/figure5/all4_pop.jasmine_augmented.txt";
	String ofn = "/home/mkirsche/jasmine_data/figures/figure5/all4_pop_jasmine_markedspecprec_support.txt";
	
	Scanner input = new Scanner(new FileInputStream(new File(augmentedtable)));
	String[] header = input.nextLine().split("\t");
	int specIndex = -1, precIndex = -1;
	for(int i = 0; i<header.length; i++)
	{
		if(header[i].equals("SPECIFIC_FLAG")) specIndex = i;
		else if(header[i].equals("PRECISE_FLAG")) precIndex = i;
	}
	HashMap<String, HashSet<String>> specificIds = new HashMap<String, HashSet<String>>();
	HashMap<String, HashSet<String>> preciseIds = new HashMap<String, HashSet<String>>();
	
	System.err.println("Reading augmented table to find specific/precise IDs");
	
	int lineNum = 0;

	while(input.hasNext())
	{
		if(lineNum > 0 && lineNum%5000 == 0)
		{
			System.err.println("Done reading " + lineNum + " lines");
		}
		lineNum++;
		String[] tokens = input.nextLine().split("\t");
		if(tokens.length == 0) continue;
		
		boolean isSpecific = tokens[specIndex].equals("1");
		boolean isPrecise = tokens[precIndex].equals("1");
		
		if(!isSpecific && !isPrecise) continue;
		
		for(int i = 0; i<header.length; i++)
		{
			if(header[i].endsWith(".vcf"))
			{
				String name = header[i];
				if(name.contains("/"))  name = name.substring(1 + name.lastIndexOf('/'));
				if(name.contains("vGRCh38"))
				{
					name = name.substring(0, name.indexOf("vGRCh38"));
				}
				
				if(!specificIds.containsKey(name))
				{
					specificIds.put(name, new HashSet<String>());
					preciseIds.put(name, new HashSet<String>());
				}
				
				String[] ids = tokens[i].split(";");
				
				for(String id : ids)
				{
					if(isSpecific)
					{
						specificIds.get(name).add(id);
					}
					if(isPrecise)
					{
						preciseIds.get(name).add(id);
					}
				}
			}
		}
		
	}
	
	input.close();
	
	input = new Scanner(new FileInputStream(new File(fn)));
	PrintWriter out = new PrintWriter(new File(ofn));
		
	String headerLine = input.nextLine();
	header = headerLine.split("\t");
	int idIndex = -1;
	int sampleIndex = -1;
	for(int i = 0; i<header.length; i++)
	{
		if(header[i].equals("ID"))
		{
			idIndex = i;
		}
		else if(header[i].equals("SAMPLE"))
		{
			sampleIndex = i;
		}
	}
	out.println(headerLine + "\tSPECIFIC_FLAG\tPRECISE_FLAG");
	
	System.err.println("Adding annotations to table");
	
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split("\t");
		String sample = tokens[sampleIndex];
		String id = tokens[idIndex];
		String specificVal = "0", preciseVal = "0";
		if(specificIds.containsKey(sample) && specificIds.get(sample).contains(id))
		{
			specificVal = "1";
		}
		if(preciseIds.containsKey(sample) && preciseIds.get(sample).contains(id))
		{
			preciseVal = "1";
		}
		out.println(line + "\t" + specificVal + "\t" + preciseVal);
	}
}
}
