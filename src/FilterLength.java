import java.util.*;
import java.io.*;
public class FilterLength {
public static void main(String[] args) throws Exception
{
	String fn = args[0];
	int length = Integer.parseInt(args[1]);
	String ofn = fn + "_lenfilter" + length + ".vcf";
	if(args.length > 2)
	{
		ofn = args[2];
	}
	
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	PrintWriter out = new PrintWriter(new File(ofn));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.startsWith("#"))
		{
			out.println(line);
			continue;
		}
		VcfEntry entry = new VcfEntry(line);
		if(entry.getNormalizedType().equals("TRA") || Math.abs(entry.getLength()) >= length)
		{
			out.println(line);
		}
	}
	out.close();
}
}
