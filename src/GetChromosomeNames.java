import java.io.File;
import java.io.FileInputStream;
import java.util.Scanner;
import java.util.TreeSet;

public class GetChromosomeNames {
public static void main(String[] args) throws Exception
{
	String fn = args[0];
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	TreeSet<String> chrs = new TreeSet<String>();
	while(input.hasNext())
	{
		String vcfFn = input.nextLine();
		if(vcfFn.length() == 0)
		{
			continue;
		}
		Scanner vcfInput = new Scanner(new FileInputStream(new File(vcfFn)));
		while(vcfInput.hasNext())
		{
			String line = vcfInput.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			String chr = line.split("\t")[0];
			chrs.add(chr);
		}
		vcfInput.close();
	}
	input.close();
	
	StringBuilder res = new StringBuilder("");
	for(String chr : chrs)
	{
		if(res.length() > 0) res.append(" ");
		res.append(chr);
	}
	System.out.println(res.toString());
}
}
