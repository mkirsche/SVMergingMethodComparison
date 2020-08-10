import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class CompareMerges {
	static String table1Fn = "";
	static String table2Fn = "";
		
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
				if(key.equalsIgnoreCase("file1"))
				{
					table1Fn = val;
				}
				else if(key.equalsIgnoreCase("file2"))
				{
					table2Fn = val;
				}
			}
			
		}
		
		if(table1Fn.length() == 0 || table2Fn.length() == 0)
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
		System.out.println("Usage: java -cp src CompareMerges [args]");
		System.out.println("Required args:");
		System.out.println("  file1   (String) - first table file (a TSV built from the BuildMergingTable script)");
		System.out.println("  file2   (String) - second table file");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		compare();
	}
	
	static void compare() throws Exception
	{
		ResultsPair rp = new ResultsPair(table1Fn, table2Fn);
		rp.print();
	}
	
	static class ResultsPair
	{
		int firstOnly = 0;
		int secondOnly = 0;
		int both = 0;
		int inter1 = 0;
		int inter2 = 0;
		int intra1 = 0;
		int intra2 = 0;
		double jaccard = 0.0;
		ResultsPair(String fn1, String fn2) throws Exception
		{
			AdjacencyList graph1 = new AdjacencyList(fn1), graph2 = new AdjacencyList(fn2);
			intra1 = graph1.intrasample;
			intra2 = graph2.intrasample;
			inter1 = graph1.intersample;
			inter2 = graph2.intersample;
			
			ArrayList<String[]> edges1 = graph1.allEdges(), edges2 = graph2.allEdges();
			
			for(String[] edge : edges1)
			{
				if(!graph2.hasEdge(edge[0], edge[1]))
				{
					firstOnly++;
				}
				else
				{
					both++;
				}
			}
			for(String[] edge : edges2)
			{
				if(!graph1.hasEdge(edge[0], edge[1]))
				{
					secondOnly++;
				}
			}
			
			jaccard = 1.0 * both / (both + firstOnly + secondOnly);
		}
		
		public void print()
		{
			System.out.println("Intersample 1: " + inter1);
			System.out.println("Intersample 2: " + inter2);
			System.out.println("Intersample 1 only: " + firstOnly);
			System.out.println("Intersample 2 only: " + secondOnly);
			System.out.println("Intersample both: " + both);
			System.out.println("Intrasample 1: " + intra1);
			System.out.println("Intrasample 2: " + intra2);
			
			System.out.println("Jaccard: " + jaccard);
		}
	}
	
	static class AdjacencyList
	{
		HashMap<String, Integer> idToNodeIndex;
		HashMap<Integer, String> nodeIndexToId;
		ArrayList<HashSet<Integer>> graph;
		String[] header;
		int intrasample = 0;
		int intersample = 0;
		AdjacencyList(String filename) throws Exception
		{
			idToNodeIndex = new HashMap<String, Integer>();
			nodeIndexToId = new HashMap<Integer, String>();
			graph = new ArrayList<HashSet<Integer>>();
			Scanner input = new Scanner(new FileInputStream(new File(filename)));
			header = input.nextLine().split("\t");
			while(input.hasNext())
			{
				String[] tokens = input.nextLine().split("\t");
				int n = tokens.length;
				for(int i = 1; i<n; i++)
				{
					String[] idsFromSampleI = tokens[i].split(";");
					intrasample += idsFromSampleI.length - 1;
					for(String id1 : idsFromSampleI)
					{
						String key1 = i + "_" + id1;
						if(!idToNodeIndex.containsKey(key1))
						{
							int nodeIndex = idToNodeIndex.size();
							idToNodeIndex.put(key1, nodeIndex);
							nodeIndexToId.put(nodeIndex, key1);
							graph.add(new HashSet<Integer>());
						}
						int node1 = idToNodeIndex.get(key1);
						for(int j = 1; j<i; j++)
						{
							String[] idsFromSampleJ = tokens[j].split(";");
							intersample += idsFromSampleJ.length;
							for(String id2 : idsFromSampleJ)
							{
								String key2 = j + "_" + id2;
								int node2 = idToNodeIndex.get(key2);
								graph.get(node2).add(node1);
							}
						}
					}
				}
			}
			input.close();
		}
		
		boolean hasEdge(String id1, String id2)
		{
			if(!idToNodeIndex.containsKey(id1) || !idToNodeIndex.containsKey(id2))
			{
				return false;
			}
			int node1 = idToNodeIndex.get(id1), node2 = idToNodeIndex.get(id2);
			return graph.get(node1).contains(node2);
		}
		
		ArrayList<String[]> allEdges()
		{
			ArrayList<String[]> res = new ArrayList<String[]>();
			for(int i = 0; i<graph.size(); i++)
			{
				String id1 = nodeIndexToId.get(i);
				for(int j : graph.get(i))
				{
					String id2 = nodeIndexToId.get(j);
					res.add(new String[] {id1, id2});
				}
			}
			return res;
		}
	}
	
	static class Edge
	{
		String from, to;
		String component;
	}
}
