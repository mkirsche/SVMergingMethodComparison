import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class CompareMerges {
	static String[] tables;	
	static String outPrefix = "";
	static String lonePrefix = "";
	/*
	 * Parse command line arguments
	 */
	static void parseArgs(String[] args)
	{
		tables = new String[2];
		tables[0] = tables[1] = "";
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
					tables[0] = val;
				}
				else if(key.equalsIgnoreCase("file2"))
				{
					tables[1] = val;
				}
				else if(key.equalsIgnoreCase("files"))
				{
					tables = val.split(",");
				}
				else if(key.equalsIgnoreCase("out_prefix"))
				{
					outPrefix = val;
				}
				else if(key.equalsIgnoreCase("lone_prefix"))
				{
					lonePrefix = val;
				}
			}
			
		}
		
		if(tables.length < 2 || tables[0].length() == 0 || tables[1].length() == 0)
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
		System.out.println("  files   (String) - comma-separated list of table files in place of file1 and file2 if there are more than two");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  out_prefix   (String) - file prefix for optionally outputting lists of joined variants for each subset of merging results");
		System.out.println("  lone_prefix  (String) - file prefix for optionally outputting lists of unjoined variants for each subset of merging results");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		compare();
	}
	
	static void compare() throws Exception
	{
		ResultsSet rp = new ResultsSet(tables);

		rp.print();
	}
	
	static class ResultsSet
	{
		HashMap<String, Integer> intersampleCountByCallerSupport;
		HashMap<String, Integer> loneCountByCallerSupport;
		HashMap<String, PrintWriter> outByCallerSupport;
		HashMap<String, PrintWriter> loneOutByCallerSupport;
		int[] interCounts;
		int[] intraCounts;
		ResultsSet(String[] files) throws Exception
		{
			intersampleCountByCallerSupport = new HashMap<String, Integer>();
			loneCountByCallerSupport = new HashMap<String, Integer>();

			outByCallerSupport = new HashMap<String, PrintWriter>();
			loneOutByCallerSupport = new HashMap<String, PrintWriter>();
			int n = files.length;
			interCounts = new int[n];
			intraCounts= new int[n];
			AdjacencyList[] graphs = new AdjacencyList[n];
			for(int i = 0; i<n; i++)
			{
				graphs[i] = new AdjacencyList(files[i]);
				interCounts[i] = graphs[i].intersample;
				intraCounts[i] = graphs[i].intrasample;
			}
			
			for(int i = 0; i<n; i++)
			{
				ArrayList<String[]> edges = graphs[i].allEdges();
				
				for(String[] edge : edges)
				{
					StringBuilder suppVec = new StringBuilder("");
					
					boolean firstSample = true;
					for(int j = 0; j<i; j++)
					{
						if(graphs[j].hasEdge(edge[0], edge[1]))
						{
							firstSample = false;
						}
						else suppVec.append("0");
					}
					
					if(!firstSample)
					{
						continue;
					}
					
					suppVec.append("1");
					
					for(int j = i+1; j<n; j++)
					{
						if(graphs[j].hasEdge(edge[0], edge[1]))
						{
							suppVec.append("1");
						}
						else suppVec.append("0");
					}
					
					String svKey = suppVec.toString();
					
					if(!intersampleCountByCallerSupport.containsKey(svKey))
					{
						if(outPrefix.length() > 0)
						{
							outByCallerSupport.put(svKey, new PrintWriter(new File(outPrefix + "_" + svKey + ".txt")));
							outByCallerSupport.get(svKey).println("SAMPLE1\tID1\tSAMPLE2\tID2");
						}
						intersampleCountByCallerSupport.put(svKey, 1);
					}
					else
					{
						intersampleCountByCallerSupport.put(svKey, 1 + intersampleCountByCallerSupport.get(svKey));
					}
					
					if(outByCallerSupport.containsKey(svKey))
					{
						String[] id1 = AdjacencyList.decode(edge[0]);
						String[] id2 = AdjacencyList.decode(edge[1]);
						outByCallerSupport.get(svKey).println(id1[0]+"\t"+id1[1]+"\t"+id2[0]+"\t"+id2[1]);
					}
				}
			}
			
			for(String s : outByCallerSupport.keySet())
			{
				outByCallerSupport.get(s).close();
			}
			
			for(int i = 0; i<n; i++)
			{
				HashMap<String, String> joinedIds = graphs[i].joinedIds;
				
				for(String id : joinedIds.keySet())
				{
					StringBuilder suppVec = new StringBuilder("");
					
					boolean firstSample = true;
					for(int j = 0; j<i; j++)
					{
						if(graphs[j].joinedIds.containsKey(id))
						{
							firstSample = false;
						}
						else suppVec.append("1");
					}
					
					if(!firstSample)
					{
						continue;
					}
					
					suppVec.append("0");
					
					for(int j = i+1; j<n; j++)
					{
						if(graphs[j].joinedIds.containsKey(id))
						{
							suppVec.append("0");
						}
						else suppVec.append("1");
					}
					
					String svKey = suppVec.toString();
					
					if(!loneOutByCallerSupport.containsKey(svKey) && lonePrefix.length() > 0)
					{
						loneOutByCallerSupport.put(svKey, new PrintWriter(new File(lonePrefix + "_" + svKey + ".txt")));
						loneOutByCallerSupport.get(svKey).println("SAMPLE1\tID1\tSAMPLE2\tID2");
						loneCountByCallerSupport.put(svKey, 1);
					}
					else
					{
						loneCountByCallerSupport.put(svKey, 1 + loneCountByCallerSupport.get(svKey));
					}
					
					if(loneOutByCallerSupport.containsKey(svKey))
					{
						String[] loneId = AdjacencyList.decode(id);
						String[] neighborId = AdjacencyList.decode(joinedIds.get(id));
						loneOutByCallerSupport.get(svKey).println(loneId[0] + "\t" + loneId[1] + "\t" + neighborId[0] + "\t" + neighborId[1]);
					}
				}
			}
			
			for(String s : loneOutByCallerSupport.keySet())
			{
				loneOutByCallerSupport.get(s).close();
			}
		}
		
		public void print()
		{
			for(int i = 0; i<intraCounts.length; i++)
			{
				System.out.println("Caller " + (i+1) + " intrasample: " + intraCounts[i]);
				System.out.println("Caller " + (i+1) + " intersample: " + interCounts[i]);
			}
			for(String s : intersampleCountByCallerSupport.keySet())
			{
				System.out.println("Merged pairs with caller vector " + s + ":" + intersampleCountByCallerSupport.get(s));
			}
			for(String s : loneCountByCallerSupport.keySet())
			{
				System.out.println("Unmerged variants with caller vector " + s + ":" + loneCountByCallerSupport.get(s));
			}
		}
	}
	
	static class AdjacencyList
	{
		HashMap<String, Integer> idToNodeIndex;
		HashMap<Integer, String> nodeIndexToId;
		ArrayList<HashSet<Integer>> graph;
		HashMap<String, String> joinedIds;
		String[] header;
		int intrasample = 0;
		int intersample = 0;
		AdjacencyList(String filename) throws Exception
		{
			joinedIds = new HashMap<String, String>();
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
					if(tokens[i].equals("."))
					{
						continue;
					}
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
							if(tokens[j].equals("."))
							{
								continue;
							}
							String[] idsFromSampleJ = tokens[j].split(";");
							intersample += idsFromSampleJ.length;
							for(String id2 : idsFromSampleJ)
							{
								String key2 = j + "_" + id2;
								int node2 = idToNodeIndex.get(key2);
								graph.get(node2).add(node1);
								
								if(!joinedIds.containsKey(key1))
								{
									joinedIds.put(key1, key2);
								}
								if(!joinedIds.containsKey(key2))
								{
									joinedIds.put(key2, key1);
								}
							}
						}
					}
				}
			}
			input.close();
		}
		
		static String[] decode(String id)
		{
			int underscoreIdx = id.indexOf('_');
			return new String[] { id.substring(0, underscoreIdx), id.substring(1 + underscoreIdx) };
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
