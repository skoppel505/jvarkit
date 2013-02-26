package fr.inserm.umr1087.jvarkit.tools.vcftabix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;


import net.sf.samtools.tabix.TabixReader;
import net.sf.samtools.util.BlockCompressedInputStream;


public class VCFTabix
	{
	private PrintStream out=System.out;
	private String tabixFile;
	private TabixReader tabixReader =null;
	private Set<String> infoIds=new LinkedHashSet<String>();
	private void run(BufferedReader in) throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			
			if(line.startsWith("#"))
				{
				if(line.startsWith("#CHROM"))
					{
					out.println("##Annotated with "+getClass()+":"+tabixFile);
					BlockCompressedInputStream src=new BlockCompressedInputStream(new File(tabixFile));
					BufferedReader r2=new BufferedReader(new InputStreamReader(src));
					String line2;
					while((line2=r2.readLine())!=null)
						{
						if(line2.startsWith("#CHROM")) break;
						if(!line2.startsWith("##INFO=")) continue;
						int i=line2.indexOf("ID=");
						if(i==-1) continue;
						int j=line2.indexOf(i+1,',');
						if(j==-1) j=line2.indexOf(i+1,'>');
						if(j==-1) continue;
						if(this.infoIds.contains(line2.substring(i+3,j)))
							{
							out.println(line2);
							}
						} 
					r2.close();
					out.println(line);
					continue;
					}
				out.println(line);
				continue;
				}
			
			String tokens[]=tab.split(line,9);
			if(tokens.length<8)
				{
				System.err.println("Error not enought columns in "+line);
				continue;
				}
			String chrom=tokens[0];
			Integer pos1=Integer.parseInt(tokens[1]);
			
		
			
			TabixReader.Iterator iter=tabixReader.query(chrom+":"+pos1);
			String line2;
			Map<String,String> map=new LinkedHashMap<String,String>(this.infoIds.size());
			while((line2=iter.next())!=null)
				{
				String tokens2[]=tab.split(line2,9);
				if(tokens2.length<8)
					{
					System.err.println("Error not enought columns in "+line2);
					continue;
					}
				if(!tokens[0].equals(tokens2[0])) continue;
				if(!tokens[1].equals(tokens2[1])) continue;
				if(!tokens[3].equalsIgnoreCase(tokens2[3])) continue;
				if(!tokens[4].equalsIgnoreCase(tokens2[4])) continue;
				if((tokens[2].isEmpty() || tokens[2].equals(".")) && !(tokens2[2].isEmpty() || tokens2[2].equals(".")))
					{
					tokens[2]=tokens2[2];
					}
				
				final String infos=tokens2[7];
				for(String id:this.infoIds)
					{
					int i=-1;
					for(;;)
						{
						i=infos.indexOf(id,i+1);
						if(i==-1) break;
						if(!(i==0 || infos.charAt(i-1)!=';')) continue;
						int j=i+id.length();
						if(j>=infos.length() || infos.charAt(j)!='=') continue;
						i=j+1;
						j=infos.indexOf(';',i);
						if(j==-1) j=infos.length();
						map.put(id, infos.substring(i,j));
 						}
					}
				}
			StringBuilder b=new StringBuilder();
			for(String k:map.keySet())
				{
				b.append(k);
				b.append("=");
				b.append(map.get(k));
				b.append(";");
				}
			
			
			String info=tokens[7];
			if(info.equals(".") || info.isEmpty())
				{
				info=b.toString();
				}
			else
				{
				info=b.toString()+info;
				}
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print('\t');
				out.print(i==7?info:tokens[i]);
				}
			out.println();
			}
		
		}
	public int run(String[] args) throws IOException
		{
		int optind=0;
		while(optind<args.length)
			{
			if(args[optind].equals("-h"))
				{
				System.out.println(" -f (tabix-file) required.");
				System.out.println(" -i (String) VCF-INFO-ID optional can be used several times.");
				return 0;
				}
			else if(args[optind].equals("-i") && optind+1< args.length)
				{
				this.infoIds.add(args[++optind]);
				}
			else if(args[optind].equals("-f") && optind+1< args.length)
				{
				this.tabixFile=args[++optind];
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unnown option: "+args[optind]);
				return -1;
				}
			else
				{
				break;
				}
			++optind;
			}
		if(tabixFile==null)
			{
			System.err.println("Undefined tabix File");
			return -1;
			}
		
		this.tabixReader=new TabixReader(this.tabixFile);
		
		if(optind==args.length)
			{
			this.run(new BufferedReader(new InputStreamReader(System.in)));
			}
		else if(optind+1==args.length)
			{
			String inputName=args[optind++];
			BufferedReader in=new BufferedReader(new FileReader(inputName));
			this.run(in);
			in.close();
			}
		else
			{
			System.err.println("Illegal Number of arguments");
			return -1;
			}
		return 0;
		}
	
	public static void main(String[] args) throws IOException
		{
		VCFTabix app=new VCFTabix();
		app.run(args);
		}
}
