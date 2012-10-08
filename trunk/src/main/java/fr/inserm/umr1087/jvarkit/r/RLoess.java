package fr.inserm.umr1087.jvarkit.r;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import fr.inserm.umr1087.jvarkit.collections.DefaultListFactory;
import fr.inserm.umr1087.jvarkit.collections.Listfactory;

public class RLoess extends AbstractRCmd
	{
	private boolean init=false;
	private int count=0;
	private Listfactory<Double> listOfDoubleFactory=new DefaultListFactory<Double>();
	
	public RLoess()
		{
		}
	
	public Listfactory<Double> getListOfDoubleFactory()
		{
		return listOfDoubleFactory;
		}
	
	public void setListOfDoubleFactory(Listfactory<Double> listOfDoubleFactory)
		{
		this.listOfDoubleFactory = listOfDoubleFactory;
		}
	
	public void add(double x,double y) throws IOException
		{
		PrintStream out=this.getOutputStream();
		if(!init)
			{
			init=true;
			out.print("T<-as.data.frame(matrix(c(");
			out.print(x);
			out.print(',');
			out.print(y);
			}
		else
			{
			out.print(',');
			out.print(x);
			out.print(',');
			out.print(y);
			}
		++count;
		}
	public List<Double> smooth() throws IOException,InterruptedException
		{
		if(!init || count==0) throw new RuntimeException();
		PrintStream out=this.getOutputStream();
		BufferedReader stdin=this.getReader();
		out.println("),ncol=2,byrow=TRUE))");
		out.println("colnames(T)<-c('x','y')");
		out.println("T2<-loess(y ~ x, T)");
		out.println("write.table(residuals(T2),'',col.names= F,row.names=F,sep='\\t')");
		out.flush();
		out.close();
		List<Double> y=this.listOfDoubleFactory.createList(count);
		for(int i=0;i< count;++i)
			{
			y.add( Double.parseDouble(stdin.readLine()) );
			}
		stdin.close();
		super.proc.waitFor();
		checkExitStatus();
		return y;
		}
	public static void main(String[] args) throws Exception
		{
		for(int i=0;i< 10;++i)
			{
			RLoess rloess=new RLoess();
			for(int j=-10;j<10;++j) rloess.add(j, j);
			System.err.println(rloess.smooth());
			}
		
		}
	}
