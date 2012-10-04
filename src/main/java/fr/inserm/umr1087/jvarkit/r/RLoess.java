package fr.inserm.umr1087.jvarkit.r;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class RLoess extends AbstractRCmd
	{
	private boolean init=false;
	private int count=0;
	public RLoess()
		{
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
	public List<Double> smooth() throws IOException
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
		List<Double> y=new ArrayList<Double>(count);
		for(int i=0;i< count;++i)
			{
			y.add( Double.parseDouble(stdin.readLine()) );
			}
		stdin.close();
		return y;
		}
	}
