package fr.inserm.umr1087.jvarkit.io;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class ConsummeInputStreamThread extends Thread	
   {
   private InputStream is;
   private OutputStream os;
   public ConsummeInputStreamThread(InputStream is,OutputStream os)
    	{
        this.is = is;
        this.os=os;
    	}
   public ConsummeInputStreamThread(InputStream is)
	{
    this(is,null);
	}
  
    @Override
    public void run()
	    {
        try
	        {
        	if(os!=null)
        		{
        		int c;
     	        while ((c=is.read())!=-1)
     	        	{
     	        	os.write(c);
     	        	}
        		}
        	else
        		{
        		while (is.read()!=-1);
        		}
	        }
        catch (IOException ioe)
              {
        	  ioe.printStackTrace();  
              }
	    }
    
	}
