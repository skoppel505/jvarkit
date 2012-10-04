package fr.inserm.umr1087.jvarkit.r;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;

import fr.inserm.umr1087.jvarkit.io.ConsummeInputStreamThread;

public abstract class AbstractRCmd
	{
	private String RScript="/usr/local/bin/Rscript";
	private File baseDir=null;
	private PrintStream outputStream;
	private BufferedReader inputStream;
	protected ConsummeInputStreamThread errStream;
	private Process proc=null;
	public AbstractRCmd()
		{
		this.baseDir=new File(System.getProperty("user.dir", "."));
		}

	private void start() throws IOException
		{
		if(this.proc!=null) return;
		this.proc=Runtime.getRuntime().exec(
				new String[]{RScript,"-"},
				new String[]{},
				baseDir
				);
		this.errStream=new ConsummeInputStreamThread(proc.getErrorStream());
		this.errStream.start();
		this.outputStream=new PrintStream(proc.getOutputStream());
		this.inputStream=new BufferedReader(new InputStreamReader(proc.getInputStream()));
		}
	public PrintStream getOutputStream() throws IOException
		{
		if(outputStream==null) start();
		return this.outputStream;
		}
	public BufferedReader getReader() throws IOException
		{
		if(inputStream==null) start();
		return this.inputStream;
		}
	}
