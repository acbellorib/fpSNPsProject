package fps;

import java.io.*;

public class RunJob_Antonio
{
	
	public int exitVal;
	
	public StreamCatcher errorCatcher;
	public StreamCatcher outputCatcher;
	
	public String [] runJob(String command)
	{		
		outputCatcher = null;

		try
		{		
			
			String[] cmd = null;
			Process proc = null;
			Runtime rt = Runtime.getRuntime();
			
			String os = getOsName();
			System.out.println("OS = " + os);
			if( os.startsWith("Windows"))
			{			
				cmd = new String[3];
				cmd[0] = "cmd.exe";
				cmd[1] = "/C";
				cmd[2] = command;
				
				proc = rt.exec(cmd);
			}
			else if(os.startsWith("Linux"))
			{
				cmd = new String[1];
				cmd[0] = command;
				
				proc = rt.exec(command.split(" "));
			}
			
			System.out.println("command is:");
			for (int i = 0; i < cmd.length; i++)
			{
				System.out.print(cmd[i] + " ");
			}
			System.out.println();
			
			// any error message?
			errorCatcher = new StreamCatcher(proc.getErrorStream(), "ERROR");            
			
			// any output?
			outputCatcher = new StreamCatcher(proc.getInputStream(), "OUTPUT");
			
			// kick them off
			errorCatcher.start();
			outputCatcher.start();
			
			// any error???
			exitVal = proc.waitFor();
			
			System.out.println("ExitValue: " + exitVal);   
					            
			if(exitVal != 0)
				System.out.println(errorCatcher.builder.toString());
		}
		catch (Throwable t)
		{
			t.printStackTrace();
		}
		
		//now return both the stdout and stderr streams
		//just use a String array for this -- not pretty but it works
		//stdout is at position [0], stderr is at [1]
		String [] outputs = new String [2];
		outputs[0] = outputCatcher.builder.toString();
		outputs[1] = errorCatcher.builder.toString();
		
		return outputs;
	}
	
	public static class StreamCatcher extends Thread
	{
		InputStream is;
		String type;
		public StringBuilder builder = new StringBuilder();
		
		StreamCatcher(InputStream is, String type)
		{
			this.is = is;
			this.type = type;
		}
		
		public void run()
		{
			try
			{
				InputStreamReader isr = new InputStreamReader(is);
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				while ((line = br.readLine()) != null)
					builder.append(line + "\n");
			}
			catch (IOException ioe)
			{
				ioe.printStackTrace();
			}
		}
	}
	
	public static String getOsName()
	{
		return System.getProperty("os.name"); 
	}
}
