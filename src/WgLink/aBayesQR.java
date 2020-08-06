package virushaplo;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import java.io.IOException;
import java.io.InputStream;

class ThreadF_aBayesQR implements Runnable {

	private String [] path;
	private String region;

	public ThreadF_aBayesQR(String [] path, String region) {
		this.path = path.clone();
		this.region= region;
	}

	public synchronized void setLog(String log, String r) {
//		System.out.println(r+"\t"+log);
		if (aBayesQR.log_dic.containsKey(r)) {
				aBayesQR.log_dic.put(r, aBayesQR.log_dic.get(r)+"\n"+log);
			}else {
				aBayesQR.log_dic.put(r, log);
		}
	}

	@Override
	public void run() {
		try {

	        Process process = Runtime.getRuntime().exec(this.path);

	        InputStream inputStream2 = process.getInputStream();
	        BufferedReader br2 = new BufferedReader(new InputStreamReader(
					inputStream2, "gb2312"));
			String line2 = null;
			while ((line2 = br2.readLine()) != null) {
				String ss =line2;
			}

	        InputStream inputStream = process.getErrorStream();
			BufferedReader br = new BufferedReader(new InputStreamReader(
					inputStream, "gb2312"));
			String line = null;
			while ((line = br.readLine()) != null) {
				setLog(line, this.region);
				
			}
//			int value = process.waitFor();
//			br.close();
			System.out.println("Region\t"+ this.region + " aBayesQR finished.");
			Thread.sleep( 10);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

public class aBayesQR {
	public String [] cmdArray;
	public String [] regionArray;
	public static String[] log_out_Arry ;
	public static String[] log_region_Arry ;
	public int N;
	public static int index =0 ;
	public static HashMap<String, String> log_dic;;
	
	public aBayesQR(String [] cmds, String [] regions) throws IOException {
		
		this.N= cmds.length*10;
		this.log_region_Arry = new String [N];
		this.log_out_Arry = new String [N];
		this.cmdArray= cmds.clone();
		this.regionArray= regions.clone();
		
	}
	
	public void run() throws IOException {
		
//		ExecutorService cachedThreadPool = Executors.newCachedThreadPool();
        ExecutorService exec = Executors.newFixedThreadPool(this.regionArray.length);
		
		this.log_dic = new HashMap<String,String>();

		for (int i = 0; i < this.cmdArray.length; i++) {
//		for (int i = this.cmdArray.length-1; i>=0  ; i--) {
//			System.out.println(this.cmdArray[i] );
			String [] tmp = this.cmdArray[i].split(" "); 
			exec.execute(new ThreadF_aBayesQR(tmp, this.regionArray[i]));
		}
//		cachedThreadPool.shutdown();
		exec.shutdown();
		while (true) {
//			if (cachedThreadPool.isTerminated()) {
			if (exec.isTerminated()) {
				for (String key : this.log_dic.keySet()) {
					if (key!=null) {
//						System.out.println("????\t"+key);
						if ((!this.log_dic.get(key).equals(null))  &&  (!this.log_dic.get(key).equals(""))) {
							System.out.println("aBayesQR Prompt the Message on Region "+key+":");
							System.out.println(this.log_dic.get(key) );
						}
					}
				}
				break;
			}
		}
		
	}
}
