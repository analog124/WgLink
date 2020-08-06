package virushaplo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Properties;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

public class Simulations {
	public int num_haps;
	public String[] haplotypes;
	public double [] coverages;
	public String reference;
	public String dwgsim;
	public int outer_dist;
	public String out_dir;
	int read_len;
	double error_rate;
	public String project_name;
	public String gunzip;
	
	public Simulations(String parameter_file) throws IOException {
		
		InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
		this.out_dir = prop.getProperty("Fastq_Output_Path")+"/";
		this.dwgsim= prop.getProperty("DWGSIM");
		this.reference= prop.getProperty("Reference_Seq");
		this.error_rate= Double.parseDouble(prop.getProperty("Error_Rate_Per_Base"));
		this.project_name= prop.getProperty("Proj_Name");
//		this.gunzip= prop.getProperty("gunzip");
		this.read_len = Integer.parseInt(prop.getProperty("Read_Len"));
		this.outer_dist = Integer.parseInt(prop.getProperty("Outer_Dist"))+ this.read_len;
		String haps = prop.getProperty("Simulation_Seq");
		haps = haps.replace("\r", "").replace("\n", "").replace(" ", "");
		String[] tmp_hap  = haps.split("\t"); 
		this.num_haps= tmp_hap.length;
		this.haplotypes= new String [this.num_haps];
		for (int i=0;i < this.num_haps;i++) {
			this.haplotypes[i]= tmp_hap[i];
		}
		String covs = prop.getProperty("Simulation_Coverage");
		covs = covs.replace("\r", "").replace("\n", "").replace(" ", "");
		String[] tmp_cov  = covs.split("\t"); 
		this.coverages= new double [this.num_haps];
		for (int i=0;i < this.num_haps;i++) {
			this.coverages[i]=  Double.parseDouble( tmp_cov[i]  );
		}
		is.close();
	}
	
	public void generate_fastq()  throws IOException, InterruptedException {
		
		new File(String.valueOf(this.out_dir) ).mkdir();
		for (int i=0;i< this.num_haps;i++) {
			String fastq_prefix= this.out_dir+"tmp_"+Integer.toString(i);
			ProcessBuilder CMDLine = new ProcessBuilder(this.dwgsim,
   	                this.haplotypes[i], 
   	                fastq_prefix, 
   	                "-e", Double.toString(this.error_rate),
   	                "-E", Double.toString(this.error_rate), "-C", Double.toString(this.coverages[i]),
   	                "-1", Integer.toString(this.read_len),
   	                "-2", Integer.toString(this.read_len),
   	                "-r", "0",
   	                "-F", "0",
   	                "-H",
   	                "-d", Integer.toString(this.outer_dist),
   	                "-o", "1",
   	                "-s", "0",
   	                "-y", "0");
			
   	            Process CMDProcess = CMDLine.start();
   	            CMDProcess.waitFor();
   	            System.out.println("Finished simulating reads for haplotype:\t" + this.haplotypes[i] );
   	            
   	            byte[] buffer = new byte[1024];
   	            FileInputStream fileIn_1 = new FileInputStream(fastq_prefix+".bwa.read1.fastq.gz");
   	         
   	            GZIPInputStream gZIPInputStream_1 = new GZIPInputStream(fileIn_1);
   	 
   	            FileOutputStream fileOutputStream_1 = new FileOutputStream(fastq_prefix+".bwa.read1.fastq");
   	 
   	            int bytes_read;
   	 
   	            while ((bytes_read = gZIPInputStream_1.read(buffer)) > 0) {
   	 
   	                fileOutputStream_1.write(buffer, 0, bytes_read);
   	            }
   	 
   	            gZIPInputStream_1.close();
   	            fileOutputStream_1.close();
   	            FileInputStream fileIn_2= new FileInputStream(fastq_prefix+".bwa.read2.fastq.gz");
	            GZIPInputStream gZIPInputStream_2 = new GZIPInputStream(fileIn_2);
	            FileOutputStream fileOutputStream_2 = new FileOutputStream(fastq_prefix+".bwa.read2.fastq");
	            while ((bytes_read = gZIPInputStream_2.read(buffer)) > 0) {
	 
	                fileOutputStream_2.write(buffer, 0, bytes_read);
	            }
	            gZIPInputStream_2.close();
	            fileOutputStream_2.close();
   	            System.out.println("The file was decompressed successfully!");
		}
		
		String output_fastq_1 = this.out_dir+this.project_name+".1.fastq";
		BufferedWriter bw_1 = new BufferedWriter(new FileWriter(output_fastq_1, false));
		for (int i=0;i< this.num_haps;i++) {
			String fastq_prefix= this.out_dir+"tmp_"+Integer.toString(i);
			String fil = fastq_prefix+".bwa.read1.fastq";
			BufferedReader br_1 = new BufferedReader(new FileReader(fil));
			String line;
			while ((line = br_1.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "");
				bw_1.write(line+"\n");
			}
			br_1.close();
		}
		bw_1.close();
		
		String output_fastq_2 = this.out_dir+this.project_name+".2.fastq";
		BufferedWriter bw_2 = new BufferedWriter(new FileWriter(output_fastq_2, false));
		for (int i=0;i< this.num_haps;i++) {
			String fastq_prefix= this.out_dir+"tmp_"+Integer.toString(i);
			String fil = fastq_prefix+".bwa.read2.fastq";
			BufferedReader br_2 = new BufferedReader(new FileReader(fil));
			String line;
			while ((line = br_2.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "");
				bw_2.write(line+"\n");
			}
			br_2.close();
		}
		bw_2.close();		
	}
	
}
