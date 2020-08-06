package virushaplo;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;


public class VariantsCall {
	public String gatk;
	public String bwa;
	public String samtools;
	public String project_name;
	public String fq_1;
	public String fq_2;
	public String reference;
	public String out_dir;
	public String java;
	public int min_read_cov= 20;
	
	public VariantsCall(String parameter_file) throws IOException {
		InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
//        this.gatk = prop.getProperty("GATK") ;
     
        this.bwa = prop.getProperty("bwa") ;
        this.samtools = prop.getProperty("samtools") ;
        this.fq_1 = prop.getProperty("Fastq_1_Path") ;
        this.fq_2 = prop.getProperty("Fastq_2_Path") ;
        this.out_dir= prop.getProperty("Output_Path")+"/" ;
        this.reference= prop.getProperty("Reference_Seq") ;
        this.project_name = prop.getProperty("Proj_Name") ;
        is.close();
	}
	
	public void bwa()  throws IOException, InterruptedException {
		new File(String.valueOf(this.out_dir) ).mkdir();
		//$bwa mem $ref $prefix\.bwa.read1.fastq $prefix\.bwa.read2.fastq 
//		Runtime r=Runtime.getRuntime();
//		String cmd = this.bwa +" mem "+this.reference+" "+this.fq_1+" "+ this.fq_2+" > "+this.out_dir+"tmp.sam";
//		System.out.println(cmd);
//		Process p=r.exec(cmd);
//		InputStream is=p.getInputStream();
//		InputStreamReader ir=new InputStreamReader(is);
//		BufferedReader br=new BufferedReader(ir);
//		String str=null;
//		while((str=br.readLine())!=null){
//			System.out.println(str);
//		}
//	    int ret=p.waitFor();
//	    int exit_v=p.exitValue();
//	    System.out.println("return value:"+ret);
//	    System.out.println("exit value:"+exit_v);
		new File(String.valueOf(this.out_dir) ).mkdir();
		new File(String.valueOf(this.out_dir)+"/intermediate" ).mkdir();
		BufferedWriter bw = new BufferedWriter(new FileWriter(this.out_dir+"/intermediate/"+ 
				this.project_name+".sam", false));
		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.bwa,
	                "mem", 
	                this.reference, 
	                this.fq_1,
	                this.fq_2
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
	            String line;
	            while ((line = br.readLine()) != null) {
	            	line= line.replace("\n", "").replace("\r", "");
	                bw.write(line+"\n");
	            }
	            CMDProcess.waitFor();
	            System.out.println("Finished aligning reads.");
	            br.close();
		bw.close();
	}
	
	public void sam2bam()  throws IOException, InterruptedException {

		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.samtools,
	                "view", 
	                "-Shub", 
	                this.out_dir+"/intermediate/"+this.project_name+".sam",
	                "-o",
	                this.out_dir+"/intermediate/"+this.project_name+".raw.bam"
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            
	            CMDProcess.waitFor();
	            System.out.println("Finished Converting sam to bam.");
	}
	
	public void sort()  throws IOException, InterruptedException {

		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.samtools,
	                "sort", 
	                "-o", 
	                this.out_dir+"/intermediate/"+this.project_name+".bam",
	                this.out_dir+"/intermediate/"+this.project_name+".raw.bam"
	                
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            
	            CMDProcess.waitFor();
	            System.out.println("Finished sorting the bam file.");
	            
	   
	}
	
	public void index()  throws IOException, InterruptedException {

		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.samtools,
	                "index", 
	                this.out_dir+"/intermediate/"+this.project_name+".bam"
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            CMDProcess.waitFor();
	            System.out.println("Finished index the bam file.");
	}
	
	
	
	public void AddOrReplaceReadGroups()  throws IOException, InterruptedException {

	    ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "AddOrReplaceReadGroups",
	                "-I",
	                this.out_dir+this.project_name+".srt.bam",
	                "-O",
	                this.out_dir+this.project_name+".rg.bam",
	                "-R",
	                this.reference,
	                "-ID",
	                this.project_name,
	                //-LB NPD -PL Illumina -PU NPD -SM 
	                "-LB",
	                "NPD",
	                "-PL",
	                "Illumina",
	                "-PU",
	                "NPD",
	                "-SM",
	                this.project_name
	                
	                
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
	            System.out.println("Finished AddOrReplaceReadGroups the bam file.");
	}	
	
	
//	public void call()  throws IOException, InterruptedException {
//
//	    ProcessBuilder CMDLine = new ProcessBuilder(this.java,
//	                "-jar", 
//	                this.gatk,
//	                "HaplotypeCaller",
//	                "-I",
//	                this.out_dir+this.project_name+".rg.bam",
//	                "-R",
//	                this.reference,
//	             
//	                "-O",
//	                this.out_dir+this.project_name+".tmp.vcf"
//
//	               );
//	            Process CMDProcess = CMDLine.start();
//	            CMDProcess.waitFor();
//	            System.out.println("Finished calling the variants using GATK.");
//	}	
	
	public void call()  throws IOException, InterruptedException {
//		$java -jar $gatk  HaplotypeCaller -R  $ref -I $inbam -ERC  GVCF -ploidy 8 --heterozygosity 0.01  --max-alternate-alleles 1 -O $outgvcf
	    ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "HaplotypeCaller",
	                "-I",
	                this.out_dir+this.project_name+".rg.bam",
	                "-R",
	                this.reference,
	                "-ERC",
	                "GVCF",
	                "-ploidy",
	                "8",
	                "--heterozygosity",
	                "0.01",
	                "-max-alternate-alleles",
	                "1",
	                "-O",
	                this.out_dir+this.project_name+".raw.g.vcf"
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
	            System.out.println("Finished calling the variants using GATK.");
	}	
	
	
	public void GenotypeGVCFs()  throws IOException, InterruptedException {
//		$java -jar $gatk   GenotypeGVCFs -R $ref -V  $prefix\.g.vcf -ploidy 8 -O $prefix\.raw.vcf	    
		ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "GenotypeGVCFs",
	                "-R",
	                this.reference,
	                "-ploidy",
	                "8",
	                "-V",
	                this.out_dir+this.project_name+".raw.g.vcf",
	                "-O",
	                this.out_dir+this.project_name+".raw.vcf"
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
//	            System.out.println("Finished calling the variants using GATK.");
	}	
	
	public void SelectVariants()  throws IOException, InterruptedException {
//		$java -jar $gatk  SelectVariants -R $ref -V  $prefix\.raw.vcf -O $prefix_vcf\.vcf    
		ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "SelectVariants",
	                "-R",
	                this.reference,

	                "-V",
	                this.out_dir+this.project_name+".raw.vcf",
	                "-O",
	                this.out_dir+this.project_name+".tmp.vcf"
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
//	            System.out.println("Finished calling the variants using GATK.");
	}	
	
	
	public double  freq_get(String ss)  throws IOException, InterruptedException {
		
		String[] tmp = ss.split("\t");  
		String tt = tmp[9];
		String[] tmp_tt = tt.split(":");
		String yy= tmp_tt[1];
		double ref =Double.parseDouble(yy.split(",")[0]);
		double alt =Double.parseDouble(yy.split(",")[1]);
		if ((ref+alt)<min_read_cov) {
			return 0.0;
		} else {
			return alt/(ref+alt);
		}
			
		
	}
	public void filter_low()  throws IOException, InterruptedException {
		int snp_end= -1;
		BufferedWriter bw = new BufferedWriter
				(new FileWriter(this.out_dir+this.project_name+".vcf", false));
		BufferedReader br = new BufferedReader(new FileReader(this.out_dir+this.project_name+".tmp.vcf" ));
		String line;
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			if (line.substring(0, 1).equals("#")) {
				bw.write(line+"\n");
			}else {
				double freq = freq_get(line);
				String []tmp =line.split("\t"); 
				if ((freq>0.01) &&  (freq<0.99) && (tmp[3].length()==1) && 
						(tmp[4].length()==1)) {
					int this_start = Integer.parseInt(tmp[1]);
					int this_end = Integer.parseInt(tmp[1]) + tmp[3].length()-1;
					if (this_start> snp_end) {
						bw.write(line+"\n");
						snp_end= this_end;
					}
				}
			}
		}
		br.close();
		bw.close();
	}	
	
	
	
	
}
