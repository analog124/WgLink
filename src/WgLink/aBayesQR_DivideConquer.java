package virushaplo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Properties;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class aBayesQR_DivideConquer {
	public int reconstruction_start;
	public int reconstruction_end;
	public int region_length;
	String out_dir;
	public int [] region_start;
	public int [] region_end;
	public int [] tiling_region_start;
	public int [] tiling_region_end;
	String project_name;
	public String reference;
	public String SNV_cutoff;
	public String min_mapping_qual;
	public String min_read_length;
	public String max_insert_length;
	public String sequence_err;
	public String MEC_improvement_cutoff;
	public int num_threads;
	public String aBayesQR;
	public String [] cmds;
	public int thread_index;
	public int Max_L0L1_Regional_Haps;
	public double min_hap_freq;
	public double mismatch_tolerance;
	public String tool; 
	public double sum_weights;
	public double maf_weights;
	public double ld_weights;
	public String R;
	public int nGamma;
	int maxSuppSize;
	public double regression_gamma_min;
    public double regression_gamma_max;
    public double regression_lambda;
	
	
	final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
	
	public aBayesQR_DivideConquer(String parameter_file) throws IOException, InterruptedException {
		
		InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
        this.out_dir= prop.getProperty("Output_Path")+"/" ;
        this.reconstruction_start = Integer.parseInt(prop.getProperty("Reconstruction_Start"));
        this.reconstruction_end = Integer.parseInt(prop.getProperty("Reconstruction_End"));
        this.region_length = Integer.parseInt(prop.getProperty("Region_Length"));
        this.project_name= prop.getProperty("Proj_Name");
        this.reference = prop.getProperty("Reference_Seq") ;
        this.SNV_cutoff = prop.getProperty("SNV_Cutoff") ;
        this.min_mapping_qual = prop.getProperty("Min_Mapping_Qual") ;
        this.min_read_length = prop.getProperty("Min_Read_Length") ;
        this.max_insert_length = prop.getProperty("Max_Insert_Length") ;
        this.sequence_err = prop.getProperty("Sequence_Err") ;
        this.MEC_improvement_cutoff= prop.getProperty("MEC_Improvement_Cutoff") ;
        this.num_threads = Integer.parseInt(prop.getProperty("Number_Threads"));
        this.aBayesQR = prop.getProperty("aBayesQR") ;
        this.min_hap_freq =  Double.parseDouble(prop.getProperty("Min_Hap_Freq"));
        this.mismatch_tolerance =  Double.parseDouble(prop.getProperty("BFS_Mismatch_Tolerance_Rate"));
        this.Max_L0L1_Regional_Haps = Integer.parseInt(prop.getProperty("Max_L0L1_Regional_Haps"));
        
        this.sum_weights = Double.parseDouble(prop.getProperty("Regression_One_Vector_Weight"));
        this.maf_weights = Double.parseDouble(prop.getProperty("Regression_Hap_MAF_Weight"));
        this.ld_weights = Double.parseDouble(prop.getProperty("Regression_Hap_LD_Weight"));
        this.R= prop.getProperty("Rscript_path") ;
        this.maxSuppSize= Integer.parseInt(prop.getProperty("Maximum_Haps_R"));
        this.nGamma= Integer.parseInt(
                prop.getProperty("Regression_n_Gamma"));
        
        this.regression_gamma_max=Double.parseDouble
        		(prop.getProperty("Regression_Gamma_Max"));
        this.regression_gamma_min=Double.parseDouble
        		(prop.getProperty("Regression_Gamma_Min"));
        
        this.regression_lambda = Double.parseDouble
        		(prop.getProperty("Regression_Lambda"));
        
//        
//        Regression_One_Vector_Weight = 5.0
//        		Regression_Hap_MAF_Weight = 2.0
//        		Regression_Hap_LD_Weight = 1.0

        
        
        is.close();
        int num_regions = (this.reconstruction_end-this.reconstruction_start+1) / this.region_length;
        if (((this.reconstruction_end-this.reconstruction_start+1) % this.region_length ) !=0) {
        	num_regions=num_regions+1;
        }
        
        this.region_start= new int [num_regions];
        this.region_end= new int [num_regions];
        this.tiling_region_start = new int [num_regions-1 ];
        this.tiling_region_end = new int [num_regions-1 ];
        
        for (int i=0; i< num_regions;i++) {
        	this.region_start[i]= this.reconstruction_start+ i* this.region_length;
        	this.region_end[i]= this.reconstruction_start+ (i+1)* this.region_length-1;
        	if (this.region_end[i]> this.reconstruction_end) {
        		this.region_end[i] = this.reconstruction_end;
        	}
        }
        
        for (int i=0; i< (num_regions-1);i++) {
        	this.tiling_region_start[i]=  (this.region_start[i] + this.region_end[i])/2+1;
        	this.tiling_region_end[i]=  (this.region_start[i+1] + this.region_end[i+1])/2;
        }
        
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.out_dir+"/intermediate/DC_PLAN.txt", false));
        String ss ="Splited REGION:";
        bw.write(ss+"\n");
        ss=Integer.toString(this.region_start[0])+":"+Integer.toString(this.region_end[0]) ;
        for (int i=1; i< num_regions;i++) {
        	ss= ss+"\t"+ Integer.toString(this.region_start[i])+":"+Integer.toString(this.region_end[i]) ;
        }
        bw.write(ss+"\n");
        ss ="TILING_REGION:";
        bw.write(ss+"\n");
        ss=Integer.toString(this.tiling_region_start[0])+":"+Integer.toString(this.tiling_region_end[0]) ;
        for (int i=1; i< (num_regions-1);i++) {
        	ss= ss+"\t"+ Integer.toString(this.tiling_region_start[i])+":"+
        			Integer.toString(this.tiling_region_end[i]) ;
        }
        bw.write(ss+"\n");	
        bw.close();
        
        this.Generate_aBayesQR_Parameter_Fil();
        
        this.aBayesQR();
        
        
	}
	
	public void Generate_aBayesQR_Parameter_Fil() throws IOException {
		for (int i=0;i< this.region_start.length;i++) {
			String start = Integer.toString(this.region_start[i]);
			String end = Integer.toString(this.region_end[i]);
			new File(String.valueOf(this.out_dir)+"/intermediate/"+start+"_"+end ).mkdir();
			BufferedWriter bw = new BufferedWriter
					(new FileWriter(String.valueOf(this.out_dir)+"/intermediate/"+start+"_"+end+"/"
							+this.project_name+".aBayesQR.config", false));
			
			bw.write("filename of reference sequence (FASTA) : "+this.reference+"\n");
			bw.write("filname of the aligned reads (sam format) : "+
					this.out_dir+"/intermediate/"+this.project_name+".sam"+"\n");
			bw.write("paired-end (1 = true, 0 = false) : 1\n");
			bw.write("SNV_thres : "+this.SNV_cutoff+"\n");
			bw.write("reconstruction_start : "+Integer.toString( region_start[i]) +"\n");
			bw.write("reconstruction_stop : "+Integer.toString( region_end[i]) +"\n");
			bw.write("min_mapping_qual : "+this.min_mapping_qual+"\n");
			bw.write("min_read_length : "+this.min_read_length+"\n");
			bw.write("max_insert_length : "+this.max_insert_length+"\n");
			bw.write("characteristic zone name : "+this.out_dir+"/intermediate/"+start+"_"+end+"/"
					+this.project_name+"\n");
			bw.write("seq_err (assumed sequencing error rate(%)) : "+this.sequence_err +"\n");
			bw.write("MEC improvement threshold : "+this.MEC_improvement_cutoff+"\n");
			bw.close();
		}
		
		for (int i=0;i< this.tiling_region_start.length;i++) {
			String start = Integer.toString(this.tiling_region_start[i]);
			String end = Integer.toString(this.tiling_region_end[i]);
			new File(String.valueOf(this.out_dir)+"/intermediate/"+start+"_"+end ).mkdir();
			BufferedWriter bw = new BufferedWriter
					(new FileWriter(String.valueOf(this.out_dir)+"/intermediate/"+start+"_"+end+"/"
							+this.project_name+".aBayesQR.config", false));
			
			bw.write("filename of reference sequence (FASTA) : "+this.reference+"\n");
			bw.write("filname of the aligned reads (sam format) : "+
					this.out_dir+"/intermediate/"+this.project_name+".sam"+"\n");
			bw.write("paired-end (1 = true, 0 = false) : 1\n");
			bw.write("SNV_thres : "+this.SNV_cutoff+"\n");
			bw.write("reconstruction_start : "+Integer.toString( tiling_region_start[i]) +"\n");
			bw.write("reconstruction_stop : "+Integer.toString( tiling_region_end[i]) +"\n");
			bw.write("min_mapping_qual : "+this.min_mapping_qual+"\n");
			bw.write("min_read_length : "+this.min_read_length+"\n");
			bw.write("max_insert_length : "+this.max_insert_length+"\n");
			bw.write("characteristic zone name : "+this.out_dir+"/intermediate/"+start+"_"+end+"/"
					+this.project_name+"\n");
			bw.write("seq_err (assumed sequencing error rate(%)) : "+this.sequence_err +"\n");
			bw.write("MEC improvement threshold : "+this.MEC_improvement_cutoff+"\n");
			bw.close();
		}
	}
	
	public void aBayesQR() throws IOException, InterruptedException {
		
		
		
		String [] regions = new String [ this.tiling_region_start.length+ this.region_start.length ];
	
		for (int i=0;i< this.region_start.length;i++) {
			regions[i]= Integer.toString(this.region_start[i])+"_" +Integer.toString(this.region_end[i]);
		}
		for (int i=0;i< this.tiling_region_start.length;i++) {
			regions[this.region_start.length+ i]=
					Integer.toString(this.tiling_region_start[i])+"_" +Integer.toString(this.tiling_region_end[i]);
		}
		
		int num_level = regions.length/ this.num_threads;
		if ((regions.length% this.num_threads )!=0) {
			num_level++;
		}
		String [][] threads_regions = new String [num_level][this.num_threads];
		for (int i=0;i< num_level;i++) {
			for (int j=0;j< this.num_threads;j++) {
				if ((this.num_threads*i+j) < regions.length) {
					threads_regions[i][j]= regions[this.num_threads*i+j];
				} else {
					threads_regions[i][j]="NULL";
				}
			}
		}
		System.out.println("To Run aBayesQR Efficiently, "+Integer.toString(this.num_threads)+" threads are applied:");
		for (int i=0;i< num_level;i++) {
			String ss="";
			for (int j=0;j< this.num_threads;j++) {
				if (!threads_regions[i][j].equals("NULL")) {
					ss= ss+ "\t"+threads_regions[i][j];
				}
			}
			System.out.println(ss);
		}
		
		for (int i=0;i< num_level;i++) {
			int real_num_threads= 0;
			for (int t=0; t< num_threads ;t++ ) {
				if (! threads_regions[i][t].equals("NULL")) {
					real_num_threads++;
				}
			}
			String ss="";
			for (int t=0;t< this.num_threads;t++) {
				if (!threads_regions[i][t].equals("NULL")) {
					ss= ss+ "\t"+threads_regions[i][t];
				}
			}
			System.out.println("------------------------------------");
			System.out.println("Running aBayesQR on Splited Regions:\t"+ss);
			String [] cmdArray= new String [real_num_threads];
			String [] regionArray= new String [real_num_threads];
			for (int t=0;t< real_num_threads;t++) {
				regionArray[t] = threads_regions[i][t];
				cmdArray[t]= this.aBayesQR+" "+ 
						String.valueOf(this.out_dir)+"/intermediate/"+
						threads_regions[i][t] +"/"+this.project_name+".aBayesQR.config";


			}
			aBayesQR ts = new aBayesQR(cmdArray,regionArray );
			ts.run();
			
			int unfinished_regions = 0;
			
			for (int t=0;t< real_num_threads;t++) {
				String outfile=this.out_dir+ "/intermediate/"+ threads_regions[i][t]+"/"+this.project_name
						+"_ViralSeq.txt";
				File file = new File(outfile);
				if (!file.exists() ) {
					System.out.println("Can not Generate aBayesQR Output File for Region:\t"+  threads_regions[i][t] );
					System.out.println("aBayesQR Error!! Please Refer to the Error Information Above." );
					unfinished_regions++;
//					System.exit(0);
				}
			}
			ts= null;
			
//			System.out.println(unfinished_regions );
//			while (unfinished_regions>0) {
//				System.out.println("Retry aBayesQR on the failed regions." );
//				
//				String [] unfinished_cmdArray= new String [unfinished_regions];
//				String [] unfinished_regionArray= new String [unfinished_regions];
//				unfinished_regions =0; 
//				for (int t=0;t< real_num_threads;t++) {
//					String outfile=this.out_dir+ "/intermediate/"+ threads_regions[i][t]+"/"+this.project_name
//							+"_ViralSeq.txt";
//					File file = new File(outfile);
//					if (!file.exists() ) {
//						unfinished_regionArray[unfinished_regions] = threads_regions[i][t];
//						unfinished_cmdArray[unfinished_regions]= this.aBayesQR+" "+ 
//								String.valueOf(this.out_dir)+"/intermediate/"+
//								threads_regions[i][t] +"/"+this.project_name+".aBayesQR.config";
//						unfinished_regions++;
//						
//						ArrayList<String  >  configArray = new ArrayList<String>();
//						BufferedReader br_conf = new BufferedReader(new FileReader(
//								String.valueOf(this.out_dir)+"/intermediate/"+
//								threads_regions[i][t] +"/"+this.project_name+".aBayesQR.config"));
//						String line ="";
//						while ((line = br_conf.readLine()) != null) {
//							line= line.replace("\n", "").replace("\r", "").replace("  ", " ");
//							if ( (line.length()> 24) &&  (line.substring(0, 23).equals("initial population size") ) ) {
//								String[] tmp = line.split(" ");
//								double int_pop =  Double.parseDouble(tmp[tmp.length-1]);
//								String last = Integer.toString( (int) (int_pop*0.9) ) ;
//								String s ="";
//								for (int k=0;k< (tmp.length-1 ); k++) {
//									s=s+ tmp[k]+ " ";
//								}
//								s=s+ last;
//								configArray.add(s); 
//								
//							} else { 
//								configArray.add(line); 
//							}
//						}
//						br_conf.close();
//						
//						BufferedWriter bw_conf = new BufferedWriter
//								(new FileWriter(
//										String.valueOf(this.out_dir)+"/intermediate/"+
//										threads_regions[i][t] +"/"+this.project_name+".aBayesQR.config", false));
//						for (int k=0;k< configArray.size();k++) {
//							bw_conf.write( configArray.get(k)+"\n");
//						}
//
//						bw_conf.close();
//					}
//				}
//
//				
//				aBayesQR ts_unfinished = new aBayesQR(unfinished_cmdArray,unfinished_regionArray );
//				ts_unfinished.run();
//				ts_unfinished= null;
//				
//				unfinished_regions =0; 
//				for (int t=0;t< real_num_threads;t++) {
//					String outfile=this.out_dir+ "/intermediate/"+ threads_regions[i][t]+"/"+this.project_name
//							+"_ViralSeq.txt";
//					File file = new File(outfile);
//					if (!file.exists() ) {
//						unfinished_regions++;
//					}
//				}
//				
//				
//			}
		}
		
		
		
//		System.exit(0);
//		
//		for (int i=0;i< regions.length;i++) {
//			System.out.println("Run aBayesQR on Region:\t"+regions[i]);
//			
//			
//			ProcessBuilder CMDLine = new ProcessBuilder(
//                this.aBayesQR, 
//                String.valueOf(this.out_dir)+"/intermediate/"+
//                		regions[i] +"/"+this.project_name+".aBayesQR.config"
//               );
//			
//            Process CMDProcess = CMDLine.start();
//            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getErrorStream()));
//            String line;
//            while ((line = br.readLine()) != null) {
//            	line= line.replace("\n", "").replace("\r", "");
//            	System.out.println(line);
//            }
//            CMDProcess.waitFor();
//            
//            String outfile=this.out_dir+ "/intermediate/"+ regions[i]+"/"+this.project_name
//					+"_ViralSeq.txt";
//			File file = new File(outfile);
//			if (!file.exists() ) {
//				System.out.println("Can not Generate aBayesQR Output File for Region:\t"+  regions[i]);
//				System.out.println("aBayesQR Error!! Please Refer to the Error Message Above." );
//			}
//			System.out.println("...... ......");
//		}
		
		
	}
	
}
