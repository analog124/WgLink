package virushaplo;

import java.util.ArrayList;
import java.io.InputStreamReader;
import java.time.temporal.TemporalAccessor;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.HashSet;
import java.util.Iterator;

import java.io.Reader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Writer;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import java.io.File;
import java.util.HashMap;


public class Main {
	
	public static void main(final String[] args) throws Exception {
		 System.out.println("Welcome to use WgLink! If you have any questions, please refer to "
		 		+ "https://github.com/theLongLab/WgLink");
		final String[] supported_functions_array = 
			{ "TenSQR",
					"aBayesQR"};
		final HashSet<String> supported_functions = new HashSet<String>();
        for (int k = 0; k < supported_functions_array.length; ++k) {
            supported_functions.add(supported_functions_array[k]);
        }
        final String function = args[0];
        if (!supported_functions.contains(function)) {
            System.out.println("Function " + function + " is not supported. A typo?");
            System.exit(0); 
        }
        if (function.equals("sim")) {
        	String config_path = args[1];
        	final Simulations sim = new Simulations(config_path);
        	sim.generate_fastq();
        } else if (function.equals("TenSQR")) {
        	final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        	String config_path = args[1];
        	System.out.println("Start mapping reads to reference Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
        	final VariantsCall varcall = new VariantsCall(config_path);
        	
        	varcall.bwa();
        	
        	System.out.println("Mapping Finished Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
        	final DivideConquer dc = new DivideConquer ( config_path );
        	final Variants vr= new Variants(config_path, dc);
        	vr.generate_inter_hap(dc);
        	System.out.println("Generating the vef file:\t" +
                	dtf.format(LocalDateTime.now()) + "\n");
        	
        	
        	BAMFormatter.generate_vef( dc.out_dir+"/intermediate/Seq_SNV.txt",dc.project_name, 
        			dc.out_dir+"/intermediate/",  dc.out_dir+"/intermediate/"+ 
        					dc.project_name+".sam");
//        	
        	final HapL0L1 hl0 = new HapL0L1 (dc, dc.out_dir+"/intermediate/"+dc.project_name+".vef" );
        	
        	vr.output_final_fil(dc.out_dir+"/intermediate/Seq_SNV.txt", dc.out_dir+"/intermediate/"+ 
        			Integer.toString(dc.region_start[0])+"_"
					+ Integer.toString(dc.region_end[dc.region_start.length -1])+"/"+ dc.project_name+"_inter.freq", 
					dc.out_dir+"Final.Haps", 
        			dc.reconstruction_start, dc.reconstruction_end);
        	
        	System.out.println("The final haplotype file has been calculated by WgLink:\t" +dc.out_dir+"Final.Haps"+"\t" +
        	dtf.format(LocalDateTime.now()) + "\n");
        	

        	
        }else if (function.equals("aBayesQR")) {
        	String config_path = args[1];
        	final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        	System.out.println("Start mapping reads to reference Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
        	final VariantsCall varcall = new VariantsCall(config_path);
        	varcall.bwa();
        	System.out.println("Mapping Finished Time:\t" + dtf.format(LocalDateTime.now()) + "\n");
        	final aBayesQR_DivideConquer dc = new aBayesQR_DivideConquer ( config_path );
        	final Variants vr= new Variants(config_path, dc);
        	vr.generate_inter_hap(dc);
        	System.out.println("Generating the vef file:\t" +
                	dtf.format(LocalDateTime.now()) + "\n");
        	BAMFormatter.generate_vef( dc.out_dir+"/intermediate/Seq_SNV.txt",dc.project_name, 
        			dc.out_dir+"/intermediate/",  dc.out_dir+"/intermediate/"+ 
        					dc.project_name+".sam");
        	System.out.println("The vef file has been generated:\t" +
                	dtf.format(LocalDateTime.now()) + "\n");
        	final aBayesQR_HapL0L1 hl0 = new aBayesQR_HapL0L1 (dc, dc.out_dir+"/intermediate/"+dc.project_name+".vef" );
        	
        	
        	vr.output_final_fil(dc.out_dir+"/intermediate/Seq_SNV.txt", dc.out_dir+"/intermediate/"+ 
        			Integer.toString(dc.region_start[0])+"_"
					+ Integer.toString(dc.region_end[dc.region_start.length -1])+"/"+ dc.project_name+"_inter.freq", 
					dc.out_dir+"Final.Haps", 
        			dc.reconstruction_start, dc.reconstruction_end);
        	System.out.println("The final haplotype file has been calculated by WgLink:\t" +dc.out_dir+"Final.Haps"+"\t" +
                	dtf.format(LocalDateTime.now()) + "\n");
        	
        }        
	}
}

