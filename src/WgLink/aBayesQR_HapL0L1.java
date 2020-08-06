package virushaplo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import java.io.File;

public class aBayesQR_HapL0L1 {
    HapConfig potential_haps; // all potential haplotypes out of GC (or other algorithms)
    public double lambda; // the parameter that sets the weight of regularizer; not used in v8
    double raw_hap_freq_cutoff;
    String memory_usage; // in bits
    double sum1_weight; // the weight for the first sample
    double[] maf_weights;
    double[][] ld_weights;
    int num_loci;
    int num_haps;
    double total_ld_covered_reads_counts;
    double total_ld_covered_paired_alleles;
    int[][] total_counts;
    int [] single_loc_counts;
    double lasso_coverage_weight;
    double lasso_distance_max_weight;
    String sequencing_technology;
    int l0l1_regression_level;
    
    HashMap<Integer, Integer> site_locations2index;
    HashMap<Integer, Integer> site_index2locations;
    HashMap <Integer,Boolean> site_indexisnull;
    double[][] l0l1_local_matrix;
	double [] l0l1_single_locationfreq ;
    double[] maf; // will be calculated differently for in-pool MAF and population MAF
    double[][] ld_matrix; // will be calculated differently for in-pool LD and population LD
    double [][] allele_matrix;
    String prefix;
    HashMap<Integer, Integer> index_2_pos_dict;
    double[] out_hap_freqs;
    String [] rscript_files;
    HashMap<Integer,String>  pos_come_from_region;
    int num_regions_from;
    aBayesQR_DivideConquer dc_param;
    double min_cov = 20;
    
	public aBayesQR_HapL0L1( aBayesQR_DivideConquer dc, String vef) throws IOException, InterruptedException {
		double total_haps = 0;
		for (int i=0;i< dc.region_start.length; i++) {
			HapConfig hc_tmp = new HapConfig (dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
					+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_inter.freq");
			total_haps= total_haps+ (double) hc_tmp.Haps.length;
		}
		int average_haps = (int)Math.round(total_haps/ (double) dc.region_start.length);		
		
		int num_merge_regions = dc.region_start.length;
		
		int start_index =0;
		int start_index_ext=0;
		int end_index =1;
				
		while (end_index <   num_merge_regions   ) {
			
			HapConfig hc1 = new HapConfig (dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[start_index_ext])+"/"+ dc.project_name+"_inter.freq");
			
			HapConfig hc2 = new HapConfig (dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[end_index])+"_"
					+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+"_inter.freq");
			
			int tiling_index =-1;
			for (int i=0;i< dc.tiling_region_start.length;i++) {
				if ( (dc.tiling_region_start[i] <=  dc.region_end[start_index_ext] )  && 
					(dc.tiling_region_end[i] >=  dc.region_start[end_index] ) ) {
					tiling_index= i;
				}
			}
			
			HapConfig hc_paste = new HapConfig (dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[tiling_index])+"_"
					+  Integer.toString(dc.tiling_region_end[tiling_index])+"/"+ dc.project_name+"_inter.freq");
			
			HapConfig hc12 =   hc2.merge2regions(hc1, hc2, hc_paste, dc.mismatch_tolerance, dc.min_hap_freq,average_haps);
			
			while ( (hc12.Haps.length< dc.Max_L0L1_Regional_Haps)  &&  (end_index<(dc.region_start.length-1) ) )   {
				end_index++;
				hc2 = new HapConfig (dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[end_index])+"_"
						+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+"_inter.freq");
				
				tiling_index =-1;
				for (int i=0;i< dc.tiling_region_start.length;i++) {
					if ( (dc.tiling_region_start[i] <=  dc.region_end[end_index-1] )  && 
						(dc.tiling_region_end[i] >=  dc.region_start[end_index] )  ) {
						tiling_index= i;
					}
				}
				hc_paste = new HapConfig (dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[tiling_index])+"_"
						+  Integer.toString(dc.tiling_region_end[tiling_index])+"/"+ dc.project_name+"_inter.freq");
				
				
				
				hc12 =   hc2.merge2regions(hc12, hc2, hc_paste, dc.mismatch_tolerance,dc.min_hap_freq, average_haps);
				
				new File(dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
						+ Integer.toString(dc.region_end[end_index])).mkdir();
				
				hc12.write2file(hc12, dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
						+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".potential.haps");
				
			}
			
			new File(dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index])).mkdir();
			
			System.out.println("Extracting the sum of alternative allele frequencies at individual sites "
					+ "and the sum of two-variant alternate allele frequencies at paired sites on region:\t"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index]) );
			
			hc12.write2file(hc12, dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".potential.haps");
			
			
			this.calcualte_LD_matrix(vef, dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".potential.haps" , dc);
			
			this.scriptR( dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".R", dc.out_dir+"/intermediate/"+ 
							Integer.toString(dc.region_start[start_index])+"_"
							+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".potential.haps.regression.in", 
							dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
							+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".potential.haps.regression.out"  );
			
			System.out.println("Running L0L1 penalty on region:\t"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index]) );
			
			this.run_R(dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".R");
			
			
			this.generate_inter_freq(dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
					+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".potential.haps", 
					dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
							+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+".potential.haps.regression.out" ,
							dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[start_index])+"_"
							+ Integer.toString(dc.region_end[end_index])+"/"+ dc.project_name+"_inter.freq" );
			
			start_index_ext = end_index;
			end_index++;
			
		}
	}
	
	public void generate_inter_freq(String potential_hap_fil, String r_result, String out_fil) 
			throws IOException, InterruptedException {
		HashMap<Integer, Double>  index2freq = new HashMap<Integer, Double>  ();
		BufferedReader br_r = new BufferedReader(new FileReader(r_result));
		String line = br_r.readLine();
		while ((line = br_r.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				String[] tmp = line.split("\t");
				int h = Integer.parseInt( tmp[0].substring(1) );
				double freq = Double.parseDouble(tmp[1]);
				index2freq.put(h+1, freq);
		}
		br_r.close();
		
		String ss="";
		BufferedReader br_ph = new BufferedReader(new FileReader(potential_hap_fil));
		BufferedWriter bw = new BufferedWriter(new FileWriter( out_fil, false));
		while ((line = br_ph.readLine()) != null) {
			if (line.startsWith("H")) {
				String[] tmp = line.split("\t");
				ss =tmp[0];
				int index =1;
				for (int i=1;i< tmp.length;i++)  {
					if (index2freq.containsKey(i) ) {
					ss=ss+"\th"+Integer.toString(index);
					index++;
					}
				}
				bw.write(ss+"\n");
			} else if (line.startsWith("F")) {
				String[] tmp = line.split("\t");
				ss =tmp[0];
				for (int i=1;i< tmp.length;i++)  {
					if (index2freq.containsKey(i) ) {
						ss=ss+"\t"+Double.toString(index2freq.get(i));
					}
				}
				bw.write(ss+"\n");
				
			} else if (line.startsWith("0")) {
				String[] tmp = line.split("\t");
				ss =tmp[0];
				for (int i=1;i< tmp.length;i++)  {
					if (index2freq.containsKey(i) ) {
					ss=ss+"\t"+tmp[i];
					}
				}
				bw.write(ss+"\n");
			}
		}
		br_ph.close();
		bw.close();
	}
	
	
	public void run_R(String r_fil) throws IOException, InterruptedException {
		ProcessBuilder CMDLine = new ProcessBuilder(this.dc_param.R,
				r_fil
               );
            Process CMDProcess = CMDLine.start();
            
            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
            String line;
            while ((line = br.readLine()) != null) {
            	line= line.replace("\n", "").replace("\r", "");
                System.out.println(line);
            }
            CMDProcess.waitFor();
            System.out.println("Finished R.");
            br.close();
	}
	
	public void scriptR(String r_fil, String infile, String outfile) throws IOException  {
		
		String rfile = r_fil;
		int minSize=1;
		BufferedWriter bw = new BufferedWriter(new FileWriter( rfile, false));
		bw.write("library(\'L0Learn\')\n");
		
		bw.write("maxSize= "+Integer.toString(this.dc_param.maxSuppSize)+"\n");
		bw.write("minSize= "+Integer.toString(minSize)+"\n");
		bw.write("mum_sites="+Integer.toString(this.num_loci)+"\n");
		bw.write("maf_weights= "+Double.toString(this.dc_param.maf_weights)+"\n");
		String tmp ="txt<-as.matrix(read.table(\'";
		
		tmp=tmp+infile+"\',sep=\' \'))";
		bw.write(tmp+"\n");
		bw.write("y<-as.numeric(txt[,2])\n");//apply(X,2,as.numeric)
		bw.write("X<-apply(txt[,3:ncol(txt)], 2,as.numeric)  \n");
		bw.write("if ((ncol(X))>minSize){\n");
//	    			cvfit = L0Learn.cvfit(X, y, nFolds=5, seed=1, penalty="L0L1", nGamma=5, gammaMin=0.001, gammaMax=0.1, maxSuppSize=maxSize)
		tmp= "cvfit = L0Learn.cvfit(X, y, nFolds=5, seed=1, penalty=\"L0L1\", "
				+ "nGamma="+Integer.toString(this.dc_param.nGamma)+", gammaMin="+ Double.toString(this.dc_param.regression_gamma_min)+","
				+ " gammaMax="+ Double.toString(this.dc_param.regression_gamma_max)+ ", maxSuppSize=maxSize"+", intercept= FALSE)";

		bw.write(tmp+"\n");
//		bw.write("optimalGammaIndex= which(unlist(lapply(cvfit$cvMeans, min)) == "
//				+ "min(unlist(lapply(cvfit$cvMeans, min))))\n");
//		bw.write("optimalLambdaIndex = which.min(cvfit$cvMeans[[optimalGammaIndex]])\n");
//		bw.write("optimalLambda  = cvfit$fit$gamma[optimalGammaIndex]\n");
		
		bw.write("optimalLambda  = "+ Double.toString(this.dc_param.regression_lambda)+ "\n");
		bw.write("optimalGammaIndex= which(unlist(lapply(cvfit$cvMeans, min)) == "
		+ "min(unlist(lapply(cvfit$cvMeans, min))))\n");
		bw.write("tmp =coef(cvfit, lambda=optimalLambda, gamma=cvfit$fit$gamma[optimalGammaIndex])\n");	
		
				tmp= "write( paste( \"#Hap_ID\", \"Freq\" ,\"Haplotype\",sep = \"\\t\" ) "
				+ ",file=\""+outfile+"\",append=FALSE)";
		bw.write(tmp+"\n");
		bw.write("for(i in 1:length(as.vector(tmp@x))){\n");
		
		bw.write("hap = gsub(\", \",\"\",toString(round(X[, "
				+ "as.vector(tmp@i)[i]+1][2:(1+ mum_sites)]/maf_weights)))\n");
		tmp= "write(paste( paste( \"h\", toString(as.vector(tmp@i)[i]),sep = \"\" ), "
				+ "as.vector(tmp@x)[i],hap, sep = \"\\t\"),file=\""+outfile
				+"\",append=TRUE)";
		bw.write(tmp+"\n");
		bw.write("}}\n");
		bw.write("print (\"FINISH\")\n");
		bw.close();
		
		
	}
	public boolean is_cross_window(int pos1, int pos2) {
		int window_index_1 = -1;
		int window_index_2 =-1;
		for (int i=0;i< dc_param.region_start.length; i++) {
			if ((pos1>= dc_param.region_start[i]) && (pos1<= dc_param.region_end[i]) ) {
				window_index_1 =i;
			}
			if ((pos2>= dc_param.region_start[i]) && (pos2<= dc_param.region_end[i]) ) {
				window_index_2 =i;
			}
		}
		if (window_index_1== window_index_2) {
			return false;
		} else {
			return true;
		}
	}
	
	public  void write_matrix (String fil ) {
		try {
			double coeff= 0.0;
            BufferedWriter bw = new BufferedWriter(new FileWriter(fil));
            bw.write(Integer.toString(0)+"| "+this.dc_param.sum_weights + "");

            // Originally, num_global_hap.
            for (int h = 0; h < this.allele_matrix.length; h++) {
//                bw.write(" " + (h + 1) + ":" + this.sum1_weight);
                bw.write(" "  + this.dc_param.sum_weights);
            }
            bw.write("\n");
//
            for (int i=0;i< this.num_loci;i++ ) {
            	if  ( this.site_indexisnull.get(i) ) {
            		int loc = this.site_index2locations.get(i);
            		bw.write( Integer.toString(loc) +"|"+
                    		" "+this.dc_param.maf_weights  * this.l0l1_single_locationfreq[i]   );
            		
                    for(int h = 0; h < this.allele_matrix.length; h++) {
                    	coeff= 0.0;
                    	if ( (this.allele_matrix[h][i]) > 0.0) {
                    		coeff = 1.0;
                    	}
                        bw.write(" "
                            + this.dc_param.maf_weights* coeff);
                    }
                    bw.write("\n");
            	} 
//            	else {
//            		int loc = this.site_index2locations.get(i);
//            		bw.write( Integer.toString(loc) +"|"+
//                    		" "+0.0 * this.l0l1_single_locationfreq[i]   );
//            		
//                    for(int h = 0; h < this.allele_matrix.length; h++) {
//                        bw.write(" "
//                            + 0.0
//                            * this.allele_matrix[h][i]);
//                    }
//                    bw.write("\n");
//            	}
            }
 //			
            double coeff1=0.0;
            double coeff2=0.0;
//            for (int i=0;i< this.num_loci;i++ ) {
//            	System.out.println(i+"\t" +this.site_index2locations.get(i) );
//            }
            for (int i=0;i< this.num_loci;i++ ) {
            	for (int j=(i+1); j< this.num_loci;j++) {
//            		System.out.println( "#\t" + this.total_counts[i][j]);
//            		System.out.println( "*\t" + this.l0l1_local_matrix[i][j]);
            		
            		double cross_window_coeff= dc_param.ld_weights*0.5;
            		if (is_cross_window(this.site_index2locations.get(i), this.site_index2locations.get(j)) ) {
            			cross_window_coeff = dc_param.ld_weights* 2.0;
            		}
            		int insert_len = Math.abs(this.site_index2locations.get(i)-this.site_index2locations.get(j) ); 
            		double allzero= 0.0;
            		for(int h = 0; h < this.allele_matrix.length; h++) {
            			if  ( (this.allele_matrix[h][i]>0.0 ) &&  (this.allele_matrix[h][j]>0.0 ) ) {
            				allzero=allzero+1.0;
            			}
            		}
            		
            		if( (cross_window_coeff> 0.0) && (allzero>0)&& (insert_len<
            				Integer.parseInt(dc_param.max_insert_length))  )  {
	            		if  ( ( this.site_indexisnull.get(i) )   && ( this.site_indexisnull.get(j)) 
	            				&& ((!Double.isNaN(this.l0l1_local_matrix[i][j]))) &&(this.total_counts[i][j]> this.min_cov) ) {
	            			bw.write(Integer.toString(this.site_index2locations.get(i))+":"+ 
	                        		Integer.toString(this.site_index2locations.get(j))+"|" 
	                        		+" "+this.dc_param.ld_weights  * this.l0l1_local_matrix[i][j]*cross_window_coeff+"");
	            			for(int h = 0; h < this.allele_matrix.length; h++) {
	            				coeff1=coeff2=0.0;
	            				if  (this.allele_matrix[h][i]>0.0 ) {
	            					coeff1 = 1.0;
	            				}
	            				if  (this.allele_matrix[h][j]>0.0 ) {
	            					coeff2 = 1.0;
	            				}
	            				bw.write(" "+ this.dc_param.ld_weights * cross_window_coeff*coeff1* coeff2);
	            			}
	            			bw.write("\n");
	            		}
            		}
            		
            	}
            }
            bw.close();
		}  catch(Exception e) {
            e.printStackTrace();
        }
		
	}
	
	
	public void calcualte_LD_matrix(String vef_file, String hap_fil, aBayesQR_DivideConquer dc) 
			throws IOException, InterruptedException {
		this.dc_param = dc;
		this.site_locations2index = new HashMap<Integer, Integer>  ();
		this.site_indexisnull = new HashMap<Integer, Boolean>  ();
		this.site_index2locations= new HashMap<Integer, Integer>  ();
		String line ="";
		BufferedReader br_hap = new BufferedReader(new FileReader(hap_fil));
		line = br_hap.readLine();
		line= line.replace("\n", "").replace("\r", "").replace(" ", "");
		int index = 0;
		while ((line = br_hap.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
			this.num_haps = tmp.length-1;
			if (line.startsWith("0")) {
				String info = tmp[0];
				int pos = Integer.parseInt(info.split(";")[1]);
				this.site_index2locations.put(index, pos); 
				this.site_locations2index.put(pos, index); 
				boolean ok=true;
				for (int i=1;i< tmp.length;i++) {
					String snv = tmp[i];
					if (snv.equals("9")) {
						ok=false;
					}
				}
				this.site_indexisnull.put(index, ok); 
				index++;
			} 
		}
		HapConfig hc = new HapConfig (hap_fil);
		this.allele_matrix = new double [hc.Haps.length] [index];
		for (int i=0 ;i< this.allele_matrix.length;i++) {
			for (int j=0; j<  this.allele_matrix[i].length;j++) {
				this.allele_matrix[i][j] = Double.parseDouble( hc.Haps[i].substring(j, j+1) );
			}
		}
		this.num_loci = index;
		br_hap.close();
		
		double[][] local_matrix = new double[this.num_loci][this.num_loci];
		double [] single_locationfreq = new double [this.num_loci];
        HashSet<String> pos_pairs_set = new HashSet<String>();
        this.total_ld_covered_paired_alleles=0;
        this.total_counts = new int[this.num_loci][this.num_loci];
        this.single_loc_counts  = new int [this.num_loci];
				
        try {
            BufferedReader br = new BufferedReader(new FileReader(vef_file));
            line = br.readLine();
            while (line!=null) {
                String[] location_alleles = line.split("\t")[1].split(";");
                if (location_alleles.length < 2) { // no multiple sites in this read
                    line = br.readLine();
                    continue;
                }

                // The search may fail, so we don't know the length of the array.
                ArrayList<Integer> indexes = new ArrayList<Integer>();
                ArrayList<String> alleles = new ArrayList<String>();
                ArrayList<Integer> positions = new ArrayList<Integer>();
                for (int i = 0; i < location_alleles.length; i++) {
                    String[] tmp_loc_allele = location_alleles[i].split("=");
                    int the_location = Integer.parseInt(tmp_loc_allele[0]);

                    // TODO: [LEFTOVER]
                    // setup_site_locations2index_map<POS_ON_REF,POS_IND_HAPCONFIG>

                    // If we do find the site in question...
                    if (this.site_locations2index.containsKey(the_location)) {
                        indexes.add(this.site_locations2index.get(the_location));
                        alleles.add(tmp_loc_allele[1]);
                        positions.add(the_location);
                    }
                }
                int num_avail_locs=indexes.size();
                for (int site1 = 0; site1 < num_avail_locs; site1++) {
                	int single_index = indexes.get(site1);
                	this.single_loc_counts [single_index]++;
                	if (alleles.get(site1).equals("1") ) {
                		single_locationfreq[single_index]++;
                	}
                    for (int site2 = site1 + 1; site2 < num_avail_locs; site2++) {
                        int index1 = indexes.get(site1);
                        int index2 = indexes.get(site2);

	                        if (index1 > index2) { // ld_matrix only record loc1 < loc2 entries
	                            int tmp = index1;
	                            index1 = index2;
	                            index2 = tmp;
	                        }
	                        if  ( (index2- index1) > 200 ) {
//	                        	System.out.println(line);
	                        }
	                        
	                        this.total_counts[index1][index2]++;
	                        this.total_ld_covered_reads_counts =
	                        		this.total_ld_covered_reads_counts+1;
	                        String tmp = Integer.toString(index1)+":"+ Integer.toString(index2);
	                        if (!pos_pairs_set.contains(tmp)) {
	                        	pos_pairs_set.add(tmp);
	                        	this.total_ld_covered_paired_alleles= 
	                        			this.total_ld_covered_paired_alleles+1;
	                        }
	                        // Both are "1".
	                        if (alleles.get(site1).equals("1") && alleles.get(site2).equals("1")) {
	                            local_matrix[index1][index2]++;  
	                        }
                    }
                }

                line=br.readLine();
            }

            br.close();

            // Divide by the total number of reads covering both sites.
            for (int loc1 = 0; loc1 < this.num_loci; loc1++) {
            	if (this.single_loc_counts[loc1]!=0) {
            		single_locationfreq[loc1] = single_locationfreq[loc1]/ this.single_loc_counts[loc1];
            	} else {
            		single_locationfreq[loc1]= Double.NaN;
            	}
                for (int loc2 = 0; loc2 < this.num_loci; loc2++) {
                    if (this.total_counts[loc1][loc2] != 0) {
                        local_matrix[loc1][loc2] = local_matrix[loc1][loc2]
                            / this.total_counts[loc1][loc2];
                        
                    } else {
                        local_matrix[loc1][loc2] = Double.NaN;
                    }
                }
            }
        } catch(Exception e) {
            e.printStackTrace();
        }
        
//        double[][] l0l1_local_matrix;
//    	double [] l0l1_single_locationfreq ;
        this.l0l1_local_matrix= local_matrix.clone();
        this.l0l1_single_locationfreq = single_locationfreq.clone();
		
        write_matrix (hap_fil+".regression.in" );
		
	}
}
