package virushaplo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

public class Variants {
	public  HashMap<Integer, Integer> SNV2index;
	public int [] SNV_indexArray;
	public int [] SNV_posArray;
	public int max_haps= 50;
	public int User_Start_Pos;
	public int User_End_Pos;
	public  HashMap<Integer, String> Pos2SNV;
	public  HashMap<Integer, String> Pos2Allele;
	public int [] SNVposArray;
	
	public Variants(String parameter_file, DivideConquer dc) throws IOException, InterruptedException {
		this.SNV2index = new HashMap<Integer, Integer>(); 
		ArrayList<Integer >  posArray = new ArrayList<Integer>();
		HashMap<Integer, String> alleleDic = new HashMap<Integer, String>();
		String[] hapArray = new String [max_haps];
		for (int i=0;i< hapArray.length;i++) {
			hapArray[i]="";
		}
		for (int i=0;i< dc.region_start.length;i++) {
			String til_hap_fil2= dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
					+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
			BufferedReader br_hap = new BufferedReader(new FileReader(til_hap_fil2));
			String line ="";
			int index =0;
			double hap_freq= 0.0;
			while ((line = br_hap.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (line.substring(0,1).equals("V") ) {
					String[] tmp = line.split(":");
					hap_freq =  Double.parseDouble(tmp[1]);
				} else {
					if (hap_freq> dc.min_hap_freq) {
						if (i==0) {
							hapArray[index] =line;
							index++;
						} else {
							if (!hapArray[index].equals("") ) {
								hapArray[index] =hapArray[index] + line;
							} else {
								hapArray[index] = hapArray[index-1];
								hapArray[index] =hapArray[index] + line;
							}
							index++;
						}
					}
				}
				
			}
			for (int j= index; j< this.max_haps;j++) {
				hapArray[j]= hapArray[index-1];
			}
			br_hap.close();
		}
		
//		BufferedWriter bw = new BufferedWriter
//				(new FileWriter(dc.out_dir+"/intermediate/tmp_seq.txt", false));
//		for (int i=0;i< hapArray.length;i++) {
//			if (!hapArray[i].equals("") ) {
//				bw.write(hapArray[i]+"\n");
//			}
//		}
//		bw.close();
		
		String refence="";
		String line = "";
		BufferedReader br_ref = new BufferedReader(new FileReader(dc.reference));
		while ((line = br_ref.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (!line.substring(0,1).equals(">") ) {
					refence=refence+ line;
				} 
		}
		br_ref.close();
		
		User_Start_Pos = dc.region_start[0];
		int count = 0;
		BufferedWriter bw_SNV = new BufferedWriter
				(new FileWriter(dc.out_dir+"/intermediate/Seq_SNV.txt", false));
		for (int p =0;p< hapArray[0].length();p++) {
			
			String allele = "*";
			for (int i=0;i< hapArray.length;i++) {
				if (!hapArray[i].equals("") ) {
					String snp = hapArray[i].substring(p, p+1);
					if (!snp.equals("*")) {
						if ( allele.equals("*" )) {
							allele= snp;
						} else {
							boolean ok= true;
							for (int l=0;l< allele.length();l++) {
								if ( allele.substring(l, l+1).equals(snp) ) {
									ok =false;
								}
							}
							if (ok) {
								allele=allele+snp;
							}
						}
					}
				}
			}
			String out_snp ="";
			String origin = refence.substring(this.User_Start_Pos+ count-1, this.User_Start_Pos+ count);
			boolean  flag1 = false ;
			boolean  flag2 = false ;
			for (int l=0;l< allele.length();l++) {
				if ( !allele.substring(l, l+1).equals(origin) ) {
					out_snp =out_snp+ allele.substring(l, l+1);
					flag1 =true;
				}
				if ( allele.substring(l, l+1).equals(origin) ) {
					flag2 =true;
				}
				
			}
			if ( out_snp.equals("") ) {
				out_snp ="*";
			}
//			if ((flag1)  && (!flag2) ) {
//				out_snp ="*";
//			}
			
//			bw_SNV.write(Integer.toString(this.User_Start_Pos+ count)+"\t"+ ref+"\t"+ allele+"\t" + 
//					refence.substring(this.User_Start_Pos+ count-1, this.User_Start_Pos+ count)+ "\n");
			
			bw_SNV.write(Integer.toString(this.User_Start_Pos+ count)+"\t"+ origin+"\t"+ out_snp+"\n");
			
			count ++ ;
		}
		bw_SNV.close();
	}
	

	
	
	
	
	public Variants(String parameter_file, aBayesQR_DivideConquer dc) throws IOException, InterruptedException {
		this.SNV2index = new HashMap<Integer, Integer>(); 
		ArrayList<Integer >  posArray = new ArrayList<Integer>();
		HashMap<Integer, String> alleleDic = new HashMap<Integer, String>();
		String[] hapArray = new String [max_haps];
		for (int i=0;i< hapArray.length;i++) {
			hapArray[i]="";
		}
		for (int i=0;i< dc.region_start.length;i++) {
			String til_hap_fil2= dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
					+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
			BufferedReader br_hap = new BufferedReader(new FileReader(til_hap_fil2));
			String line ="";
			int index =0;
			double hap_freq= 0.0;
			while ((line = br_hap.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "").replace("N", "*");
				if  (line.substring(0,1).equals("V") ) {
					String[] tmp = line.split(":");
					hap_freq =  Double.parseDouble(tmp[1]);
				} else {
					if (hap_freq> dc.min_hap_freq) {
						if (i==0) {
							hapArray[index] =line;
							index++;
						} else {
							if (!hapArray[index].equals("") ) {
								hapArray[index] =hapArray[index] + line;
							} else {
								hapArray[index] = hapArray[index-1];
								hapArray[index] =hapArray[index] + line;
							}
							index++;
						}
					}
				}
				
			}
			for (int j= index; j< this.max_haps;j++) {
				hapArray[j]= hapArray[index-1];
			}
			br_hap.close();
		}
		
		BufferedWriter bw = new BufferedWriter
				(new FileWriter(dc.out_dir+"/intermediate/tmp_seq.txt", false));
		for (int i=0;i< hapArray.length;i++) {
			if (!hapArray[i].equals("") ) {
				bw.write(hapArray[i]+"\n");
			}
		}
		bw.close();
		
		String refence="";
		String line = "";
		BufferedReader br_ref = new BufferedReader(new FileReader(dc.reference));
		while ((line = br_ref.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (!line.substring(0,1).equals(">") ) {
					refence=refence+ line;
				} 
		}
		br_ref.close();
		
		User_Start_Pos = dc.region_start[0];
		int count = 0;
		BufferedWriter bw_SNV = new BufferedWriter
				(new FileWriter(dc.out_dir+"/intermediate/Seq_SNV.txt", false));
		for (int p =0;p< hapArray[0].length();p++) {
			
			String allele = "*";
			for (int i=0;i< hapArray.length;i++) {
				if (!hapArray[i].equals("") ) {
					String snp = hapArray[i].substring(p, p+1);
					if (!snp.equals("*")) {
						if ( allele.equals("*" )) {
							allele= snp;
						} else {
							boolean ok= true;
							for (int l=0;l< allele.length();l++) {
								if ( allele.substring(l, l+1).equals(snp) ) {
									ok =false;
								}
							}
							if (ok) {
								allele=allele+snp;
							}
						}
					}
				}
			}
			String out_snp ="";
			String origin = refence.substring(this.User_Start_Pos+ count-1, this.User_Start_Pos+ count);
			for (int l=0;l< allele.length();l++) {
				if ( !allele.substring(l, l+1).equals(origin) ) {
					out_snp =out_snp+ allele.substring(l, l+1);
				}
			}
			if ( out_snp.equals("") ) {
				out_snp ="*";
			}
				
			bw_SNV.write(Integer.toString(this.User_Start_Pos+ count)+"\t"+ origin+"\t"+ out_snp+"\n");
			
			count ++ ;
		}
		bw_SNV.close();
	}
	
	public void output_final_fil(String snv_fil, String hap_fil, String out_fil, 
			int start, int end) throws IOException, InterruptedException {
		HapConfig hc = new HapConfig (hap_fil);
		BufferedReader br = new BufferedReader(new FileReader( snv_fil));
		String line="";
		HashMap <Integer, String> pos2ref = new HashMap<Integer, String>  ();
		HashMap <Integer, String> pos2allele = new HashMap<Integer, String>  ();
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace("-", "*");
			String [] tmp =line.split("\t");
			int pos= Integer.parseInt(tmp[0]);
			pos2ref.put(pos, tmp[1]);
			pos2allele.put(pos, tmp[2]);
		}
		
		br.close();
		
		BufferedWriter bw= new BufferedWriter
				(new FileWriter(out_fil, false));
		
		for (int h=0;h< hc.Haps.length;h++) {
			String ss = "Haplotype_"+Integer.toString(h+1)+" Freq:\t"+ Double.toString(hc.Freqs[h]);
			bw.write(ss+"\n");
			ss= "";
			String snv ="";
			for (int p =start; p<=end;p++) {
				String key = Integer.toString(p)+":"+Integer.toString(h);
				if ( hc.Pos2SNV.containsKey(key) ) {
					if  (hc.Pos2SNV.get(key).equals("9") ) {
						snv="N"; 
					} else if (hc.Pos2SNV.get(key).equals("0") ) {
						snv= pos2ref.get(p); 
					} else { 
						snv= pos2allele.get(p).substring(Integer.parseInt( hc.Pos2SNV.get(key))-1 ,
								Integer.parseInt( hc.Pos2SNV.get(key))); 
					}
//					boolean flag = false ;
//					for (int ht=0;ht< hc.Haps.length;ht++) {
//						String key_tmp  = Integer.toString(p)+":"+Integer.toString(ht);
//						if  (hc.Pos2SNV.get(key_tmp).equals("9") ) {
//							flag= true;
//						}
//					}
//					if (flag) {
//						snv= "N";
//					}
				} else { 
					snv= pos2ref.get(p); 
				}
				ss=ss+snv;
			}
			bw.write(ss+"\n");
		}
		
		bw.close();
		
		
		
	}
	

	
	public void rewrite_haplos(String in_fil, String out_fil, int start, int end) throws IOException {
//		System.out.println(in_fil+ out_fil+  start+  end);
		ArrayList<Double >  freqArray = new ArrayList<Double>();
		ArrayList<String >  seqArray = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader( in_fil));
		String line="";
		String hap="";
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace("-", "*");
//			System.out.println(hap);
			if  (line.substring(0,1).equals(">") ) {
				String[] tmp = line.split("_");
				double hap_freq =  Double.parseDouble(tmp[1].replace("*", "-") );
				freqArray.add(hap_freq);
				if (!hap.equals("")) {
					seqArray.add(hap);
					hap="";
				}
			} else {
				hap =hap+ line;
			}
		}
		seqArray.add(hap);
		br.close();
		
		BufferedWriter bw= new BufferedWriter
				(new FileWriter(out_fil, false));
		for (int i=0;i< freqArray.size();i++ ) {
			String ss ="Viral Quasispecies - strain1_fre : "+ Double.toString(freqArray.get(i));
			bw.write(ss+"\n");
			
			bw.write(seqArray.get(i).substring(start-1, end)+"\n");
		}
		
		bw.close();
		
		
		
	}
	
	public void generate_inter_hap(DivideConquer dc) throws NumberFormatException, IOException {
		BufferedReader br_snv = new BufferedReader(new FileReader(dc.out_dir+"/intermediate/Seq_SNV.txt"));
		String line="";
		Vector<Integer  > posVec = new Vector<Integer  >();
		this.Pos2SNV = new HashMap<Integer, String>(); 
		this.Pos2Allele = new HashMap<Integer, String>(); 
		while ((line = br_snv.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				String[] tmp = line.split("\t");
				int pos =  Integer.parseInt(tmp[0]);
				if  (! tmp[2].equals("*") ) {
					String ref = tmp[1];
					this.Pos2SNV.put(pos, ref); 
					posVec.add(pos) ; 
					this.Pos2Allele .put(pos, tmp[2]); 
				}
		}
		br_snv.close();
		this.SNVposArray =new int [posVec.size()];
		for (int i=0;i< posVec.size();i++) {
			this.SNVposArray[i] = posVec.get(i);
		}
		
		
		double hap_freq = 0.0;
		for (int i=0;i< dc.region_start.length;i++) {
			HashMap<String, Integer> PosIndex2geno =  new HashMap<String, Integer>(); 
			Vector<Double  > freqVec = new Vector<Double  >();
			String seq_fil = dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
					+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
			BufferedReader br_seq = new BufferedReader(new FileReader(seq_fil));
			int hap_count =0;
			while ((line = br_seq.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (line.substring(0,1).equals("V") ) {
					String[] tmp = line.split(":");
					hap_freq =  Double.parseDouble(tmp[1]);
				} else {
					if (hap_freq > dc.min_hap_freq) {
						freqVec.add(hap_freq); 
						
						for (int p=0; p< line.length(); p++) {
							int pos = dc.region_start[i] + p;
							String snp = line.substring(p,p+1);
							if  (this.Pos2SNV.containsKey(pos) ) {
								String key = Integer.toString(pos)+":"+Integer.toString(hap_count);
								int value = 9;
								if (snp.equals("*") ) {
									value = 9;
								} else {
									if (snp.equals(this.Pos2SNV.get(pos) ) ) {
										value = 0;
									} else {
										for (int s=0;s< this.Pos2Allele.get(pos).length();s++) {
											if ( Pos2Allele.get(pos).substring(s, s+1).equals(snp) ) { 
												value = s+1;
											}
										}
									}
								}
								PosIndex2geno.put(key, value);
							}
						}
						hap_count++;
					}
				}
			}
			br_seq.close();
			BufferedWriter bw_inter_freq = new BufferedWriter
					(new FileWriter(dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
							+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_inter.freq", false));
			String ss ="Hap_ID";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\th"+ Integer.toString(j+1);
			}
			bw_inter_freq.write(ss+"\n");
			ss ="Freq";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\t"+ Double.toString(freqVec.get(j));
			}
			bw_inter_freq.write(ss+"\n");
			for (int j=0; j< this.SNVposArray.length; j++) {
				if  (( this.SNVposArray[j]>=dc.region_start[i] ) && (this.SNVposArray[j]<=dc.region_end[i])) {
					ss= "0;"+ Integer.toString(this.SNVposArray[j])+";" + Integer.toString(this.SNVposArray[j])
						+";0:1";
					for (int h=0;h< freqVec.size();h++) {
						String key = Integer.toString(this.SNVposArray[j])+":"+ Integer.toString(h);
						ss =ss+"\t" + Integer.toString(PosIndex2geno.get(key));
					}
					bw_inter_freq.write(ss+"\n");
				}
			}
			
			bw_inter_freq.close();
		}
		
		for (int i=0;i< dc.tiling_region_start.length;i++) {
			HashMap<String, Integer> PosIndex2geno =  new HashMap<String, Integer>(); 
			Vector<Double  > freqVec = new Vector<Double  >();
			String seq_fil = dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
					+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
			BufferedReader br_seq = new BufferedReader(new FileReader(seq_fil));
			int hap_count =0;
			while ((line = br_seq.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (line.substring(0,1).equals("V") ) {
					String[] tmp = line.split(":");
					hap_freq =  Double.parseDouble(tmp[1]);
				} else {
					if (hap_freq > dc.min_hap_freq) {
						freqVec.add(hap_freq); 
						
						for (int p=0; p< line.length(); p++) {
							int pos = dc.tiling_region_start[i] + p;
							String snp = line.substring(p,p+1);
							if  (this.Pos2SNV.containsKey(pos) ) {
								String key = Integer.toString(pos)+":"+Integer.toString(hap_count);
								int value = 9;
								if (snp.equals("*") ) {
									value = 9;
								} else {
									if (snp.equals(this.Pos2SNV.get(pos) ) ) {
										value = 0;
									} else {
										for (int s=0;s< this.Pos2Allele.get(pos).length();s++) {
											if ( Pos2Allele.get(pos).substring(s, s+1).equals(snp) ) { 
												value = s+1;
											}
										}
									}
								}
								PosIndex2geno.put(key, value);
							}
						}
						hap_count++;
					}
				}
			}
			br_seq.close();
			BufferedWriter bw_inter_freq = new BufferedWriter
					(new FileWriter(dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
							+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_inter.freq", false));
			String ss ="Hap_ID";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\th"+ Integer.toString(j+1);
			}
			bw_inter_freq.write(ss+"\n");
			ss ="Freq";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\t"+ Double.toString(freqVec.get(j));
			}
			bw_inter_freq.write(ss+"\n");
			for (int j=0; j< this.SNVposArray.length; j++) {
				if  (( this.SNVposArray[j]>=dc.tiling_region_start[i] ) && (this.SNVposArray[j]<=dc.tiling_region_end[i])) {
					ss= "0;"+ Integer.toString(this.SNVposArray[j])+";" + Integer.toString(this.SNVposArray[j])
						+";0:1";
					for (int h=0;h< freqVec.size();h++) {
						String key = Integer.toString(this.SNVposArray[j])+":"+ Integer.toString(h);
						ss =ss+"\t" + Integer.toString(PosIndex2geno.get(key));
					}
					bw_inter_freq.write(ss+"\n");
				}
			}
			
			bw_inter_freq.close();
		}

		
	}
	
	
	public void generate_inter_hap(aBayesQR_DivideConquer dc) throws NumberFormatException, IOException {
		BufferedReader br_snv = new BufferedReader(new FileReader(dc.out_dir+"/intermediate/Seq_SNV.txt"));
		String line="";
		Vector<Integer  > posVec = new Vector<Integer  >();
		this.Pos2SNV = new HashMap<Integer, String>(); 
		this.Pos2Allele = new HashMap<Integer, String>(); 
		while ((line = br_snv.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				String[] tmp = line.split("\t");
				int pos =  Integer.parseInt(tmp[0]);
				if  (! tmp[2].equals("*") ) {
					String ref = tmp[1];
					this.Pos2SNV.put(pos, ref); 
					posVec.add(pos) ; 
					this.Pos2Allele .put(pos, tmp[2]); 
				}
		}
		br_snv.close();
		this.SNVposArray =new int [posVec.size()];
		for (int i=0;i< posVec.size();i++) {
			this.SNVposArray[i] = posVec.get(i);
		}
		
		
		double hap_freq = 0.0;
		for (int i=0;i< dc.region_start.length;i++) {
			HashMap<String, Integer> PosIndex2geno =  new HashMap<String, Integer>(); 
			Vector<Double  > freqVec = new Vector<Double  >();
			String seq_fil = dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
					+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
			BufferedReader br_seq = new BufferedReader(new FileReader(seq_fil));
			int hap_count =0;
			while ((line = br_seq.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (line.substring(0,1).equals("V") ) {
					String[] tmp = line.split(":");
					hap_freq =  Double.parseDouble(tmp[1]);
				} else {
					if (hap_freq > dc.min_hap_freq) {
						freqVec.add(hap_freq); 
						
						for (int p=0; p< line.length(); p++) {
							int pos = dc.region_start[i] + p;
							String snp = line.substring(p,p+1);
							if  (this.Pos2SNV.containsKey(pos) ) {
								String key = Integer.toString(pos)+":"+Integer.toString(hap_count);
								int value = 9;
								if (snp.equals("*") ) {
									value = 9;
								} else {
									if (snp.equals(this.Pos2SNV.get(pos) ) ) {
										value = 0;
									} else {
										for (int s=0;s< this.Pos2Allele.get(pos).length();s++) {
											if ( Pos2Allele.get(pos).substring(s, s+1).equals(snp) ) { 
												value = s+1;
											}
										}
									}
								}
								PosIndex2geno.put(key, value);
							}
						}
						hap_count++;
					}
				}
			}
			br_seq.close();
			BufferedWriter bw_inter_freq = new BufferedWriter
					(new FileWriter(dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
							+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_inter.freq", false));
			String ss ="Hap_ID";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\th"+ Integer.toString(j+1);
			}
			bw_inter_freq.write(ss+"\n");
			ss ="Freq";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\t"+ Double.toString(freqVec.get(j));
			}
			bw_inter_freq.write(ss+"\n");
			for (int j=0; j< this.SNVposArray.length; j++) {
				if  (( this.SNVposArray[j]>=dc.region_start[i] ) && (this.SNVposArray[j]<=dc.region_end[i])) {
					ss= "0;"+ Integer.toString(this.SNVposArray[j])+";" + Integer.toString(this.SNVposArray[j])
						+";0:1";
					for (int h=0;h< freqVec.size();h++) {
						String key = Integer.toString(this.SNVposArray[j])+":"+ Integer.toString(h);
						ss =ss+"\t" + Integer.toString(PosIndex2geno.get(key));
					}
					bw_inter_freq.write(ss+"\n");
				}
			}
			
			bw_inter_freq.close();
		}
		
		for (int i=0;i< dc.tiling_region_start.length;i++) {
			HashMap<String, Integer> PosIndex2geno =  new HashMap<String, Integer>(); 
			Vector<Double  > freqVec = new Vector<Double  >();
			String seq_fil = dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
					+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
			BufferedReader br_seq = new BufferedReader(new FileReader(seq_fil));
			int hap_count =0;
			while ((line = br_seq.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (line.substring(0,1).equals("V") ) {
					String[] tmp = line.split(":");
					hap_freq =  Double.parseDouble(tmp[1]);
				} else {
					if (hap_freq > dc.min_hap_freq) {
						freqVec.add(hap_freq); 
						
						for (int p=0; p< line.length(); p++) {
							int pos = dc.tiling_region_start[i] + p;
							String snp = line.substring(p,p+1);
							if  (this.Pos2SNV.containsKey(pos) ) {
								String key = Integer.toString(pos)+":"+Integer.toString(hap_count);
								int value = 9;
								if (snp.equals("*") ) {
									value = 9;
								} else {
									if (snp.equals(this.Pos2SNV.get(pos) ) ) {
										value = 0;
									} else {
										for (int s=0;s< this.Pos2Allele.get(pos).length();s++) {
											if ( Pos2Allele.get(pos).substring(s, s+1).equals(snp) ) { 
												value = s+1;
											}
										}
									}
								}
								PosIndex2geno.put(key, value);
							}
						}
						hap_count++;
					}
				}
			}
			br_seq.close();
			BufferedWriter bw_inter_freq = new BufferedWriter
					(new FileWriter(dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
							+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_inter.freq", false));
			String ss ="Hap_ID";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\th"+ Integer.toString(j+1);
			}
			bw_inter_freq.write(ss+"\n");
			ss ="Freq";
			for (int j=0; j< freqVec.size(); j++) {
				ss =ss+"\t"+ Double.toString(freqVec.get(j));
			}
			bw_inter_freq.write(ss+"\n");
			for (int j=0; j< this.SNVposArray.length; j++) {
				if  (( this.SNVposArray[j]>=dc.tiling_region_start[i] ) && (this.SNVposArray[j]<=dc.tiling_region_end[i])) {
					ss= "0;"+ Integer.toString(this.SNVposArray[j])+";" + Integer.toString(this.SNVposArray[j])
						+";0:1";
					for (int h=0;h< freqVec.size();h++) {
						String key = Integer.toString(this.SNVposArray[j])+":"+ Integer.toString(h);
						ss =ss+"\t" + Integer.toString(PosIndex2geno.get(key));
					}
					bw_inter_freq.write(ss+"\n");
				}
			}
			
			bw_inter_freq.close();
		}
	}
	
	
	
	public void write_hap(String snv_fil, String hap_freq_fil, String out_fil) 
			throws NumberFormatException, IOException, InterruptedException {
		HapConfig hc_tmp = new HapConfig (hap_freq_fil);
		
		
	}
	
	
	public void Variants2(String parameter_file, DivideConquer dc) throws IOException, InterruptedException {
		this.SNV2index = new HashMap<Integer, Integer>(); 
		ArrayList<Integer >  posArray = new ArrayList<Integer>();
		HashMap<Integer, String> alleleDic = new HashMap<Integer, String>();
		for (int i=0;i< dc.region_start.length;i++) {
			if (i==0) {
				
				int hap_len = dc.region_end[i] - dc.region_start[i]+1;
				ArrayList<String >  haploArray = new ArrayList<String>(); 
				String line ="";
				double hap_freq = 0.0;
				String til_hap_fil2= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
						+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
				BufferedReader br_til_hap2 = new BufferedReader(new FileReader(til_hap_fil2));
				while ((line = br_til_hap2.readLine()) != null) {
						line= line.replace("\n", "").replace("\r", "").replace(" ", "");
						if  (line.substring(0,1).equals("V") ) {
							String[] tmp = line.split(":");
							hap_freq =  Double.parseDouble(tmp[1]);
						} else {
							if (hap_freq> dc.min_hap_freq) {
								if (line.length() ==  ( dc.tiling_region_end[i]+1- dc.tiling_region_start[i])) {
									String ss= line.substring(0,  dc.region_end[i] -
											dc.tiling_region_start[i]+1);
									while (ss.length()< hap_len) {
										ss ="-"+ ss;
									}
									haploArray.add(ss);
								}
							}
						}
				}
				br_til_hap2.close();
				
				String hap_fil= dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
						+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
				BufferedReader br_hap = new BufferedReader(new FileReader(hap_fil));
				while ((line = br_hap.readLine()) != null) {
						line= line.replace("\n", "").replace("\r", "").replace(" ", "");
						if  (line.substring(0,1).equals("V") ) {
							String[] tmp = line.split(":");
							hap_freq =  Double.parseDouble(tmp[1]);
						} else {
							if (hap_freq> dc.min_hap_freq) {
								if (line.length() ==  ( dc.region_end[i]+1- dc.region_start[i])) {
									haploArray.add(line);
								}
							}
						}
				}
				br_hap.close();
				HashMap<Integer, Integer> pos_dic = new HashMap<Integer, Integer>();
		        String tiling_fil_2= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
						+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_SNV_pos.txt" ;
				BufferedReader br_2 = new BufferedReader(new FileReader(tiling_fil_2));
		        while ((line = br_2.readLine()) != null) {
		        	line= line.replace("\n", "").replace("\r", "");
		        	String[] tmp = line.split(" ");
		        	for (int t=0; t< tmp.length;t++) {
		        		int pos =Integer.parseInt(tmp[t] );
		        		boolean flag = true;;
		        		for (int j=0; j<haploArray.size();j++) {
		        			int p_hap= pos +dc.region_end[i] -
									dc.tiling_region_start[i] ;
		        			if (p_hap< hap_len) {
			        			if (haploArray.get(j).substring(p_hap, p_hap+1).equals("*")) {
			        				flag =false;
			        			}
		        			}
		        		}
		        		if (flag==true) {
		        			pos_dic.put(Integer.parseInt(tmp[t] )+ dc.tiling_region_start[i],1);
		        		}
		        	}
		        }
		        br_2.close();
		        
		        String fil = dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
						+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_SNV_pos.txt" ;
		        BufferedReader br = new BufferedReader(new FileReader(fil));
		        while ((line = br.readLine()) != null) {
		        	line= line.replace("\n", "").replace("\r", "");
		        	String[] tmp = line.split(" ");
		        	for (int t=0; t< tmp.length;t++) {
		        		int pos =Integer.parseInt(tmp[t] );
		        		
		        		if (pos_dic.containsKey(Integer.parseInt(tmp[t] )+ dc.region_start[i])) {
		        			boolean flag = true;
		        			String allele ="";
		        			for (int j=0; j<haploArray.size();j++) {
			        			int p_hap= pos ;
			        			if (p_hap< hap_len) {
			        				allele= allele+haploArray.get(j).substring(pos, pos+1);
				        			if (haploArray.get(j).substring(pos, pos+1).equals("*")) {
				        				flag =false;
				        			}
			        			}
			        		}
			        		if (flag==true) {
			        			posArray.add(Integer.parseInt(tmp[t] )+ dc.region_start[i] );
			        			alleleDic.put(Integer.parseInt(tmp[t] )+ dc.region_start[i], allele);
			        		} 
		        		}else {
		        			if  ((Integer.parseInt(tmp[t] )+ dc.region_start[i]) < 
		        					dc.tiling_region_start[i]) {
		        				
		        				boolean flag = true;
		        				String allele ="";
			        			for (int j=0; j<haploArray.size();j++) {
				        			int p_hap= pos ;
				        			if (p_hap< hap_len) {
				        				allele= allele+haploArray.get(j).substring(pos, pos+1);
					        			if (haploArray.get(j).substring(pos, pos+1).equals("*")) {
					        				flag =false;
					        			}
				        			}
				        		}
			        			if (flag==true) {
				        		alleleDic.put(Integer.parseInt(tmp[t] )+ dc.region_start[i], allele);
				        		posArray.add(Integer.parseInt(tmp[t] )+ dc.region_start[i] );
			        			}
	        				
		        			}
	        			}
		        	}
		        }
		        br.close();

		        
			}else if (i==(dc.region_start.length-1)) {
				int hap_len = dc.region_end[i] - dc.region_start[i]+1;
				ArrayList<String >  haploArray = new ArrayList<String>(); 
				String til_hap_fil1= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i-1])+"_"
						+ Integer.toString(dc.tiling_region_end[i-1])+"/"+ dc.project_name+"_ViralSeq.txt";
					String line =null;
					double hap_freq = 0.0;
					BufferedReader br_til_hap1 = new BufferedReader(new FileReader(til_hap_fil1));
					while ((line = br_til_hap1.readLine()) != null) {
						line= line.replace("\n", "").replace("\r", "").replace(" ", "");
						if  (line.substring(0,1).equals("V") ) {
							String[] tmp = line.split(":");
							hap_freq =  Double.parseDouble(tmp[1]);
						} else {
							if (hap_freq> dc.min_hap_freq) {
								if (line.length() ==  ( dc.tiling_region_end[i-1]+1- dc.tiling_region_start[i-1])) {
									String ss= line.substring( dc.region_start[i] - 
											dc.tiling_region_start[i-1],dc.tiling_region_end[i-1]- 
											dc.tiling_region_start[i-1] +1);
									while (ss.length()< hap_len) {
										ss = ss +"-";
									}
									haploArray.add(ss);
								}
							}
						}
					}
					br_til_hap1.close();
					
					String hap_fil= dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
							+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
					BufferedReader br_hap = new BufferedReader(new FileReader(hap_fil));
					while ((line = br_hap.readLine()) != null) {
							line= line.replace("\n", "").replace("\r", "").replace(" ", "");
							if  (line.substring(0,1).equals("V") ) {
								String[] tmp = line.split(":");
								hap_freq =  Double.parseDouble(tmp[1]);
							} else {
								if (hap_freq> dc.min_hap_freq) {
									if (line.length() ==  ( dc.region_end[i]+1- dc.region_start[i])) {
										haploArray.add(line);
									}
								}
							}
					}
					br_hap.close();
					
					
					HashMap<Integer, Integer> pos_dic = new HashMap<Integer, Integer>();
					String tiling_fil_1= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i-1])+"_"
							+ Integer.toString(dc.tiling_region_end[i-1])+"/"+ dc.project_name+"_SNV_pos.txt" ;
					BufferedReader br_1 = new BufferedReader(new FileReader(tiling_fil_1));
			        while ((line = br_1.readLine()) != null) {
			        	line= line.replace("\n", "").replace("\r", "");
			        	String[] tmp = line.split(" ");
			        	for (int t=0; t< tmp.length;t++) {
			        		int pos =Integer.parseInt(tmp[t] );
			        		boolean flag = true;;
			        		for (int j=0; j<haploArray.size();j++) {
			        			int p_hap= pos - (dc.region_start[i] - 
										dc.tiling_region_start[i-1]);
			        			if (p_hap>=0 ) {
				        			if (haploArray.get(j).substring(p_hap, p_hap+1).equals("*")) {
				        				flag =false;
				        			}
			        			}
			        		}
			        		if (flag==true) {
			        			pos_dic.put(Integer.parseInt(tmp[t] )+ dc.tiling_region_start[i-1], 1);
			        		} 
			        	}
			        }
			        br_1.close();
			        								
			        String fil = dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
							+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_SNV_pos.txt" ;
			        BufferedReader br = new BufferedReader(new FileReader(fil));
			        while ((line = br.readLine()) != null) {
			        	line= line.replace("\n", "").replace("\r", "");
			        	String[] tmp = line.split(" ");
			        	for (int t=0; t< tmp.length;t++) {
			        		int pos =Integer.parseInt(tmp[t] );
			        		if (pos_dic.containsKey(Integer.parseInt(tmp[t] )+ dc.region_start[i] )) {
			        			boolean flag = true;;
			        			String allele ="";
			        			for (int j=0; j<haploArray.size();j++) {
				        			int p_hap= pos ;
				        			if (p_hap< hap_len) {
				        				allele= allele+haploArray.get(j).substring(pos, pos+1);
					        			if (haploArray.get(j).substring(pos, pos+1).equals("*")) {
					        				flag =false;
					        			}
				        			}
				        		}
				        		if (flag==true) {
				        			posArray.add(Integer.parseInt(tmp[t] )+ dc.region_start[i] );
				        			alleleDic.put(Integer.parseInt(tmp[t] )+ dc.region_start[i], allele);
				        		} 
			        		}else {
			        			if ((Integer.parseInt(tmp[t] )+ dc.region_start[i]) > dc.tiling_region_end[i-1]) {
			        				boolean flag = true;
			        				String allele ="";
				        			for (int j=0; j<haploArray.size();j++) {
					        			int p_hap= pos ;
					        			if (p_hap< hap_len) {
					        				allele= allele+haploArray.get(j).substring(pos, pos+1);
						        			if (haploArray.get(j).substring(pos, pos+1).equals("*")) {
						        				flag =false;
						        			}
					        			}
					        		}
				        			if (flag==true) {
						        		alleleDic.put(Integer.parseInt(tmp[t] )+ dc.region_start[i], allele);
						        		posArray.add(Integer.parseInt(tmp[t] )+ dc.region_start[i] );
					        		}
					        		
					        		
			        			}
			        		}
			        	}
			        }
			        br.close();
			}else {
				int hap_len = dc.region_end[i] - dc.region_start[i]+1;
				ArrayList<String >  haploArray = new ArrayList<String>(); 
				String til_hap_fil1= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i-1])+"_"
					+ Integer.toString(dc.tiling_region_end[i-1])+"/"+ dc.project_name+"_ViralSeq.txt";
				String line =null;
				double hap_freq = 0.0;
				BufferedReader br_til_hap1 = new BufferedReader(new FileReader(til_hap_fil1));
				while ((line = br_til_hap1.readLine()) != null) {
					line= line.replace("\n", "").replace("\r", "").replace(" ", "");
					if  (line.substring(0,1).equals("V") ) {
						String[] tmp = line.split(":");
						hap_freq =  Double.parseDouble(tmp[1]);
					} else {
						if (hap_freq> dc.min_hap_freq) {
							if (line.length() ==  ( dc.tiling_region_end[i-1]+1- dc.tiling_region_start[i-1])) {
								String ss= line.substring( dc.region_start[i] - 
										dc.tiling_region_start[i-1], dc.tiling_region_end[i-1]- 
										dc.tiling_region_start[i-1] +1);
								while (ss.length()< hap_len) {
									ss =ss+"-";
								}
								haploArray.add(ss);
							}
						}
					}
				}
				br_til_hap1.close();
				
				String til_hap_fil2= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
						+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
				BufferedReader br_til_hap2 = new BufferedReader(new FileReader(til_hap_fil2));
				while ((line = br_til_hap2.readLine()) != null) {
						line= line.replace("\n", "").replace("\r", "").replace(" ", "");
						if  (line.substring(0,1).equals("V") ) {
							String[] tmp = line.split(":");
							hap_freq =  Double.parseDouble(tmp[1]);
						} else {
							if (hap_freq> dc.min_hap_freq) {
								if (line.length() ==  ( dc.tiling_region_end[i]+1- dc.tiling_region_start[i])) {
									String ss= line.substring(0,  dc.region_end[i] -
											dc.tiling_region_start[i]+1);
									while (ss.length()< hap_len) {
										ss ="-"+ ss;
									}
									haploArray.add(ss);
								}
							}
						}
				}
				br_til_hap2.close();
				
				
				String hap_fil= dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
						+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_ViralSeq.txt";
				BufferedReader br_hap = new BufferedReader(new FileReader(hap_fil));
				while ((line = br_hap.readLine()) != null) {
						line= line.replace("\n", "").replace("\r", "").replace(" ", "");
						if  (line.substring(0,1).equals("V") ) {
							String[] tmp = line.split(":");
							hap_freq =  Double.parseDouble(tmp[1]);
						} else {
							if (hap_freq> dc.min_hap_freq) {
								if (line.length() ==  ( dc.region_end[i]+1- dc.region_start[i])) {
									haploArray.add(line);
								}
							}
						}
				}
				br_hap.close();
				
				HashMap<Integer, Integer> pos_dic = new HashMap<Integer, Integer>();
				
				String tiling_fil_1= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i-1])+"_"
						+ Integer.toString(dc.tiling_region_end[i-1])+"/"+ dc.project_name+"_SNV_pos.txt" ;
				BufferedReader br_1 = new BufferedReader(new FileReader(tiling_fil_1));
		        while ((line = br_1.readLine()) != null) {
		        	line= line.replace("\n", "").replace("\r", "");
		        	String[] tmp = line.split(" ");
		        	for (int t=0; t< tmp.length;t++) {
		        		int pos =Integer.parseInt(tmp[t] );
		        		
		        		boolean flag = true;;
		        		for (int j=0; j<haploArray.size();j++) {
		        			int p_hap= pos - (dc.region_start[i] - 
									dc.tiling_region_start[i-1]);
		        			if (p_hap>=0 ) {
			        			if (haploArray.get(j).substring(p_hap, p_hap+1).equals("*")) {
			        				flag =false;
			        			}
		        			}
		        		}
		        		if (flag==true) {
		        			pos_dic.put(Integer.parseInt(tmp[t] )+ dc.tiling_region_start[i-1], 1);
		        		}
		        	}
		        }
		        br_1.close();
		        
		        String tiling_fil_2= dc.out_dir+"/intermediate/"+ Integer.toString(dc.tiling_region_start[i])+"_"
						+ Integer.toString(dc.tiling_region_end[i])+"/"+ dc.project_name+"_SNV_pos.txt" ;
				BufferedReader br_2 = new BufferedReader(new FileReader(tiling_fil_2));
		        while ((line = br_2.readLine()) != null) {
		        	line= line.replace("\n", "").replace("\r", "");
		        	String[] tmp = line.split(" ");
		        	for (int t=0; t< tmp.length;t++) {
		        		int pos =Integer.parseInt(tmp[t] );
		        		boolean flag = true;;
		        		for (int j=0; j<haploArray.size();j++) {
		        			int p_hap= pos +dc.region_end[i] -
									dc.tiling_region_start[i] ;
		        			if (p_hap< hap_len) {
			        			if (haploArray.get(j).substring(p_hap, p_hap+1).equals("*")) {
			        				flag =false;
			        			}
		        			}
		        		}
		        		if (flag==true) {
		        			pos_dic.put(Integer.parseInt(tmp[t] )+ dc.tiling_region_start[i], 1);
		        		}
		        	}
		        }
		        br_2.close();
		        
		        String fil = dc.out_dir+"/intermediate/"+ Integer.toString(dc.region_start[i])+"_"
						+ Integer.toString(dc.region_end[i])+"/"+ dc.project_name+"_SNV_pos.txt" ;
		        BufferedReader br = new BufferedReader(new FileReader(fil));
		        while ((line = br.readLine()) != null) {
		        	line= line.replace("\n", "").replace("\r", "");
		        	String[] tmp = line.split(" ");
		        	for (int t=0; t< tmp.length;t++) {
		        		int pos =Integer.parseInt(tmp[t] );
		        		if (pos_dic.containsKey(Integer.parseInt(tmp[t] )+ dc.region_start[i] )) {
		        			boolean flag = true;
		        			String allele= "";
		        			for (int j=0; j<haploArray.size();j++) {
			        			int p_hap= pos ;
			        			if (p_hap< hap_len) {
			        				allele=allele+ haploArray.get(j).substring(pos, pos+1);
				        			if (haploArray.get(j).substring(pos, pos+1).equals("*")) {
				        				flag =false;
				        			}
			        			}
			        		}
			        		if (flag==true) {
			        			posArray.add(Integer.parseInt(tmp[t] )+ dc.region_start[i] );
			        			alleleDic.put(Integer.parseInt(tmp[t] )+ dc.region_start[i] , allele); 
			        		}
		        		}
		        	}
		        }
		        br.close();
			}
		}
		
		this.SNV_posArray = new int [posArray.size()];
		for (int i=0;i< posArray.size(); i++) {
			this.SNV2index.put(posArray.get(i), i);
			this.SNV_posArray[i]= posArray.get(i);
		}
		
		
		String ref="";
		String line = null;
		
		BufferedReader br_ref = new BufferedReader(new FileReader(dc.reference));
		while ((line = br_ref.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				if  (!line.substring(0,1).equals(">") ) {
					ref=ref+ line;
				} 
		}
		br_ref.close();
//		System.out.println(ref);
		
		BufferedWriter bw = new BufferedWriter
				(new FileWriter(String.valueOf(dc.out_dir)+"/intermediate/"
						+dc.project_name+".Variants.txt", false));
		String ss ="START\tEND\tREF\tALT";
		bw.write(ss+"\n");
		
		for (int i=0; i<this.SNV_posArray.length;i++ ) {
			if (alleleDic.containsKey(this.SNV_posArray[i])){
				ss= alleleDic.get(this.SNV_posArray[i]); 
			}
			String ref_allele = ref.substring(this.SNV_posArray[i]-1, this.SNV_posArray[i]);
			bw.write(Integer.toString(this.SNV_posArray[i])+"\t"+ss+"\t"+ 
					ref.substring(this.SNV_posArray[i]-1, this.SNV_posArray[i])+"\n");
			
		}
		bw.close();
	}
}
