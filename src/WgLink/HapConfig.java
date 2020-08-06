package virushaplo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.HashSet;
import java.time.temporal.TemporalAccessor;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.io.InputStreamReader;

public class HapConfig {
	public  HashMap<String, String> Index2SNV;
	public  HashMap<Integer, Integer> Index2Pos;
	public  HashMap<String, String> Pos2SNV;
	public String [] Haps;
	public double [] Freqs;
	public String [] Haps_names;
	public  HashSet<Integer> AcceptPos;
	
	
	public HapConfig( ) throws IOException, InterruptedException {
		
	}

	public HapConfig( String fil) throws IOException, InterruptedException {
		this.Pos2SNV = new HashMap<String, String> ();
		this.Index2Pos = new HashMap<Integer, Integer> ();
		this.Index2SNV = new HashMap<String, String> ();
		this.AcceptPos = new HashSet <Integer>  ();
		String line ="";
		BufferedReader br_hap = new BufferedReader(new FileReader(fil));
		line = br_hap.readLine();
		line= line.replace("\n", "").replace("\r", "").replace(" ", "");
		String[] tmp2 = line.split("\t");
		this.Haps = new String [tmp2.length-1];
		this.Freqs = new double [tmp2.length-1];
		this.Haps_names = new String [tmp2.length-1];
		
		for (int i=0;i< this.Haps.length;i++) {
			this.Haps[i]="";
			this.Freqs[i]=0.0;
			this.Haps_names[i] = tmp2[i+1];
		}
		int index = 0;
		while ((line = br_hap.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
			if (line.startsWith("0")) {
				String info = tmp[0];
				int pos = Integer.parseInt(info.split(";")[1]);
				boolean ok =true;
				for (int i=1;i< tmp.length;i++) {
					this.Haps[i-1]= this.Haps[i-1] + tmp[i];
					String snv = tmp[i];
					this.Index2Pos.put(index, pos);
					this.Index2SNV.put(Integer.toString(index)+":"+ Integer.toString(i-1), snv);
					this.Pos2SNV.put(Integer.toString(pos)+":"+ Integer.toString(i-1), snv);
					if ( snv.equals("9") ) {
						ok =false;
					}
				}
				if (ok) {
				this.AcceptPos.add(pos);
				}
				index++;
			} else if (line.startsWith("F")) {
				for (int i=1;i< tmp.length;i++) {
					this.Freqs[i-1]= Double.parseDouble(tmp[i]); 
					
				}
			}
		}
		br_hap.close();
	}
	
	public void write2file( HapConfig hc, String fil) throws IOException, InterruptedException {	
		BufferedWriter bw_hap = new BufferedWriter
				(new FileWriter(fil, false));
		String ss ="Hap_ID";
		for (int j=0; j< hc.Freqs.length; j++) {
			ss =ss+"\th"+ Integer.toString(j+1);
		}
		bw_hap.write(ss+"\n");
		ss ="Freq";
		for (int j=0; j< hc.Freqs.length ; j++) {
			ss =ss+"\t"+ Double.toString(hc.Freqs[j]);
		}
		
		bw_hap.write(ss+"\n");
		for (int p=0;p< hc.Haps[0].length();p++ ) {
			ss = "0;" + Integer.toString(hc.Index2Pos.get(p))+ ";" + Integer.toString(hc.Index2Pos.get(p)) 
				+";0:1"; 
			for (int h=0; h< hc.Haps.length ;h++ ) {
				String key = Integer.toString(p)+":"+   Integer.toString( h);
				ss= ss+ "\t"+ hc.Index2SNV.get(key);
			}
			bw_hap.write(ss+"\n");
		}
		bw_hap.close();
	}
	
	
	
	@SuppressWarnings("unlikely-arg-type")
	public HapConfig merge2regions( HapConfig hc1,HapConfig hc2, HapConfig hc3 , double cutoff, double min_freq,
			int ave_haps)
			throws IOException, InterruptedException {
		
		double [] hc2_tmp = new double[hc2.Freqs.length];
		for (int i=0;i< hc2.Freqs.length; i++ )  {
			hc2_tmp[i]= hc2.Freqs[i];
		}
		for (int i=0;i< hc2_tmp.length;i++) { 
			for (int j=i;j< hc2_tmp.length;j++) {
				if (hc2_tmp[j]> hc2_tmp[i]) { 
					double tm_value=  hc2_tmp[i];
					hc2_tmp[i]= hc2_tmp[j];
					hc2_tmp[j]= tm_value;
					
				}
			}
		}
		if (ave_haps>= hc2_tmp.length) {
			min_freq =0;
		} else {
			min_freq= hc2_tmp[ave_haps-1];
		}
		
		
		
		HapConfig hc =new HapConfig();
		boolean [] ismergedArray = new boolean [hc1.Haps.length* hc2.Haps.length];
		double [] hc1_misratio = new double [hc1.Haps.length* hc2.Haps.length];
		double [] hc2_misratio = new double [hc1.Haps.length* hc2.Haps.length];
		
		double [] hc1_misratio_sec = new double [hc1.Haps.length* hc2.Haps.length];
		double [] hc2_misratio_sec = new double [hc1.Haps.length* hc2.Haps.length];
		
		
		
		double [][][] ratio1_arr = new double [hc1.Haps.length] [hc2.Haps.length][ hc3.Haps.length];
		double [][][] ratio2_arr = new double [hc1.Haps.length][ hc2.Haps.length][ hc3.Haps.length];
		
		
		
		int [] h2_count = new int[hc2.Haps.length];
		for (int h2=0;h2< hc2.Haps.length; h2++) {
			h2_count[h2]=0;
		}
		
		int count=-1;
		int count_extention1 = 0;
		int count_extention2 = 0;
		for (int h1=0; h1< hc1.Haps.length ;h1++) {
			int num_extention= 0;
			
			for (int h2=0;h2< hc2.Haps.length; h2++) {
				count++;
				hc1_misratio[count]=1.0;
				hc2_misratio[count]=1.0;
				hc1_misratio_sec[count]=1.0;
				hc2_misratio_sec[count]=1.0;
				
				boolean ismerged =false;
				
				for (int h3=0; h3< hc3.Haps.length; h3++) {
					
					double mismatch_1= 0.0;
					double overlap_snp_1=0.01;
					double mismatch_2= 0.0;
					double overlap_snp_2=0.01;
					for (int i=0;i< hc1.Haps[h1].length();i++ ) {
						if  ( hc3.Pos2SNV.containsKey( Integer.toString( hc1.Index2Pos.get(i)) +":"
								+ Integer.toString(h3) ) ) {
							String hc3_snv = hc3.Pos2SNV.get(Integer.toString( hc1.Index2Pos.get(i)) +":"
									+ Integer.toString(h3) ) ;
							
							String hc1_snv = hc1.Index2SNV.get( Integer.toString(i) + ":"+ 
									Integer.toString(h1) ) ;
							if ( ( ! hc3_snv.equals("9") )  &&  ( ! hc1_snv.equals("9") ) ) {
								if (  Integer.parseInt( hc3_snv )  >0 ) {
									hc3_snv="1";
								} 
								if (  Integer.parseInt( hc1_snv )  >0 ) {
									hc1_snv="1";
								} 
								overlap_snp_1=overlap_snp_1+1;
								if  ( ! hc3_snv.equals(hc1_snv ) ) {
									mismatch_1=mismatch_1+1.0;
								}
							}
						}
					}
					
					for (int i=0;i< hc2.Haps[h2].length();i++ ) {
						if  ( hc3.Pos2SNV.containsKey( Integer.toString( hc2.Index2Pos.get(i)) +":"
								+ Integer.toString(h3) ) ) {
							String hc3_snv = hc3.Pos2SNV.get(Integer.toString( hc2.Index2Pos.get(i)) +":"
									+ Integer.toString(h3) ) ;
							
							String hc2_snv = hc2.Index2SNV.get( Integer.toString(i) + ":"+ 
									Integer.toString(h2) ) ;
							if ( ( ! hc3_snv.equals("9") )  &&  ( ! hc2_snv.equals("9") ) ) {
								if (  Integer.parseInt( hc3_snv )  >0 ) {
									hc3_snv="1";
								} 
								if (  Integer.parseInt( hc2_snv )  >0 ) {
									hc2_snv="1";
								} 
								overlap_snp_2=overlap_snp_2+1;
								if  ( ! hc3_snv.equals(hc2_snv ) ) {
									mismatch_2=mismatch_2+1.0;
								}
							}
						}
					}
					
					
					double ratio_1 = mismatch_1/  overlap_snp_1;
					double ratio_2 = mismatch_2/  overlap_snp_2;
					ratio1_arr[h1][h2][h3] = ratio_1;
					ratio2_arr[h1][h2][h3] = ratio_2;
					
					if (ratio_2< hc2_misratio_sec[count] ) {
						hc1_misratio_sec[count]=ratio_1;
						hc2_misratio_sec[count]=ratio_2;
					}
					if (hc2_misratio_sec[count]< hc2_misratio[count] ) {
						double tmp_value = hc1_misratio_sec[count];
						hc1_misratio_sec[count]= hc1_misratio[count];
						hc1_misratio[count]= tmp_value;
						tmp_value = hc2_misratio_sec[count];
						hc2_misratio_sec[count]= hc2_misratio[count];
						hc2_misratio[count]= tmp_value;
						
					}
					
						
					
					
					if ((ratio_1 < cutoff) && (ratio_2<cutoff) )  {
						double hap_proportion = hc1.Freqs[h1]/ hc2.Freqs[h2];
						if (hap_proportion>0) {
							ismerged= true;
							num_extention++;
							break;
						}
					}
				}
				ismergedArray[count] = ismerged;
				if (ismerged) {
					h2_count[h2]= h2_count[h2]+1;
				}
			}
			
			if (num_extention==0)  {
				double min_h1=1.0;
				double min_h2=1.0;
				double min_h1_sec=1.0;
				double min_h2_sec=1.0;
				
				for (int h2=0; h2<hc2.Haps.length;h2++ ) {
					for (int h3=0;h3< hc3.Haps.length;h3++) {
						if ( ratio1_arr[h1][h2][h3]< min_h1 ) {
							min_h1=ratio1_arr[h1][h2][h3];
						}
					}
				}
				

				
				for (int h2=0; h2<hc2.Haps.length;h2++ ) {
					for (int h3=0;h3< hc3.Haps.length;h3++) {
						if (Double.compare(ratio1_arr[h1][h2][h3],min_h1 ) == 0 ) { 
							if (ratio2_arr[h1][h2][h3]< min_h2) {
								min_h2=ratio2_arr[h1][h2][h3];
							}
						}
					}
				}
				
				
				for (int h2=0; h2<hc2.Haps.length;h2++ ) {
					for (int h3=0;h3< hc3.Haps.length;h3++) {
						if ( (Double.compare(ratio1_arr[h1][h2][h3],min_h1 ) == 0 )   && 
						(Double.compare(ratio2_arr[h1][h2][h3],min_h2 ) == 0 )  ) { 
							ismergedArray[count-h2] =true;
							h2_count[h2]= h2_count[h2]+1;
						}
					}
				}
				
				for (int h2=0; h2<hc2.Haps.length;h2++ ) {
					for (int h3=0;h3< hc3.Haps.length;h3++) {
						if (( ratio1_arr[h1][h2][h3]< min_h1_sec ) &&
						( ratio1_arr[h1][h2][h3]> min_h1  )  ){
							min_h1_sec=ratio1_arr[h1][h2][h3];
						}
					}
				}
				
				for (int h2=0; h2<hc2.Haps.length;h2++ ) {
					for (int h3=0;h3< hc3.Haps.length;h3++) {
						if (Double.compare(ratio1_arr[h1][h2][h3],min_h1_sec ) == 0 ) { 
							if (ratio2_arr[h1][h2][h3]< min_h2_sec) {
								min_h2_sec=ratio2_arr[h1][h2][h3];
							}
						}
					}
				}
				
				for (int h2=0; h2<hc2.Haps.length;h2++ ) {
					for (int h3=0;h3< hc3.Haps.length;h3++) {
						if ( (Double.compare(ratio2_arr[h1][h2][h3],min_h2_sec ) == 0 )   && 
						(Double.compare(ratio1_arr[h1][h2][h3],min_h1_sec ) == 0 )  ) { 
							ismergedArray[count-h2] =true;
							h2_count[h2]= h2_count[h2]+1;
						}
					}
				}	
			}
			
		}
		
		
		boolean [] [] flag_h1h2=  new boolean [ hc1.Haps.length][hc2.Haps.length];
		for (int h2=0; h2<hc2.Haps.length;h2++ ) {
			for (int h1=0; h1<hc1.Haps.length;h1++ ) {
				flag_h1h2[h1][h2] = false;
			}
		}
		
		for (int h2=0; h2<hc2.Haps.length;h2++ ) {
			double min_h1=1.0;
			double min_h2=1.0;
			double min_h1_sec=1.0;
			double min_h2_sec=1.0;
			for (int h1=0; h1<hc1.Haps.length;h1++ ) {
				for (int h3=0;h3< hc3.Haps.length;h3++) {
					if ( ratio2_arr[h1][h2][h3]< min_h2 ) {
						min_h2=ratio2_arr[h1][h2][h3];
					}
				}
			}
			
			
			
			
			for (int h1=0; h1<hc1.Haps.length;h1++ ) {
				for (int h3=0;h3< hc3.Haps.length;h3++) {
					if (Double.compare(ratio2_arr[h1][h2][h3],min_h2 ) == 0 ) { 
						if (ratio1_arr[h1][h2][h3]< min_h1) {
							min_h1=ratio1_arr[h1][h2][h3];
						}
					}
				}
			}
			
			
			for (int h1=0; h1<hc1.Haps.length;h1++ ) {
				for (int h3=0;h3< hc3.Haps.length;h3++) {
					if ( (Double.compare(ratio2_arr[h1][h2][h3],min_h2 ) == 0 )   && 
					(Double.compare(ratio1_arr[h1][h2][h3],min_h1 ) == 0 )  ) { 
						flag_h1h2[h1][h2]=true;
					}
				}
			}
			
			for (int h1=0; h1<hc1.Haps.length;h1++ ) {
				for (int h3=0;h3< hc3.Haps.length;h3++) {
					if (( ratio2_arr[h1][h2][h3]< min_h2_sec ) &&
					( ratio2_arr[h1][h2][h3]> min_h2  )  ){
						min_h2_sec=ratio2_arr[h1][h2][h3];
					}
				}
			}
			
			for (int h1=0; h1<hc1.Haps.length;h1++ ) {
				for (int h3=0;h3< hc3.Haps.length;h3++) {
					if (Double.compare(ratio2_arr[h1][h2][h3],min_h2_sec ) == 0 ) { 
						if (ratio1_arr[h1][h2][h3]< min_h1_sec) {
							min_h1_sec=ratio1_arr[h1][h2][h3];
						}
					}
				}
			}
			
			for (int h1=0; h1<hc1.Haps.length;h1++ ) {
				for (int h3=0;h3< hc3.Haps.length;h3++) {
					if ( (Double.compare(ratio2_arr[h1][h2][h3],min_h2_sec ) == 0 )   && 
					(Double.compare(ratio1_arr[h1][h2][h3],min_h1_sec ) == 0 )  ) { 
						flag_h1h2[h1][h2]=true;
					}
				}
			}
			
			
		}
		
		
		
		
		count=-1;
		for (int h1=0; h1< hc1.Haps.length ;h1++) {
			for (int h2=0;h2< hc2.Haps.length; h2++) {
				count++;
				if ( (h2_count[h2]==0)  && ( hc2.Freqs[h2]>= (min_freq) ) )  {
					if (flag_h1h2[h1][h2] ) {
						ismergedArray[count] = true;
					}
				}
			}
		}
		
		int num_haps = 0;
		for (int i=0;i< ismergedArray.length;i++) {
			if ( ismergedArray[i] ) {
				num_haps++;
			}
		}
		final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
		System.out.println("Region "+ hc1.Index2Pos.get(0) +"_"+ hc1.Index2Pos.get(hc1.Index2Pos.size()-1 ) +" & "+ 
				"region "+ hc2.Index2Pos.get(0) +"_"+ hc2.Index2Pos.get(hc2.Index2Pos.size()-1 ) +" are merged by:\t"+ 
				"region "+ hc3.Index2Pos.get(0) +"_"+ hc3.Index2Pos.get(hc3.Index2Pos.size()-1 ) + "\t"+
				dtf.format(LocalDateTime.now()) );
		
//		System.out.println(num_haps);
//		System.out.println("Extention1:\t"+ count_extention1);
//		System.out.println("Extention2:\t"+ count_extention2);
		
		hc.Haps = new String [num_haps];
		hc.Freqs = new double [num_haps];
		hc.Index2Pos = new HashMap<Integer, Integer> ();
		hc.Index2SNV =new HashMap<String, String> ();
		hc.Pos2SNV  =new HashMap<String, String> ();
		
		int index = -1;
		count =-1;
		for (int h1=0; h1< hc1.Haps.length ;h1++) {
			for (int h2=0;h2< hc2.Haps.length; h2++) {
				count++;
				if  (ismergedArray[count] ) {
					
					index++;
					hc.Freqs[index] = 1.0/(double) num_haps;
					hc.Haps[index] =  hc1.Haps[h1]+ hc2.Haps[h2];
					
					for (int i=0;i< hc1.Haps[h1].length();i++ ) {
						String snv = hc1.Haps[h1].substring(i, i+1);
						int pos= hc1.Index2Pos.get( i );
						String pindex = Integer.toString(i);
						String hindex = Integer.toString(index);
						hc.Index2Pos.put(i ,pos);
						hc.Index2SNV.put(pindex +":"+hindex,snv);
						hc.Pos2SNV.put(Integer.toString(pos)+":"+Integer.toString(index) ,snv);
					}
					
					for (int i=0;i< hc2.Haps[h2].length();i++ ) {
						String snv = hc2.Haps[h2].substring(i, i+1);
						int pos= hc2.Index2Pos.get( i);
						String pindex = Integer.toString(i+ hc1.Haps[0].length());
						String hindex = Integer.toString(index);
						hc.Index2Pos.put(i+ hc1.Haps[0].length(),pos);
						hc.Index2SNV.put(pindex +":"+hindex ,snv);
						hc.Pos2SNV.put(Integer.toString(pos)+":"+Integer.toString(index) ,snv);
					}
					
				}
			}
		}
		
		return hc;
	}

}

