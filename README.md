## Users’ Manual of WgLink (Version 1.0)
### Preamble
The WgLink program reconstructs haplotypes.


The workflow of WgLink is briefed as follows: (A) WgLink first determines local region & tiling region haplotypes by applying tools: TenSQR or aBayesQR. 
(B) Regional haplotypes are joined together using breadth first search which could tolerate specified SNV mismatch when extending regions. C) The
haplotypes is estimated by L0L1 regularized regression to solve for the haplotype distribution problem.

### Installation
WgLink.jar is a batteries-included JAR executable. All needed external jar packages are included in the downloadable, PoolHapX.jar. To download all necessary files, users can use the command 
`git clone https://github.com/theLongLab/WgLink.git`
As we used an R package L0Learn, the users have to install R and L0Learn (https://cran.r-project.org/web/packages/L0Learn/index.html). The versions of R and R package L0Learn that we have used on our platform are: version 1.2.0 for L0Learn and version 3.6.1 for R. 

Several other tools are prerequisites for running WgLink. Users can download and install them from the websites:
* bwa: https://github.com/lh3/bwa
* python3
* ExtractMatrix (if using TenSQR calculating local region haplotype): https://github.com/SoYeonA/TenSQR
* TenSQR.py and its prerequisites (if using TenSQR calculating local region haplotype): https://github.com/SoYeonA/TenSQR
* aBayesQR (if using aBayesQR calculating local region haplotype) : https://github.com/SoYeonA/aBayesQR

### Functions
* TenSQR: Using TenSQR as local region haplotype calculation tool
* aBayesQR: Using aBayesQR as local region haplotype calculation tool


### Quick start with included example data

Example data is provided in the "Example" folder. Users can see the reference, fastq files under the “Example” folder. After updating absolute paths of executable (such as bwa etc) and parent folder in the TenSQR.config file, users can run WgLink by a simple commands:

Usage:

`java -jar WgLink.jar TenSQR TenSQR.config`

WgLink will generate the final haplotype results at the “output” folder: Final.Haps .

If users want use abayesQR as local region haplotype calculation tool. Users can updating absolute paths of executable (such as bwa etc) and parent folder in the aBayesQR.config file, users can run WgLink by a simple commands:

`java -jar WgLink.jar aBayesQR aBayesQR.config`

**#WgLink Parameters**<br>
##########<br>
## The name of the project.<br>
**Proj_Name** = Project_Name<br> 
## Absolute paths for fastq files.
Fastq_1_Path = /Path/to/Test.1.fastq
Fastq_2_Path = /Path/to/Test.2.fastq
## Absolute path for reference file.
Reference_Seq= /Path/to/reference.fa
## File locations: output files directory.
Output_Path = /Path/to/Output_Folder
## Reconstruction Start Position
Reconstruction_Start = 501 
## Reconstruction End Position
Reconstruction_End = 9219 
## The length of each region for divide and conquer
Region_Length = 500 
## Mapping quality cutoff
Min_Mapping_Qual = 40
## Read length cutoff
Min_Read_Length = 100 
## Maximum insert read length
Max_Insert_Length = 1000 
## Estimated sequencing error
Sequence_Err = 0.1
## MEC_Improvement_Cutoff, please refer to TenSQR user manual
MEC_Improvement_Cutoff = 0.0312 
## Initial_Population_Size, please refer to TenSQR user manual
Initial_Population_Size= 5 
## The weight of the constraint Sigma freq_i = 1, where freq_i is in the in-pool frequency for haplotype_i.
Regression_One_Vector_Weight = 50.0 
## The weight of the constraints Sigma freq_i * h_ij = MAF_j (j is the SNP index and i is the haplotype index)
Regression_Hap_MAF_Weight = 5.0  
## The weight for LD (specifically, the probability for both SNP_k and SNP_j being the alternate allele)
Regression_Hap_LD_Weight = 1.0  
## Maximum SNP mismatch ratio tolerance in region extextion
BFS_Mismatch_Tolerance_Rate = 0.1
Number_Threads= 8 
## Number of maximum selected haplotypes to generate higher level potential haplotypes for following L0L1 regression.
Maximum_Haps_R = 20 
## The minimum regularization gamma penalty for L0L1 regression.
Regression_Gamma_Min=0.0001 
## The maximum regularization gamma penalty for L0L1 regression.
Regression_Gamma_Max=0.1 
## The number of gamma values beween Regression_Gamma_Min and Regression_Gamma_Max for L0L1 regression.
Regression_n_Gamma = 10 
## The lambda penalty for L0L1 regression.
Regression_Lambda = 0.002
## The minimum frequency cutoff for haplotype. 
Min_Hap_Freq = 0.01 
## Number of maximum potential haplotypes for L0L1 regression.
Max_L0L1_Regional_Haps = 1000 
## Absolute paths for bwa.
bwa = /Path/to/bwa
## Absolute paths for Rscript.
Rscript_path= /Path/to/Rscript 
## Absolute paths for ExtractMatrix.
ExtractMatrix = /Path/to/ExtractMatrix
## Absolute paths for TenSQR.py.
TenSQR = /Path/to/TenSQR.py
## Absolute paths for PYTHON3.
PYTHON = /Path/to/PYTHON3


### Citations
If you use the TenSQR function, please cite:<br>
Ahn S, Ke Z, Vikalo H. Viral quasispecies reconstruction via tensor factorization with successive read removal[J]. Bioinformatics, 2018, 34(13): i23-i31.<br>
If you use the aBayesQR function, please cite:<br>
Ahn S, Vikalo H. aBayesQR: A bayesian method for reconstruction of viral populations characterized by low diversity[C]//International Conference on Research in Computational Molecular Biology. Springer, Cham, 2017: 353-369. <br>


### Contacts
Chen Cao, chen.cao@ucalgary.ca<br>
Quan Long, quan.long@ucalgary.ca<br>


### Copyright License (MIT Open Source)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.








