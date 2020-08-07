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








