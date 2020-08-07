## Users’ Manual of WgLink (Version 1.0)
### Preamble
The WgLink program reconstructs haplotypes.


The workflow of WgLink is briefed as follows: (A) WgLink first determines local region & tiling region haplotypes by applying tools: TenSQR or aBayesQR. 
(B) Regional haplotypes are joined together using breadth first search which could tolerate specified SNV mismatch when extending regions. C) The
haplotypes is estimated by L0L1 regularized regression to solve for the haplotype distribution calculation.

### Installation
WgLink.jar ibatteries-included JAR executable. All needed external jar packages are included in the downloadable, PoolHapX.jar. To download all necessary files, users can use the command 
`git clone https://github.com/theLongLab/PoolHapX.git`
Please notice that "git lfs" should be installed before downloading PoolHapX since there exists several very large files which can not be downloaded by "git clone" directly. Users can install "git lfs" according to https://github.com/git-lfs/git-lfs/wiki/Installation. As we used an R package L0Learn, the users have to install R and L0Learn (https://cran.r-project.org/web/packages/L0Learn/index.html). The versions of R and R package L0Learn that we have used on our platform are: version 1.2.0 for L0Learn and version 3.6.1 for R. Other versions are not tested, although they may work. Users are also expected to have java (version: 1.8) on their platform.Longranger (version: 2.2.2) should be installed if processing 10x linked-reads.
