# ISIC - Incremental Clustering Software Based on Stochastic Inexact Nonsmooth Optimization

**ISIC** is a nonsmooth optimization-based clustering algorithm for solving the minimum sum-of-squares clustering (MSSC) problem on very large and big datasets. It employs an incremental clustering framework in which clusters are added one at a time. At each incremental step, the *greedy starting-point algorithm* (GSPA) constructs a promising starting center for the new cluster, and the *stochastic inexact limited-memory bundle method* (ISLMB) is used as the optimization engine for solving the resulting clustering problem.

Due to its incremental structure, **ISIC** obtains solutions not only to the *k*-partition problem but also to all intermediate *l*-partition problems, where *l*=1,...,*k*-1.

**ISIC** builds on the incremental framework of [**bigClust**](https://github.com/napsu/bigClust), but modifies it substantially. In particular, GSPA directly constructs a starting center for each new cluster, eliminating the need for the auxiliary clustering problems and costly multi-start procedure used in earlier incremental approaches. In addition, **ISIC** uses ISLMB, a stochastic variant of [**InexactLMBM**](https://github.com/napsu/InexactLMBM), as its optimization engine. By explicitly accommodating inexact function and subgradient information, iSLMB is well suited to the inaccuracies introduced by stochastic computations.

## Files included
* `isic.f95`          
  - Mainprogram for clustering software.
* `initisic.f95`  
  - Initialization of clustering parameters and InexactLMBM. Includes modules:
    + `initclust` - Initialization of parameters for clustering.
    + `initislmb` - Initialization of ISLMB.
* `clusteringmod.f95`     
  - Subroutines for clustering software.
* `islmb.f95`              
  - ISLMB - Stochastic inexact limited memory bundle method.
* `objfun.f95`            
  - Computation of the cluster function and subgradients values.
* `subpro.f95`            
  - Subprograms for ISLMB.
* `parameters.f95`        
  - Parameters. Inludes modules:
    + `r_precision` - Precision for reals,
    + `param` - Parameters,
    + `exe_time` - Execution time.

* Makefile              
  - `Makefile`: requires a Fortran compiler (gfortran) to be installed.

* Auxiliary Python scripts
  - `run_all_datasets.py`
    + Runs __ISIC__ automatically for all datasets.
  - `run_all_seeds.py`
    + Runs all datasets with multiple random seeds and summarizes the results.
  - `run_all_seeds_subsample.py`
    + As `run_all_seeds.py`, but preserves the subsampling definition in `initisic.f95`.

## Installation and usage

To use the code:

1) Modify `initisic.f95` as needed. At minimum, select the dataset and specify the number of data points, the number of features, and the maximum number of clusters, "_nclust_".
2) Run `Makefile` by typing `make`. The `Makefile` uses `gfortran` as default.
3) Finally, run __ISIC__ by typing `./isic`.

The algorithm returns a text file containing the clustering function values, Dunn and Davies-Bouldin validity indices, and elapsed CPU times for solutions with up to _nclust_ clusters.

In addition, it returns a separate text file containing the final cluster centers for the _nclust_-partition problem, as well as the solutions to all intermediate _l_-partition problems, where
 _l_ = 1,...,_nclust_-1.

 For batch experiments, the auxiliary Python scripts can be used to run __ISIC__ automatically over all datasets, multiple random seeds, or subsampled datasets. See the comments at the beginning of each script for usage details.

## References:
* __ISIC and ISLMB:__
  - N. Karmitsa, J. Lampainen, K. Joki, M.M. Mäkelä, "ISIC - Incremental Clustering Software Based on Stochastic Inexact Nonsmooth Optimization", 2026.

* __BigClust and SLMBA:__
  - N. Karmitsa, V.-P. Eronen, M.M. Mäkelä, T. Pahikkala, A. Airola, "[Stochastic limited memory bundle algorithm for clustering in big data](https://www.sciencedirect.com/science/article/pii/S0031320325003140?via%3Dihub)", _Pattern Recognition_, Vol. 165,  111654, 2025.
    
* __Clustering:__
  - A. Bagirov, N. Karmitsa, S Taheri, "[Partitional Clustering via Nonsmooth Optimization](https://link.springer.com/book/10.1007/978-3-030-37826-4)", _Springer_, 2020.

* __InexactLMBM:__
  - J. Lampainen, K. Joki, N. Karmitsa, M.M. Mäkelä, "[Inexact limited memory bundle method](https://arxiv.org/abs/2604.08067)", 2026. 

* __LMBM:__
  - N. Haarala, K. Miettinen, M.M. Mäkelä, "[Globally Convergent Limited Memory Bundle Method for Large-Scale Nonsmooth Optimization](https://link.springer.com/article/10.1007/s10107-006-0728-2)", _Mathematical Programming_, Vol. 109, No. 1, pp. 181-205, 2007.
  - M. Haarala, K. Miettinen, M.M. Mäkelä, "[New Limited Memory Bundle Method for Large-Scale Nonsmooth Optimization](https://www.tandfonline.com/doi/abs/10.1080/10556780410001689225)", _Optimization Methods and Software_, Vol. 19, No. 6, pp. 673-692, 2004.

* __Nonsmooth optimization:__
  - A. Bagirov, N. Karmitsa, M.M. Mäkelä, "[Introduction to nonsmooth optimization: theory, practice and software](https://link.springer.com/book/10.1007/978-3-319-08114-4)", _Springer_, 2014.

## Acknowledgements
This work was part of the [PRIVIDEMA](https://prividema.eu/) project (Privacy Preserving Identity Management for
Digital Wallet and Secure Data Sharing and Processing for Cyber Threat
Intelligence Data) funded by the European Union under Grant Agreement No. 101167964.
Views and opinions expressed are, however, those of the author(s) only and do
not necessarily reflect those of the European Union or the European Cybersecurity Competence Centre. Neither the European Union nor the European
Cybersecurity Competence Centre can be held responsible for them.

   
   
