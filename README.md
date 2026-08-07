# ISIC - Incremental Clustering Software Based on Stochastic Inexact Nonsmooth Optimization (bigClust 2.0)

__ISIC__ is a nonsmooth optimization-based clustering algorithm for solving the minimum sum-of-squares clustering (MSSC) problem on very large-scale and big datasets. __ISIC__ consists of two components: an incremental algorithm designed to produce high-quality, near-global solutions to the clustering problem, and the _stochastic inexact limited-memory bundle method_ (iSLMB), which is used at each iteration of the incremental algorithm to solve both the clustering and auxiliary clustering problems from different starting points.

Due to its incremental structure, __ISIC__ obtains solutions not only to the _k_-partition problem but also to all intermediate _l_-partition problems where _l_=1,…,_k_-1.

__ISIC__ differs from the original [__bigClust__](https://github.com/napsu/bigClust) in two main respects. First, it uses iSLMB, a stochastic variant of [__InexactLMBM__](https://github.com/napsu/InexactLMBM), as its optimization engine. This allows the optimization method to better account for the inexactness introduced by stochastic computations. Second, __ISIC__ employs a _greedy starting-point selection procedure_ that avoids the costly computations required by the original starting-point strategy.

## Files included
* `isic.f95`          
  - Mainprogram for clustering software.
* `initisic.f95`  
  - Initialization of clustering parameters and InexactLMBM. Includes modules:
    + `initclust` - Initialization of parameters for clustering.
    + `initislmb` - Initialization of iSLMB.
* `clusteringmod.f95`     
  - Subroutines for clustering software.
* `islmb.f95`              
  - iSLMB - Stochastic inexact limited memory bundle method.
* `objfun.f95`            
  - Computation of the cluster function and subgradients values.
* `subpro.f95`            
  - Subprograms for iSLMB.
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
* __ISIC and iSLMB:__
  - N. Karmitsa, J. Lampainen, K. Joki, M.M. Mäkelä, "Incremental Clustering Software Based on Stochastic Inexact Nonsmooth Optimization", 2026.

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

   
   
