# Data

The datasets are not included in this repository. The datasets used in the
numerical experiments of the accompanying paper are listed there together
with links to their original sources.

By default, `initisic.f95` assumes that the input data file is located in
this directory (`../data/datafile.txt`). The data file and its path can be
modified in `initisic.f95` as needed.

## Data format

The input file should contain numerical values only, with one data point per
row and features given in columns. Missing values are not supported.

The number of features to be used in clustering is specified by `nft` in
`initisic.f95`. If the dataset contains a class label in the last column,
the label can be excluded from clustering by setting `nft` to the number of
columns minus one. In this case, the last column is not read by __ISIC__.