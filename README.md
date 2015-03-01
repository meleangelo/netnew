# netnew: a replication package for Mele (2015)
This is an R package for replication of the methods and analysis in __A Structural Model of Segregation in Social Networks__. 

### Installation

To install the package from the Github repository, you need to have `devtools` already installed.   
If you don't have the `devtools`, install it using the following command

```r
install.packages("devtools")
```

Then open `devtools` and use the `install_github` command to download and install `netnew`:
```r
library(devtools)
install_github("meleangelo/netnew")
```

This should install the package. This works well in Windows, using RStudio.   
To use the package just open it as an R library
```r
library(netnew)
```

Now you should be able to run the codes.


### Simulating a network model

We want to simulate a model with $n=100$ players and only direct links and indirect
links utility, starting the simulations at the empty network. We want to use a local chain and get a sapmle of $R=10000$ networks obtained sampling every 15 iterations of the sampler.


```r
library(netnew)
set.seed(1977)
n <- 100 # size of network
a <- c(-3,0,3/n) # parameters for simulations

g <- matrix(data = rbinom(n*n,1,0), nrow = n, ncol = n)
diag(g)<-0
net <- simulatemodel(model = "3params", g, a, skip = 15, samplesize = 10000)
```

The function `simulatemodel` is used to simulate some specific model. `3params` has 3 parameters: direct utility, reciprocity, and indirect links. Here we set the parameter for reciprocity to zero. The result is the sample of simulated network statistics 

```r
head(net)
  links mutual indirect
1     0      0        0
2     1      0        0
3     1      0        0
4     2      0        0
5     4      0        0
6     4      0        0
tail(net)
      links mutual indirect
9995    722     26     5300
9996    721     26     5282
9997    720     26     5266
9998    720     26     5261
9999    718     26     5227
10000   719     26     5244
```


### General method for simulating models
TBA

### Estimating a model using the exchange algorithm
TBA

### Simulating a model, getting a network adjacency matrix as output
TBA