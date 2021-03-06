# ACST (Aggregated Conditional Score Test for Rare Variant Detection)
Original article:

Sugasawa, S., Noma, H., Otani, T., Nishino, J. and Matsui, S. (2017).  An efficient and flexible test for rare variant effects.  *European Journal of Human Genetics* 25, 752-757.  https://dx.doi.org/10.1038/ejhg.2017.43

This is a tutorial on applying aggregated conditional score test (ACST) for rare variant effects.
Required input for this testing procedure is a vector of disease status (1: case, 0: control) and a genotype matrix with each low corresponding to a vector of numbers of minor alleles (0,1 or 2) in each individual.

## Example data
In this tutorial, we use a simulated dataset contained in *example.RData*.
The file consists of disease status `y` and genotype matrix `X` simulated from a case-control study based on the same method used in Section 3.2 of Sugasawa et al. (2017). 


```{r cars}
load("example.RData")
```

Note: this dataset was generated under the following assumptions;

 * 1000 cases and 1000 controls
 * Disease prevalence of 20%
 * 20 variants
 * MAF of <1% 
 * Effect sizes (log odds ratios) of all variants are 0.3 (All variants are associated with disease status) 
 * There is no correlation among genotypes


## Testing procedure
Functions used in this tutorial are implemented in *ACST.R*.
Read R code from the file.

```{r read}
source("ACST.R")
```


Calculate p-values from two types of ACST, independent ACST (ACST-I) and correlated ACST (ACST-C), based on 500 permutations. 
```{r input}
# p-value of independent ACST (ACST-I).
ACST(y,X,B=500,C=0)

# p-value of correlated ACST (ACST-C).
ACST(y,X,B=500,C=1)
```
