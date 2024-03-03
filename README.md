# GOF (Goodness Of Fit) - Feature selection for clustering scRNA-seq data

This repository contains a R package for selecting biological informative features for clustering single cell RNA-seq data. In [our study](), we developed a novel univariate distribution-oriented suite of feature selection methods, called GOF, for clustering scRNA-seq data.

The main idea of GOF is to select features based on the goodness of fit of raw UMI count data to a mixture of Negative Binomial (NB) distributions, termed “Average Negative Binomial” (ANB). The ANB distribution models generic variation such as cell library effects, so departures from ANB indicate important further structure such as cell types. We develop four variants of GOF in terms of how the goodness of fit is quantified. 


### Please cite:


![GOF framework](https://github.com/siyao-liu/GOF/blob/main/figure.png)


## Installation

```{r}
library(GOF)
```


## Example

For this vignette, we use a 3 cell line mixture dataset published from Dong et al. 2019 (https://pubmed.ncbi.nlm.nih.gov/31925417/) to demonstrate the framework of GOF and how different variants of GOF . The 3 cell line mixture dataset contains ~2,600 cells and is included in the GOF package.

```{r}
data(p3cl)
```

### Step 1: Fit ANB model. 

```{r}
myfit <- apply(p3cl, 1, function(x) { FitDist(vdata=x, family="average negative binomial") })
```

### Step 2: Quantify the goodness of fit.

#### PP method

```{r}
pp.list <- RunGOF(counts=p3cl, 
                  countsFit=myfit, 
                  method="PP", top.n=2000)
head(pp.list$PP.abc)
head(pp.list$topPP)
```

#### QQ method

```{r}
qq.list <- RunGOF(counts=p3cl, 
                  countsFit=myfit, 
                  method="QQ", top.n=2000)
length(qq.list)
names(qq.list)
```

#### 1-Wasserstein distance adjusted by gene mean

```{r}
wdist.mn.list <- RunGOF(counts=p3cl, 
                  countsFit=myfit, 
                  method="Wdist.mean", top.n=2000)
head(wdist.mn.list$Wdist.mn)
head(wdist.mn.list$topWdist.mn)
```

#### 1-Wasserstein distance adjusted by gene median

```{r}
wdist.med.list <- RunGOF(counts=p3cl, 
                        countsFit=myfit, 
                        method="Wdist.med", top.n=2000)
head(wdist.med.list$Wdist.med)
head(wdist.med.list$topWdist.med)
```


### Diagnostic Plot

We will take gene as an example and show the diagnostic plots from the 4 GOF methods.

#### P-P plot

```{r}
g <- "COL8A1"
fit.data <- myfit[[g]]
ppPlot <- ppplot(P=cumsum(fit.data$vpi), 
                 Q=cumsum(fit.data$vqi), 
                 prob1="Empirical probabilities", 
                 prob2="Theoretical probabilities")
ppPlot
```


#### Q-Q plot

```{r}
qqPlot <- qqplot_small_test(P=fit.data$Data, 
                 Q=fit.data$FittedData, 
                 sample1="Sample quantiles", 
                 sample2="Theoretical quantiles")
qqPlot
```

#### 1-Wasserstein distance adjusted by mean diagnostic plot

```{r}
wasser <- calcDiscWasser1(vdata=fit.data$Data, vpi=fit.data$vpi, vqi=fit.data$vqi, nmax=fit.data$nmax, method=c("mean"))
wasserPlot <- wasser1_plot(vigrid=fit.data$vigrid, 
                              vpi=fit.data$vpi, 
                              vqi=fit.data$vqi, 
                              vDj=wasser$vDj, 
                              Wdist=wasser$Wdist.adj, 
                              nmax=fit.data$nmax, 
                              GeneName=g)
wasserPlot
```

#### 1-Wasserstein distance adjusted by median diagnostic plot

```{r}
wasser2 <- calcDiscWasser1(vdata=fit.data$Data, vpi=fit.data$vpi, vqi=fit.data$vqi, nmax=fit.data$nmax, method=c("median"))
wasserPlot2 <- wasser1_plot(vigrid=fit.data$vigrid, 
                              vpi=fit.data$vpi, 
                              vqi=fit.data$vqi, 
                              vDj=wasser2$vDj, 
                              Wdist=wasser2$Wdist.adj, 
                              nmax=fit.data$nmax, 
                              GeneName=g)
wasserPlot2
```


## License

This software is licensed under GPL-2 License.


## Contact

If you have any questions, please contact: Siyao Liu (Siyao_Liu@med.unc.edu)

