# How to Compare Latent Factor Models and Psychometric Network Models 

This is a short tutorial on how to compare latent factor models and psychomentric network models using the R package Psychonetrics.

The tutorial accompanies the paper Kan, KJ, de Jonge, H, van der Maas, HLJ, Levine, SZ, & Epskamp, S. (2020). How to Compare Latent Factor Models and Psychomentric Network Models. A comment on McFarland                                   #

We illustrate how a network can be extracted from the one sample and fitted on another sample. In other words we illustrate how to test if a network replicates.

The data concern WAIS US data and WAIS Hungary data.  

Next we illustrate how the statistics of that network can be compared to the fit statistics of (various) factor models. 

These include:
- a measurement model                                           
- a 2nd order g model                                                    
- a bifactor model                                                       
- a networks extracted from the US standardization sample                

# Preparation

Let's clear our workspace first (only run if you really want that)

```{r}
# clear work space, not run
rm(list = ls())

```

Next to package Psychonetrics, we need to load a few more packages, e.g. to make a plot of the networks.

```{r}
# load required R packages
library( "psychonetrics" )
library( "qgraph"        )
library( "dplyr"         )
```

Of course we also need our data, the US and Hungarian WAIS correlation matrices

```{r}
# WAIS correlation matrix in the US standardization sample
load( url ( "https://github.com/KJKan/mcfarland/blob/master/WAIS_Hungary.Rdata?raw=true" ) )
# WAIS correlation matrix in the Hungarian standardization sample
load( url( "https://github.com/KJKan/mcfarland/blob/master/WAIS_US.Rdata?raw=true" ) )
```

According to the manuals, these are the sample sizes:

```{r} 
# sample sizes
n_US      <- 1800 
n_Hungary <- 1112 
```

For our information, what variables were actually measured?

```{r}
# observed variables
( yvars <- colnames( WAIS_US ) )
```

Ok, so how many variables are that?

```{r}
# number of observed variables
ny <- length( yvars )
```


# Build the statistical models

Now we 'know' our data, let's start with the rela job, building the network and latent variable models

## Theoretical model

Theoretically speaking the WAIS measures the following latent variables, Verbal ability, Perceptual Organization, Working Memory Capacity, and Cognitive Speed.

```{r}
# latent constructs to be measured (etas)
lvars <- c( 
  "P", # Perceptual
  "V", # Verbal
  "W", # Working Memory
  "S"  # Speed
)
```

So how many latent variables is the WAIS supposed to measure? (4)
```{r}
# number latent constructs to be measured
( ne <- length( lvars ) )
```

These 4 latent variables are indexed (measured) by the 16 subtests as follows:

```{r}
# theoretical pattern of factor loadings
lambda <- matrix( c (
  #P  V  W  S   
  1, 0, 0, 0, # BD
  0, 1, 0, 0, # SI
  0, 0, 1, 0, # DS
  1, 0, 0, 0, # MR
  0, 1, 0, 0, # VC
  0, 0, 1, 0, # AR 
  0, 0, 0, 1, # SS
  1, 0, 0, 0, # VP
  0, 1, 0, 0, # IN
  0, 0, 0, 1, # CD
  0, 0, 1, 0, # LN
  0, 0, 1, 0, # FW  
  0, 1, 0, 0, # CO
  0, 0, 0, 1, # CA
  1, 0, 0, 0  # PC
), 
ncol = ne, 
byrow = TRUE,
dimnames = list( yvars, lvars ) 
)
```

Theory is theory, they say. In practice, some indicators indicate more than one latent variable, see for example the WAIS-IV manual, Figure 5.2.

The WAIS-IV measurement model is therefore actually slightly different from the theoretical model. It contains two additional factor loadings.

```{r}
# this model contains two additional (cross-)loadings as compared to the theoretical model
lambda_measurement          <- lambda
lambda_measurement[  6, 2 ] <- 1 # Working Memory indicator Arithmetic on the Verbal factor
lambda_measurement[ 12, 1 ] <- 1 # Speed indicator Figure Weights on the Perceptual factor
```

Now we know how the measurement model looks like, we can define it in Psychonetrics using the function lvm (which stands for 'latent variable model').

The input of lvm is 
- a data matrix, in our case we want to analyze a correlation matrix, which is covariance matrix (in standardized form); covs
- the pattern of factor loadings; in LISREL notation this matrix is termed lambda
- the number of observations - obs - which is the numbe rof participants in the sample
- how to identify the model, e.g. by standardizing the latent variances

```{r}
measurementModel <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                         lambda = lambda_measurement,
                         nobs = n_Hungary,
                         identification = "variance" )
```

According to g theory these latent variables correlate positively because they all depend on a common source of variance, g. In the statistical model, the representation of g is a general factor.   

Let's build this model - the second order g factor model


```{r}
# this model explains the covariance structure among the measured constructs
# by adding one more latent variable, g
```

g itself has no direct indicators

```{r}
lambda_g               <- cbind( lambda_measurement, g = 0 )
beta_g                 <- cbind( matrix( 0, ne + 1, ne + 1 ) ) 
beta_g[ 1:ne, ne + 1 ] <- 1

gModel    <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                  lambda = lambda_g, 
                  beta = beta_g,
                  sigma_zeta = 'empty',
                  nobs = n_Hungary,
                  identification = "variance" )


# --- Bifactor model 

# this model describes the variance covariance among the observed variables
# by decomposing their variance into general, specific, and shared components
lambda_bifactor    <- cbind( lambda, g = 1 )

bifactorModel    <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                         lambda = lambda_bifactor, 
                         sigma_zeta = 'empty',
                         nobs = n_Hungary,
                         identification = "variance" )


# --- Network model

# this model will be extracted from the US sample, using the following steps

# estimate the partial correlation matrix from the US sample correlation matrix
saturatedModel_US <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                          omega = "Full",
                          nobs = n_US )

# remove insignificant partial correlations ('edges')
prunedModel <- saturatedModel_US %>% prune( alpha = 0.01, recursive = TRUE )

# aim for further improvement of the model
finalModel  <- prunedModel %>% stepup

# extract the skeleton (adjacency matrix)
adjacency <- 1*( getmatrix( finalModel, "omega" ) !=0 )

# use the skeleton to estimate the parameters in the Hungarian data
nwModel <- ggm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary,
                omega = adjacency,
                nobs = n_Hungary )


# --- Saturated models (defined as network models, so partial correlation matrices)

# US
saturatedModel_US <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                          omega = "Full",
                          nobs = n_US )

# Hungary
saturatedModel    <- ggm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary,
                          omega = "Full",
                          nobs = n_Hungary )


# ------------- Fit statistics for the US network model

results_saturatedModel_US   <- saturatedModel_US %>% runmodel
results_nwModel_US          <- finalModel        %>% runmodel

compare( saturated   = results_saturatedModel_US,
         network     = results_nwModel_US )



# ------------- Fit the factor and network models in the Hungarian sample

results_saturatedModel   <- saturatedModel   %>% runmodel
results_measurementModel <- measurementModel %>% runmodel
results_bifactorModel    <- bifactorModel    %>% runmodel
results_gModel           <- gModel           %>% runmodel
results_nwModel          <- nwModel          %>% runmodel

# print their fit statistics
fit( results_saturatedModel )
fit( results_measurementModel )
fit( results_bifactorModel )
fit( results_gModel )
fit( results_nwModel )

compare( saturated   = results_saturatedModel,
         measurement = results_measurementModel )

compare( saturated   = results_saturatedModel,
         bifactor    = results_bifactorModel )

compare( saturated   = results_saturatedModel,
         gmodel      = results_gModel )

compare( saturated   = results_saturatedModel,
         network     = results_nwModel )


# ------------- Plot the network model 

qgraph( getmatrix( nwModel, "omega" ), 
        labels = yvars,
        groups = list( Perceptual    = which( lambda[ , 1 ] == 1 ),
                       Verbal        = which( lambda[ , 2 ] == 1 ),
                       WorkingMemory = which( lambda[ , 3 ] == 1 ),
                       Speed         = which( lambda[ , 4 ] == 1 )),
        layout = "spring" )

```
