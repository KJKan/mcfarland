# How to Compare Latent Factor Models and Psychometric Network Models 

This is a short tutorial on how to compare latent factor models and psychometric network models using the R package Psychonetrics.

The tutorial accompanies the paper Kan, K.J., de Jonge, H., van der Maas, H.L.J., Levine, S.Z., & Epskamp, S. (2020). How to Compare Latent Factor Models and Psychometric Network Models. A comment on McFarland. *Journal of Intelligence*.

We illustrate here
- How a network can be extracted from the data from one sample and fitted on the data of another sample. 
     - In other words, we test if a gieven network replicates in a second sample.
- How the statistics of that network can be compared to the fit statistics of (various) factor models. 

Factor models fitted are those that were considered by McFarland (2020):
- A measurement model                                           
- A (second order) g model                                                    
- A bifactor model                                                       

The data concern WAIS US data and WAIS Hungary data (used in the network analyses of Kan, van der Maas & Levine, 2019 and Schmank et al. 2019, to which McFarland, 2020 referred to).  

# Preparation

Let's clear our workspace first (run this line only if you really want that)

```{r}
# clear work space, not run
rm(list = ls())

```

Next to package Psychonetrics, we load a few more packages, needed to make a plot of the networks, for instance.

```{r}
# load required R packages
library( "psychonetrics" )
library( "qgraph"        )
library( "dplyr"         )
```

We also need our data, the US and Hungarian WAIS correlation matrices, so let's read them in.

```{r}
# US sample
load( url ( "https://github.com/KJKan/mcfarland/blob/master/WAIS_Hungary.Rdata?raw=true" ) )
# Hungarian sample
load( url( "https://github.com/KJKan/mcfarland/blob/master/WAIS_US.Rdata?raw=true" ) )
```

According to the WAIS manuals, these are the sample sizes:

```{r} 
# sample sizes
n_US      <- 1800 
n_Hungary <- 1112 
```

For our information, what variables does the WAIS assess?

```{r}
# observed variables
( yvars <- colnames( WAIS_US ) )
```

Block Design (BD)
Similarities (SI)
Digit Span (DS)
Matrix Reasoning (MR or MA)
Vocabulary (VC or VO)
Arithmetic (AR) 
Symbol Search (SS)
Visual Puzzles (VP) 
Information (IN)
Coding (CO)
Letter Number Sequencing (LN)
Figure Weights (FW)
Comprehension (CO)
Cancellation (CA)
Picture Completion (PC)

So 16 in total.

```{r}
# number of observed variables
( ny <- length( yvars ) )
```


# Build the statistical models

Now we know how our data looks like, let's start with the real job, building our statistical models.

## Theoretical model

In theory, the WAIS measures the following latent variables, Verbal ability, Perceptual Organization, Working Memory Capacity, and Cognitive Speed.

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/TheoreticalModel.jpg)

```{r}
# latent constructs to be measured (etas)
lvars <- c( 
  "P", # Perceptual
  "V", # Verbal
  "W", # Working Memory
  "S"  # Speed
)
```

So, that's 4 latent variables
```{r}
# number latent constructs to be measured
( ne <- length( lvars ) )
```

They are indexed (measured) by the 16 subtests as follows:

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

"Theory is theory", they say. This also holds for the WAIS theoretical factor model. 

In practice, some indicators turn out to indicate more than one latent variable, see for example the WAIS-IV manual, Figure 5.2. In addition to Working Memory Capacity Arithmetic also loads on the factor Verbal Ability; In addition to Speed Figure Weights loads also laods on the Perceptual factor.

The so-called WAIS-IV 'measurement model' thus differs (slightly)from the theoretical model. (It contains two additional factor loadings.)

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/MeasurementModel.jpg)

```{r}
# this model contains two additional (cross-)loadings as compared to the theoretical model
lambda_measurement          <- lambda
lambda_measurement[  6, 2 ] <- 1 # Working Memory indicator Arithmetic on the Verbal factor
lambda_measurement[ 12, 1 ] <- 1 # Speed indicator Figure Weights on the Perceptual factor
```

Now we have established how the measurement model looks like, we can ask ourselves how to define this model in Psychonetrics. 

It's quite easy, with the function lvm (which stands for 'latent variable model').

The input of lvm is 
- a data matrix; in our case we want to analyze a correlation matrix, which is covariance matrix (in standardized form)
    - covs
- the pattern of factor loadings; in LISREL notation, this matrix is termed lambda
    - lambda
- the number of observations, so the number of participants in the sample,
    - nobs
- the way we choose to identify the model, e.g. by standardizing the latent variances
    - identification

```{r}
measurementModel <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                         lambda = lambda_measurement,
                         nobs = n_Hungary,
                         identification = "variance" )
```

According to g theory the measured latent variables correlate positively because they all depend on a common source of variance, that is 'g'. In the statistical model, g is represented by the most general factor.   

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/SecondOrdergModel.jpg)

Let's build this g model (termed the 'second order g factor model'). The matrix beta defines how the latent variables in the measurement model are influenced by g.
g hitself has no indicators (i.e. measured variables that load on it).


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
```

In the alternative model, the bifactor model, measures are not unidimensional, and would indicate g directly. Note that this goes against g theory (Jensen, 1998), in which it is stated g itself is not a cognitive ability itself. 

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/BifactorModel.jpg)

The model can also be considered a means to simply decompose the variance of a variable into variance components of course.  

Anyhow, this is how the model looks like in Psychonetrics.


```{r}
lambda_bifactor    <- cbind( lambda, g = 1 )

bifactorModel    <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                         lambda = lambda_bifactor, 
                         sigma_zeta = 'empty',
                         nobs = n_Hungary,
                         identification = "variance" )
```

For the WAIS we could also come up with a model.

Let's extracted one from the US sample

This entails we first estblish the partial correlation matrix from the US sample correlation matrix

Note: this can be considered a saturated model as it is nothing more than a rediscription of the correlation matrix.

```{r}
saturatedModel_US <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                          omega = "Full",
                          nobs = n_US )
```

More interesting is of course a more sparse network model that explains this structure.

Therefore we remove partial correlations that are insignificant or are spurious. This process is called pruning.
                          
```{r}
# remove insignificant partial correlations ('edges')
prunedModel <- saturatedModel_US %>% prune( alpha = 0.01, recursive = TRUE )

# aim for further improvement of the model
finalModel  <- prunedModel %>% stepup
```

The sketelon of this matrix is to be fitted confirmatory in the Hungarian sample

```{r}
# extract the skeleton (adjacency matrix)
adjacency <- 1*( getmatrix( finalModel, "omega" ) !=0 )

# use the skeleton to estimate the parameters in the Hungarian data
nwModel <- ggm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary,
                omega = adjacency,
                nobs = n_Hungary )
```

Now we have established how all our models look like, let's fit them all (in the Hungarian sample) against the saturated model (of the Hungarian sample)

```{r}
# Hungary
saturatedModel    <- ggm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary,
                          omega = "Full",
                          nobs = n_Hungary )
```

# Fit the models, obtain their fit statistics, and compare these statistics 

By the way, this is how we run models in Psychonetrics:

```{r}
results_saturatedModel_US   <- saturatedModel_US %>% runmodel
results_nwModel_US          <- finalModel        %>% runmodel
```
And this is how we can obtain some abslute and relative fit statistics

```{r}

fit( results_saturatedModel_US )
fit( results_nwModel_US )

compare( saturated   = results_saturatedModel_US,
         network     = results_nwModel_US )
```
So let's do that for all models we aimed to fit in the Hungarian sample

```{r}

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
  
  ```        

According to standard fit criteria ()Schermelleh et al), we would conclude that the network model fits best.

# Plot the favored model, which is the network model

This is how this model looks like

```{r}
qgraph( getmatrix( nwModel, "omega" ), 
        labels = yvars,
        groups = list( Perceptual    = which( lambda[ , 1 ] == 1 ),
                       Verbal        = which( lambda[ , 2 ] == 1 ),
                       WorkingMemory = which( lambda[ , 3 ] == 1 ),
                       Speed         = which( lambda[ , 4 ] == 1 )),
        layout = "spring" )
```

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/NetworkModel.jpg)
