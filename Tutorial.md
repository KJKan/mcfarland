# Table of Contents
- [How to Compare Latent Factor Models and Psychometric Network Models](#Intro)
- [Preparation](#Preparation)
- [WAIS-IV Factor Theoretical Model of Intelligence](#TheoreticalModel)  
- [Statistical Factor models](#Building)
     - [Measurement Model](#Measurement)
     - [g Model](#gModel)
     - [Bifactor Model](#BiModel)
- [Network Model](#Network)
     - [Exploratory Network Analysis](#Extract)
- [Fit the models (confirmatively) and compare the fit statistics](#Fits)
     - [Confirmatory Factor Analysis](FitFactorModels)
     - [Confirmatory Network Analysis](FitNetworkModel)
- [Conclusion](#Conclusion)


# How to Compare Latent Factor Models and Psychometric Network Models <a name="Intro"></a>

This is a short tutorial on how to compare latent factor models and psychometric network models using the R package Psychonetrics.

The tutorial accompanies the paper Kan, K.J., de Jonge, H., van der Maas, H.L.J., Levine, S.Z., & Epskamp, S. (2020). How to Compare Latent Factor Models and Psychometric Network Models. In Memory of John Dennis McFarland. *Journal of Intelligence*.

We illustrate here hoe
-  A network can be extracted from the data from one sample and fitted on the data of another sample. 
     - In other words, we provide the means how to test if a given network *replicates*
- The fit statistics of that network can be compared witj the fit statistics of (various) factor models. 

The gactor models being fitted are those that were considered by McFarland (2020):
- A measurement model                                           
- A (second order) g model                                                    
- A bifactor model                                                       

The data concern WAIS US data and WAIS Hungary data. These were used in the network analyses of Kan, van der Maas and Levine (2019) and Schmank et al. (2019), to which McFarland (2020) referred to.  

# Preparation <a name="Preparation"></a>

Let's clear our workspace first (run this line only if you really want that; you will lose everything you had in the workspace)

```{r}
rm(list = ls())
```

Next to package Psychonetrics, we load a few more packages, needed to make a plot of the networks, for instance.

```{r}
library( "psychonetrics" )
library( "qgraph" )
library( "dplyr" )
```

We also need our data, the US and Hungarian WAIS correlation matrices, so let's read them in.

```{r}
load( url ( "https://github.com/KJKan/mcfarland/blob/master/WAIS_Hungary.Rdata?raw=true" ) )
load( url( "https://github.com/KJKan/mcfarland/blob/master/WAIS_US.Rdata?raw=true" ) )
```

According to the WAIS manuals, these are the sample sizes:

```{r} 
n_US      <- 1800 
n_Hungary <- 1112 
```

For our information, what variables does the WAIS-IV assess again?

```{r}
( yvars <- colnames( WAIS_US ) )
```

These are:
- Block Design (BD)
- Similarities (SI)
- Digit Span (DS)
- Matrix Reasoning (MR or MA)
- Vocabulary (VC or VO)
- Arithmetic (AR) 
- Symbol Search (SS)
- Visual Puzzles (VP) 
- Information (IN)
- Coding (CO)
- Letter Number Sequencing (LN)
- Figure Weights (FW)
- Comprehension (CO)
- Cancellation (CA)
- Picture Completion (PC)

So 16 in total. Let's store that information.

```{r}
ny <- length( yvars ) 
```

# WAIS-IV Factor Theoretical Model of Intelligence <a name="TheoreticalModel"></a>

In theory, the WAIS measures the following latent variables, Verbal ability, Perceptual Organization, Working Memory Capacity, and Cognitive Speed.

```{r}
# latent constructs to be measured (etas)
lvars <- c( 
  "P", # Perceptual
  "V", # Verbal
  "W", # Working Memory
  "S"  # Speed
)
```

So, that's 4 latent variables. Let's store that information again.
```{r}
ne <- length( lvars ) 
```
So graphically, the model looks like this:

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/TheoreticalModel.jpg)

So apparently, the latent variables are indexed (measured) as follows:

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

# Build the statistical models <a name="Building"></a>

## Measurement Model <a name="Measurement"></a>

"Theory is theory", they say. The WAIS-IV ought aims to measures 4 unidimensional latent factors, but in practice unidimensionality does not hold, at least for some measures.  
- Working Memory Capacity indicator Arithmetic also loads on the factor Verbal Ability.
- Speed indicator Figure Weights also loads on the Perceptual factor.

The WAIS-IV measurement model thus deviates somewhat from the theoretical model. 

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/MeasurementModel.jpg)

Henxe, the pattern of loadings is:

```{r}
lambda_measurement          <- lambda
lambda_measurement[  6, 2 ] <- 1 # Working Memory indicator Arithmetic on the Verbal factor
lambda_measurement[ 12, 1 ] <- 1 # Speed indicator Figure Weights on the Perceptual factor
```

Now we know how the measurement model looks like, we can ask ourselves how to define this model in Psychonetrics. 

It's quite easy: With the function lvm (which stands for 'latent variable model').

The input of lvm is 
- a data matrix; in our case we want to analyze a correlation matrix, which is covariance matrix (in standardized form)
    - argument 'covs'
- the pattern of factor loadings; in LISREL notation, this matrix is termed lambda
    - argument 'lambda'
- the number of observations, so the number of participants in the sample,
    - argument 'nobs'
- the way we choose to identify the model, e.g. by standardizing the latent variances
    - argument 'identification'

```{r}
measurementModel <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                         lambda = lambda_measurement,
                         nobs = n_Hungary,
                         identification = "variance" )
```

## g Model <a name="gModel"></a>

According to *g* theory the measured latent variables correlate positively because they all depend on a common source of variance, that is '*g*'. In the statistical model, *g* is represented by the most general factor.   

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/SecondOrdergModel.jpg)

Let's build this *g* model (which we termed the 'second order g factor model'). 
The matrix beta defines how the latent variables in the measurement model are influenced by *g*.
Note that *g* itself has no indicators (i.e. measured variables that load on it).


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

## Bifactor Model <a name="BiModel"></a>

Often, a bifactor model is considered an alternative for the second order *g* factor model. 
Note that the bifactor model implies tha:
- none of the cognitive ability measures are unidimensional
- *g* has direct indicators
     - this goes against g theory (Jensen, 1998), in which it is stated *g* itself is not a cognitive ability itself. 

![](https://raw.githubusercontent.com/KJKan/mcfarland/master/BifactorModel.jpg)

However, the model can also be considered simply as a means to decompose the variance of a variable into variance components. In that sense, there is nothing wrong with it. 

Anyhow, this is how the bifactor model looks like in Psychonetrics.


```{r}
lambda_bifactor    <- cbind( lambda, g = 1 )

bifactorModel    <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                         lambda = lambda_bifactor, 
                         sigma_zeta = 'empty',
                         nobs = n_Hungary,
                         identification = "variance" )
```

# Network Model <a name="Network"></a>

For the WAIS we could also come up with a psychometric network model, in line with the mutualism model of van der Maas et al. (2006), for example.

Let's extracted one from the US sample data.

## Extract a Network Model (Exploratory Network Analysis ) <a name="Extract"></a>

We first establish the partial correlation matrix.

```{r}
saturatedModel_US <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                          omega = "Full",
                          nobs = n_US )
```

Next we remove those partial correlations that are insignificant or are deemed spurious. This process of removing is called 'pruning'.

We take $\alpha$ $0.01$.
                          
```{r}
prunedModel <- saturatedModel_US %>% prune( alpha = 0.01, recursive = TRUE )

# aim for further improvement of the model
finalModel  <- prunedModel %>% stepup
```

The 'skeleton' or 'adjacency matrix' can be considered 'the network', which we aim to be fit (truly confirmatory) in the Hungarian sample. 
It contains which edges are present (1) and which are absent (0).

```{r}
adjacency <- 1*( getmatrix( finalModel, "omega" ) !=0 )
```

# Fit the models, obtain their fit statistics, and compare these statistics <a name="Fits"></a>

Now we have established how all our models look like, let's fit them (in a different sample than we extracted the network from).

## Saturated model

Or wait, we can first establish a saturated model (explicitly), in which all other models are nested.

We can parameterisize it as the partial correlation matrix. This shows factor models are actually nested within network models. Paths are constrained by the mediating effects of the latent variables. 

```{r}
saturatedModel    <- ggm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary,
                          omega = "Full",
                          nobs = n_Hungary )
```

In psychonetrics, this is how we would 'run' a model, the saturated model in this case:

```{r}
results_saturatedModel   <- saturatedModel   %>% runmodel
```

## Factor models

Let's run all factor models first.

```{r}
results_measurementModel <- measurementModel %>% runmodel
results_bifactorModel    <- bifactorModel    %>% runmodel
results_gModel           <- gModel           %>% runmodel
```

To obtain their fit statistics one can use the function fit() 

```{r}
fit( results_measurementModel )
fit( results_bifactorModel )
fit( results_gModel )
```
And (or) function compare()

```{r}
compare( saturated   = results_saturatedModel,
         measurement = results_measurementModel )
compare( saturated   = results_saturatedModel,
         bifactor    = results_bifactorModel )
compare( saturated   = results_saturatedModel,
         gmodel      = results_gModel )

```

## Network model

Now let's do the same for the network model. (Recall that we stored this in the object 'adjacency'. 

```{r}
nwModel <- ggm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary,
                omega = adjacency,
                nobs = n_Hungary )

results_nwModel <- nwModel %>% runmodel

fit( results_gModel )

compare( saturated   = results_saturatedModel,
         network     = results_nwModel )
```        
# Conclusion <a name="Conclusion"></a>

According to standard fit criteria (Schermelleh et al), we would conclude that the network model fits best:
- The RMSEA, TLI, and CFI, for example, point at excellent fit absolute fit. 
- In addition, it has the lowest AIC and BIC
     - The confidence intervals of the RSMSEA of the network model and *g* theoretical model do not overlap
     - So if we had to choose between the network interpretation of general intelligence or g theory, we would favor the network interpretation

# Plot the favored model <a name="Plot"></a>

Let's plot our favored model!

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
