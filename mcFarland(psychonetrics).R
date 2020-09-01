##############################################################################
#                                                                            #
# Kan, KJ, de Jonge, H, van der Maas, HLJ, Levine, SZ, & Epskamp, S. (2020). #
#                     comment on McFarland                                   #
#                                                                            #
# A re-analysis of the US and Hugarian WAIS IV correlation matrices          #
#                                                                            #
# In Psychonetrics we fitted                                                 #
#  - the WAIS-IV measurement model                                           #
#  - a 2nd order g model                                                     #
#  - a bifactor model                                                        #
#  - a networks extracted from the US standardization sample                 #
#                                                                            #
# Author: Kees-Jan Kan                                                       #
#                                                                            #
##############################################################################


# ------------- Preparation

# clear work space, not run
rm(list = ls())

# load required R packages
library( "psychonetrics" )
library( "qgraph"        )
library( "dplyr"         )


# ------------- Read in the data

# WAIS correlation matrix in the US standardization sample
load( url ( "https://github.com/KJKan/mcfarland/blob/master/WAIS_Hungary.Rdata?raw=true" ) )
# WAIS correlation matrix in the Hungarian standardization sample
load( url( "https://github.com/KJKan/mcfarland/blob/master/WAIS_US.Rdata?raw=true" ) )

# sample sizes
n_US      <- 1800 
n_Hungary <- 1112 

# observed variables
yvars <- colnames( WAIS_US )

# number of observed variables
ny <- length( yvars )



# ------------- Build the statistical models

# - Theoretical model

# latent constructs to be measured (etas)
lvars <- c( 
  "P", # Perceptual
  "V", # Verbal
  "W", # Working Memory
  "S"  # Speed
)

# number latent constructs to be measured
ne <- length( lvars )

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



# --- Measurement model (see WAIS-IV manual, Figure 5.2)

# this model contains two additional (cross-)loadings as compared to the theoretical model
lambda_measurement          <- lambda
lambda_measurement[  6, 2 ] <- 1 # Working Memory indicator Arithmetic on the Verbal factor
lambda_measurement[ 12, 1 ] <- 1 # Speed indicator Figure Weights on the Perceptual factor

measurementModel <- lvm( covs = ( n_Hungary - 1 )/n_Hungary*WAIS_Hungary, 
                         lambda = lambda_measurement,
                         nobs = n_Hungary,
                         identification = "variance" )


# --- Second order g factor model

# this model explains the covariance structure among the measured constructs
# by adding one more latent variable, g
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

