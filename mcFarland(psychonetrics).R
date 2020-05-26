##############################################################################
#                                                                            #
# Kan, KJ, de Jonge, H, van der Maas, HLJ, Levine, SZ, & Epskamp, S. (2020). #
#                     comment on McFarland                                   #
#                                                                            #
# A re-analysis of the US and Hugarian WAIS IV correlation matrices          #
#                                                                            #
# In Psychonetrics                                                           #
#  - an explicit saturated model is fitted                                   #
#  - the WAIS-IV measurement model is fitted                                 #
#  - a 2nd order g model is fitted                                           #
#  - McFarland's bifactor model is fitted                                    #
#  - networks extracted by qgraph are being cross-validated                  #
#                                                                            #
# Author: Kees-Jan Kan                                                       #
#                                                                            #
##############################################################################

# clear work space, not run
rm(list = ls())

# install development version of psychonetrics
# devtools::install_github("sachaepskamp/psychonetrics")

# load required R packages
library( "psychonetrics" )
library( "qgraph"        )
library( "dplyr"         )
library( "lavaan"        )

# read in data
con <- url( "https://github.com/KJKan/mcfarland/blob/master/WAIS_Hungary.Rdata?raw=true" )
load( con )
con <- url( "https://github.com/KJKan/mcfarland/blob/master/WAIS_US.Rdata?raw=true" )
load( con )

# sample sizes
n_US      <- 1800 
n_Hungary <- 1112 

# observed variables
yvars <- colnames( WAIS_US )

# number of observed variables
ny <- length( yvars )

# ------------- Build the factor models

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
                      0, 0, 1, 0, # AR # 0, 1, 1, 0
                      0, 0, 0, 1, # SS
                      1, 0, 0, 0, # VP
                      0, 1, 0, 0, # IN
                      0, 0, 0, 1, # CD
                      0, 0, 1, 0, # LN
                      0, 0, 1, 0, # FW # 1, 0, 1, 0 
                      0, 1, 0, 0, # CO
                      0, 0, 0, 1, # CA
                      1, 0, 0, 0  # PC
                     ), 
                  ncol = ne, 
                  byrow = TRUE,
                  dimnames = list( yvars, lvars ) 
                  )

# build the measurement model (see WAIS-IV manual, Figure 5.2)
# the contains two additional (cross-)loadings as compared to the theoretical model
lambda_measurement          <- lambda
lambda_measurement[  6, 2 ] <- 1 # WM indicator Arithmetic on the Verbal factor
lambda_measurement[ 12, 1 ] <- 1 # Speed indicator Figure Weights on the Perceptual factor

measurementModel <- lvm( covs = ( n_US - 1 )/n_US*WAIS_US, 
                         lambda = lambda_measurement,
                         nobs = n_US,
                         identification = "variance" )

# build a second order g factor model
# this model explains the covariance structure among the measured constructs
lambda_g               <- cbind( lambda_measurement, g = 0 )
beta_g                 <- cbind( matrix( 0, ne + 1, ne + 1 ) ) 
beta_g[ 1:ne, ne + 1 ] <- 1

gModel    <- lvm( covs = ( n_US - 1 )/n_US*WAIS_US, 
                  lambda = lambda_g, 
                  beta = beta_g,
                  sigma_zeta = 'empty',
                  nobs = n_US,
                  identification = "variance" )


# build a bifactor model 
# this model describes the variance covariance among the observed variables
# by decomposing their variance into general, specific, and shared components
lambda_bifactor    <- cbind( lambda, g = 1 )

bifactorModel    <- lvm( covs = ( n_US - 1 )/n_US*WAIS_US, 
                         lambda = lambda_bifactor, 
                         sigma_zeta = 'empty',
                         nobs = n_US,
                         identification = "variance" )

# build a network model
# estimate the partial correlation matrix from data set i
saturatedModel_i <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                         omega = "Full",
                         nobs = n_US )

# remove insignificant edges
prunedModel <- saturatedModel_i %>% prune( alpha = 0.01, recursive = TRUE )

# aim for further improvement
model_stepup <- prunedModel %>% stepup

# extract the adjacency matrix
adjacency <- 1*( getmatrix( model_stepup, "omega" ) !=0 )

# use this matrix to estimate the network in another data set (j)
nwModel <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                omega = adjacency,
                nobs = n_US )

# the saturated model (defined as network model)
saturatedModel_j <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                         omega = "Full",
                         nobs = n_US )

# run all models
results_saturatedModel   <- saturatedModel_j %>% runmodel
results_measurementModel <- measurementModel %>% runmodel
results_bifactorModel    <- bifactorModel    %>% runmodel
results_gModel           <- gModel           %>% runmodel
results_nwModel          <- nwModel          %>% runmodel

# obtain their fit statistics
compare( saturated   = results_saturatedModel,
         measurement = results_measurementModel, 
         bifactor    = results_bifactorModel,
         gmodel      = results_gModel,
         network     = results_nwModel )

# plot the network model (the winner)
qgraph( getmatrix( model_stepup, "omega" ), 
        labels = yvars,
        groups = list( Perceptual    = which( lambda[ , 1 ] == 1 ),
                       Verbal        = which( lambda[ , 2 ] == 1 ),
                       WorkingMemory = which( lambda[ , 3 ] == 1 ),
                       Speed         = which( lambda[ , 4 ] == 1 )),
        layout = "spring" )


