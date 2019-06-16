#########################################################################################
#########################################################################################
############################ PARAMETER FILE TEMPLATE ####################################
#########################################################################################
#########################################################################################
# This file is the template into which users can set parameter values
# Each value should be set after the assignment operator of R -- i.e. the   <-   symbol
# If uncertain how to code the value of the parameter look at the original template or read
# the hashed out comments above the parameter assignment line.

####### NOTE: Don't include commas in any string literals (words) typed as parameters

#########################################################################################
#########################################################################################
############################ FILE-PATHING AND RUN NAME ##################################
#########################################################################################
#########################################################################################
#### NOTE: All of these parameters should be input a single values as that is how they'll be used.

# This is a logical toggle (TRUE/FALSE) of whether or not you'r working on a server (which best use bash shell!)
shape_serverFarm <- FALSE
# SHAPE generates a processed output of each simulation, this logical toggle (TRUE/FALSE)
# controls if the individual steps get removed (select TRUE to reduce drive space required)
shape_results_removeSteps <- TRUE
# This is a logical toggle (TRUE/FALSE) asking if you want SHAPE to perform replicate simulations
# by submitting independent jobs after each completes or using a for loop and running the same script
# RECOMMENDED setting is FALSE, use true only if there is a maximum run-time permitted for single executions
#	(eg: TRUE used if server enforced caps on wall-time of jobs)
shape_externalSelfing <- FALSE
# This is a logical toggle (TRUE/FALSE) of whether or not you want the job to produce an error and stop if
# any error message causes a particular replicate of the same fitness landscape to have failed.
shape_toggle_forceCompletion <- FALSE

# This is the root directory for this "experiment" within which all files input/output will be found
### NOTE: Each experiment should have their own unique workDir.  BEWARE.
shape_workDir <- "~/defaultSHAPE/"

# This is the unique string literal name you want applied to jobs related to this experiment
### WARNING: DO NOT INCLUDE any underscore "_" characters, and try to avoid special symbols
shape_save_batchBase <- "yourJob"

# This is a starting value assigned to parameter combinations, should generally be 1, but could be changed
# for reasons such as you are manually combining built sets so that parameters entered are not fully factorial across combinations.
shape_save_batchJob <- 1

# These parameters control the number of replicates you want performed, where:
# "maxReplicates" is for the number of replicates which re-use a unique fitness landscape (for measuring repeatability of dynamics)
# "uniqueReplicates" is how many uniquely seeded fitness landscapes you want created (for statistical power)
# THEREFORE: total replicates for a set of parameters will be: shape_maxReplicates * uniqueReplicates
shape_maxReplicates <- 100
uniqueReplicates <- 3

### NOTE: All the parameters in the below sections can be given more than one value, separated by a comma (and any number of spaces)
###			For example:    someParameter <- 1 , 2 , 7 ,     21         ,   9
###			Will result in your experiment having different jobs each using one of the values {1,2,7,21,9}
### NOTE: Each set of values will form part of a fully factorial experimental design, DON'T DO TOO MUCH!

#########################################################################################
#########################################################################################
############################# EXPERIMENTAL CONDITIONS ###################################
#########################################################################################
#########################################################################################

# This is a string literal for the type of disturbance events that are simulated and must be one of:
# "bottleneck"   OR   "random"   ; where the first simulates fixed size disturbance and the other simulates stochastic events.
shape_const_distType <- "bottleneck"
# This is a values for the expected disturbance factor applied to population reduction.  The first should be a value >= 1 and
# for shape_const_distType == "bottlenecks" this is the dilution factor $D$, for "random" this is the mean of a normal distribution from which $D$ is drawn.
# the second term "random" only has meaning if shape_const_distType is "random", in which case it is the shape parameter of the distribution used to draw $D$
shape_init_distPars <- c(factor = 100, random = 1)
# This value controls the number of generations between disturbance events, it can be an integer value > 0
# OR set to NULL to allow SHAPE to calculate this value based on the disturbanceSize and growth rate "r" (set below)
shape_const_growthGenerations <- NULL


#########################################################################################
#########################################################################################
############################# EVOLUTIONARY PARAMETERS  ##################################
#########################################################################################
#########################################################################################

# How many interations of the death-birth-mutation function calls do you want simulated
### (like a generation, but not necessarily depending on your parameters)
shape_numGenerations <- 2000
# This should be an integer value for the length of the simulated genomes, otherwise the number of sites at which mutations can occur.
shape_genomeLength <- 100
# This should be an integer value that represents the target number of individuals.  That target has different meaning
# depending on the parameters of your experiment but most often will relate to the carrying capacity of the environment
shape_const_focal_popValue <- 500000
# This should be a value between 0 and 1 which represents the per generation probability of an individual giving birth.
shape_const_birthProb <- 1
# This is the intrinsic growth rate, "r", which expresses the number of offspring an average individual giving birth will produce
shape_const_growthRate <- 2
# This should be a value between 0 and 1 which represents the maximum per generation probability of an individual dying
shape_const_deathProb <- 0.1
# This is a logical toggle (TRUE/FALSE) of whether or not the probability of deaths is affected by population density
shape_death_byDensity <- TRUE
# This is a scaling factor, used as an exponential value, that affects the shape of a density dependent death rate
shape_death_densityCorrelation <- 4
# This is the value at which for a density dependent death rate, the probability is maximised.
### NOTE: This can be set as a numeric value OR (RECOMMENDED) you can leave it as NULL in which case SHAPE uses - shape_const_focal_popValue
shape_death_densityCap <- NULL
# This is a logical toggle (TRUE/FALSE) as to whether or not the number of deaths which occur are replaced by a number of individuals
# calculated from a nested call to the birht function with a "Constant" growth model.
### Setting TRUE means that realised growth will follow the trend of the growth model and parameters, and that death simply exists as a
###		means of stochastic loss for genotype.
### Setting FALSE means that the realised growth will follow the trend of the model, but not necessarily the pace since deaths affect net gain.
shape_scaleGrowth_byDeaths <- TRUE
# This should be a value between 0 and 1 which represents the per generation probability of an individual mutation event occuring.
shape_const_mutProb <- 1e-5, 1e-6
# This is a logical toggle (TRUE/FALSE) of whether mutations occur only from birth events or from all living individuals (see Desai 2007)
shape_muts_onlyBirths <- FALSE
# This is a logical toggle (TRUE/FALSE) to control if mutations can cause loss of mutation, otherwise known as a revertant (ie: state 1 -> 0)
shape_allow_backMutations <- TRUE
# This should be a string literal for the growth model used to calculate births.  Currently it can be one of:
# Poisson, Constant, Exponential, Logistic
shape_const_growthForm <- "logistic"
# This should be a positive value, nominally between 0 and 1, representing the proportion of a standard biological generation
# which is to be simulated in a single discrete time step of a SHAPE run
shape_size_timeStep <- 1



###################################################################################################################################
######################################## TAKE CARE SETTING THESE FITNESS PARAMETERS ###############################################
###################################################################################################################################
# The following five parameters control the distribution of random effects used when calculating fitness of genotypes.
### NOTE: The parameters values will not be applied in a factorial way to one another, which reults
###			in the values being used as a set when building factorial parameter combinations with the simModel parameters.
###			if you do not set all combinations manually the script will try to recycle parameters.  This is fine under many circumstances but not all.
###			eg:  xsimModel <- "RMF"
###					xconst_ancestFitness <- 0, 1
###					xconstDist <- "exp", "Normal"
###					xconst_distParameters <- c(10), c(0,1)
###					xconst_distAsS <- TRUE,FALSE
###			Results in:  two parameter combinations of:    RMF 0 exp(10) TRUE    ;     RMF 1 Normal(0,1) FALSE

# This should be a string literal which defines the fitness landscape model to be implemented in the simulations
# Current models include: HoC, Additive, NK, or RMF
shape_simModel <- "RMF" , "NK" , "Additive", "HoC", "RMF" , "NK" , "Additive", "HoC"
# This should be a numeric value for what the ancestral genotpe fitness value should initially be set to (For RMF this will be recalculated....)
shape_const_ancestFitness <- 1, 0, 0, 1, 1, 0, 0, 1
# This should be the string literal which describes the probability generation function from which random fitness elements will be drawn.
# Currently accepted values: Fixed,Gamma,Uniform,Normal,Chi2,beta,exp,evd,rweibull,frechet,skewNorm
shape_constDist <- "Normal","Normal","Normal","Normal","skewNorm","skewNorm","skewNorm","skewNorm"
# The should be a concatenated vector of the parameter(s) appropriate to the distribution(s) selected.
# Each set of parameters should be in a   c()    function call, but different sets separated by comma's and spaces (like other parameters)

shape_const_distParameters <- c(0,0.2), c(0,0.2), c(0,0.2), c(0,0.2), c(0.92,0.32,-10,0.5), c(0.92,0.32,-10,0.5), c(0.92,0.32,-10,0.5), c(0.92,0.32,-10,0.5)
# This is a logical toggle (TRUE/FALSE) to define if we subtract 1 from the parameterised draws.  This is likely only needed to be TRUE if the model is "RMF"
# and you've passed a parameterised distribution which makes draws of large ( > 1) values.
shape_const_distAsS <- FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE

# RECALL: These values are used 1:1, it is recommended you supply an equal number of each parameter or a single value for any of these which can be constantly re-used.
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################



###################################################################################################################################
############################# THESE CREATE FACTORIAL COMBINATIONS ONLY WHEN APPROPRIATE MODEL IS INCLUDED #########################
###################################################################################################################################

# These next two parameters only affect a RMF fitness landscape model (See Neidhart, 2014), where the initial distance is the Hamming
# distance of the ancestral genotype from the additive effect optimum in the landscape and the theta value represents the ruggedness
# of the landscape as it controls how much wight the additive component of the landscape has.  Smaller values mean more rugged landscapes.
# The Neidhart (2014) study suggesting that theta = 0.25 was very rugged, so this may be an appropriate lower boundary.
# FYI:  theta = c/(sqrt(var(random_component))), and "c" is what is used to control the scalar of the additive component
shape_const_RMF_initiDistance <- 10
shape_const_RMF_theta <- 0.5


# This parameter is only used in the NK fitness landcape and it is the "K" of the NK model.  It controls the number of correlations between sites.
# This value should be a positive integer
shape_const_numInteractions <- 10

# This parameter is only important if the landscape model is "Fixed" and it should be a .csv file name for a table, saved in the workDir,
# which contains two columns one of which is named "binaryString" and contains the genotype binary strings (eg:   000,    010,   100,     101, etc....)
# while the other should be named "fitness" and contains numeric values to be used as the fitness of genotypes.
shape_const_fixedFrame <- "someFile_with_twoColumns.csv"

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################




#########################################################################################
#########################################################################################
######################### HOUSE KEEPING - OPTIONAL PARAMETERS  ##########################
#########################################################################################
#########################################################################################

####### These can be changed, but they are the least likely parameters you'd ever need to change.  Your choice. ##########

# This is a logical toggle (TRUE/FALSE) to control if SHAPE takes the deterministic values calculated for number of births
# and gets a stochastic value by making draws from a Poisson distribution.  RECOMMENDED: TRUE
shape_includeDrift <- TRUE
# This is a logical toggle (TRUE/FALSE) as to whether or not you want only whole individuals to be tracked,
### NOTE: setting as TRUE will require that decimal values be rounded (handled stochastically to not penalise small populations)
shape_track_asWhole <- FALSE
# This value is used to establish if a genotype has a number of individuals large enough to be considered as established and so worthy of
# tracking during post analysis.  Making this value too low will be computationally costly.  This value can be an integer relating to the exact
# number of individuals required, a value between 0 and 1 for the proportion for the total number of individuals required, or the string literal
# "Desai" to use the equation he outlines in Desai et al. 2007 - which will only ever track beneficial mutants and is based on their selective coefficient.
shape_const_estProp <- 250
# This is a value which when divided by the mutationRate (set above) gives the number of individuals of a genotype that must exist before we add it
# to a list of genotypes for which the mutational space is tracked during a run for quick reference.  This value does not affect output, only computational
# efficiency for which intermediate values are optimal
# RECOMMENDED: what you enter divided by the mutation rate and the targetNumber [ie: # / (mutationRate * targetNumber) ]  should be around 0.001
shape_const_hoodThresh <- 0.01

