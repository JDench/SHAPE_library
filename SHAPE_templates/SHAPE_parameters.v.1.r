# This is the actual file that can be directly changed by users who want to directly write new parameters for a run.
# This is the file that will contain the parameters used in common between SHAPE scripts - EXCEPTION filepathing can be independently set for the plotting
# Alternatively, and recommended, is that the SHAPE_farmerParms.v.#.r file be used to set parameters in ranges and then
# The SHAPE_farmer.v.#.r be run to build all scripts and files for a run.
### NOTE: There are parameters that can be set in this file that cannot be set elsewhere
###			BUT that's because it's not expect anyone should ever need to unless they're re-coding SHAPE.
################################################################################################
####################################### PARAMETER VALUES #######################################
################################################################################################


###################################################################################################
################################# BEGIN OF CONSTANT DEFINITIONS ###################################
###################################################################################################
# This section includes those objects/variables that the user should predefine before running the script
shape_serverFarm <- FALSE
# This is a secondary logical toggle as to whether or not this script will delete it's Steps file after pre-processing
shape_results_removeSteps <- TRUE
# This is a logical toggle as to whether the script replicates by external job creation or internally with for loops.
shape_externalSelfing <- FALSE

####### FILE AND SYSTEM LOCATIONS #######
shape_workDir <- "A_Folder/Some_subDir/"

# This is a string that can be used to uniquely identify the batch of jobs, something must be defined for funcBase,
# nothing is defined for the rep as those are within a job's directory.  The jobID and setID values can be left as NULL
# and if so they will be later interpreted as such to indicated nothing was used in naming the file
#### NOTE: DO NOT INCLUDE THE shape_sepString STRING PATTERN ANYWHERE IN THE FOLLOWING THREE OBJECTS!
shape_save_batchBase <- "yourJob"
shape_save_batchSet <- 1
shape_save_batchJob <- 1
# This overwrites the default batchString and outDirectory values based on the parameters defined
shape_save_batchString <- name_batchString(funcBase = shape_save_batchBase,
                                           func_setID = shape_save_batchSet,
                                           func_jobID = shape_save_batchJob,
                                           func_sepString = getOption("shape_sepString"))
# We calculate what the outDir should be
shape_outDir <- paste(shape_workDir,shape_save_batchString,"/",sep="")

# Now, if we're doing processing on a remote server with SLURM submissions, then we may have been
# passed an shape_outDir argument that is meant for the compute node location, in which case we'll need to save
# a new object for the final repository.  shape_outDir mya be over-written by the outside arguments.
shape_finalDir <- shape_outDir

# I will consider a run to have an active replicate of at least 1 to start a call, this
# may be modified by arguments passed through the outside and processed below.
shape_thisRep <- 1
# Now inherit arguments which have been supplied at the command line, this is before the general pre-definitions thus it requires that I have hashed out their call
# in the main body of this script.
### NOTE: the command line must be called as R CMD BATCH '--args <arg1> <arg2> <arg3> ...' <filepath>.r where each each argument
# must not have space within, as they are space separated.
#### FOR THIS SCRIPT: I assume a shape_thisRep argument is being passed.
shape_outsideArgs <- commandArgs(trailingOnly = TRUE)
# This checks if any arguments were passed, if not we warn that script defaults are used,
if(length(shape_outsideArgs)==0){
     replicate(5,print("No arguments supplied, if you called SHPAE via a batch call be certain this was intentional."))
} else {
     # If arguments are passed, then they are set as per the exact argument passed along through batch call
     for(i in 1:length(shape_outsideArgs)){
          # This evaluates the arguments being passed, which means as long as they were assignments the values should be written as objects.
          eval(parse(text= shape_outsideArgs[[i]]))
     }
}

print(shape_outDir)

# We query if the shape_outDir has a trailing "/" or "\\" value
if(!grepl("/$",shape_outDir)){
     shape_outDir <- paste(gsub("[\\]+$","",shape_outDir),"/",sep="")
}

if(!dir.exists(shape_outDir)){ dir.create(shape_outDir,recursive = TRUE) }

# This is a logical toggle to have the simulations stop if a run does not complete
# NOTE: A run may not be expected to complete if we're simulating something such as evolutionary rescue.
shape_toggle_forceCompletion <- FALSE

# This is a logcial to define if this run is recycled, steps likely ought not to be recycled... only the landscape
# I've not set this up to automatically recycle parameters -- EXCEPT: some NK and RMF model parameters which must be recycled to be equivalent
#### NOTE: If recycling then please view comments at the bottom as the script will resubmit it's next version and there are expectations of how
####		the "shape_save_batchIndex" line is to have been written.
shape_run_isRecycling <- c("Landscape" = TRUE, "Steps" = FALSE, "Parameters"=TRUE, "Neighbourhood"=FALSE)
# We define what is the max index of our recycled calls
#### NOTE: shape_thisRep value should be anything but 1 when recycling is not desired.
shape_maxReplicates <- 5
# I introduce a shape_recycle_repStart object so that I can track what is the startingIndex value proposed by
# a recycled script, this allows me to have multiple runs working toward the same end batch of replicates
#### NOTE: The fitnessLandscape will still point to the 1st replicate and thus this must be initialised.....
shape_recycle_repStart <- 1

# This is an object to control which strains we get deep neighbourhood information for
# It should be one of "none","limited","priority","full"
# setting this higher will cost more and more in post analysis runtime.
shape_const_hoodDepth <- "limited"

################################################################################
########################## POPULATION AND SIMULATION ###########################
################################################################################
# We have a focal population value which is interpreted differently based on the growth form used
# if exponential this is a starting value and the value to which the population will be disturbed during perturbations
# for logistic values this is the carrying capacity, for constant this is the population size at all times.
shape_const_focal_popValue <- 4e+06
# The size of the genotype will define the number of sites in which mutations can occur
shape_genomeLength <- 100
# This is the probability of there being a mutation per generation of our individuals (i.e. - mu_g)
shape_const_mutProb <- 1.11e-4
# This is a logical to define if mutations occur only in those newly born, or across the population at large
shape_muts_onlyBirths <- FALSE
# This is a logical toggle which controls if we force rounding of values so that individuals
# are tracked as whole, ie integer values.
shape_track_asWhole <- FALSE
# We consider the size of the time step to review, this should be a vlue between {0,1} and represents
# the proportion of the life of an individual that passes in a step, this is because: death rate = avg(birth rate) = 1
##### NOTE: This value must follow 0 < shape_size_timeStep <= 1
shape_size_timeStep <- 1
# This is the number of generations that we want to simulate
shape_numGenerations <- 10
# The death probability is a constant and gets passed as the product of our death rate (1) and our timeStep size
shape_const_deathProb <- 0.1
# This is a logical toggle for setting if death rate should be denisty dependent
shape_death_byDensity <- TRUE
# These next two values have no meaning if the death denisty logical toggle is FALSE....
# This is a value for affecting the strength of the the death rates' density dependence (where 1 means linear)
shape_death_densityCorrelation <- 4
# This value determines what is the "carrying cappacity", which is basically the unit for scalling denisty
shape_death_densityCap <- shape_const_focal_popValue
# This is the probability of birth, it is a basal value, makes most sense to be 1, as it is passed along with fitness and shape_size_timeStep
shape_const_birthProb <- 1
# We can set the ancestral fitness, but if we're using a RMF model we should let the distance to optima set this value
# Thus this will actually be reset below....
# NOTE: It is most intuitive that this value be 1 since the initial genotypeFrame uses this value for fitness and not a computed one.
###		However, provided the shape_const_relativeFitness = TRUE, this can be any value, and acts as a baseline to which the distribution's effects are added (outside RMF)
shape_const_ancestFitness <- 1
# We define what are the states possible at each site of our genotype, this ought to be only 0,1 as for the moment I only simulate binary states.
### NOTE: The program handles this for Additive, NK and RMF models only (where it has meaning!)
shape_const_siteStates <- c(0,1)
# This can take several values.  If it's numeric and less than or equal to 1 then is the proportion of total population a population must have reached,
# If this is a value greater than 1 then it is assumed to be the exact size of the lineage and I'll strip any decimal values SO 1.01 will == to any lineage that exists!
# It can also take the form of a type of calculation which has been coded in querryEstablished(), see the function for more,
# if this is anything else it's assumed to be an expression that can be evaluated in the querryEstablished() function call.
# only affects reporting and tracking not the actual growth dynamics. AND I've not got much crash proofing here, so user beware with expressions.
# Implemented functions: "Desai"
shape_const_estProp <- 1e-4
# This value is the number of individuals a lineage must have before we track its nearest neighbours in our neighbourhood reference database
# This can be any value but I've prefered to use a constant fraction divided by the mutation rate * stepSize so that the fraction controls
# what "probability" a lineage has of generating a mutant in a time step.
#### NOTE: Setting this too low will mean the reference becomes so large you may run out of disk space!  USER BEWARE!
shape_const_hoodThresh <- ceiling(0.005 / (shape_const_mutProb * shape_size_timeStep))


################################################################################
############################## GROWTH FUNCTION #################################
################################################################################
# This is the suite of parameters which affect how population growth occurs, this is a combination of disturbances to size,
# and the form of growth to carrying capacity.

# This is the disturbance type to the population, I've only coded for "bottleneck" or "random",
shape_const_distType <- "bottleneck"
# This sets the initial size of disturbance, which should be considered as the factor by which to reduce the population, to occur in the population
# as well we can include a random component, this only has meaning when the type of disturbance in "random" .. at least for now.
shape_init_distPars <- c("factor"=100,"random"=1)
# This is initialising an object which will be updated to track each time a disturbacne occurs
shape_track_distSize <- NULL
# This is the type of growth function that we'll employ, at present I've only coded for:
# constant (this will overide disturbance mechanics and shape_const_growthForm), exponential (unbounded), logistic growth (bounded) and poisson
shape_const_growthForm <- "logistic"
# We set the basal growth rate of the individuals - this represents the expected number of offspring generated by a single birth event
# this value represents (1+r) as seen in the literature of ODE's dealing with evolution or growth functions.  Another way to see this value
# is as the growth rate of a deterministic system, or the 1+growth rate of a continuous time set of equations.
shape_const_growthRate <- 2
# If a value is set, this is the specific number of generations between disturbance events.  If left as NULL then a value
# will be calculated assuming that growth is logistic and will be based on the disturbance size and growth rate
shape_const_growthGenerations <- NULL
# This is a logical toggle which when true, the number of births is scaled by the number of deaths.  This is intended for functions
# other than constant where this is assumed, since deaths are the only source of space for births.  The addition of this function is
# mostly meant for situations where prob_d > 0 and we're trying to replciate a scenario where the amount of growth is not affected by
# deaths.  Basically implemented to help replicate analytical results when prob_b == prob_d == 1, and since my functions calculate offpsring
# as the difference between current and sequential growth predicted by the growth funcitons, this adjustment becomes necessary.
shape_scaleGrowth_byDeaths <- TRUE
# This is a logical toggle asking if we want drift to be considered or not.  Drift will be simulated by passing the expected number of offpsring
# through a poisson distribuiton.  Otherwise these are deterministic growth functions (except for the stochastic rounding).
shape_includeDrift <- TRUE

# This called function defines the number of steps between disturbances, it can be used in the script to allow fluctuating disturbance times
shape_init_distSteps <- compute_distGrowth(func_distFactor = if(shape_const_growthForm == "logistic"){
																	shape_init_distPars
																}else{
																	0
																},
                                           func_growthType = shape_const_growthForm,
                                           func_distType = shape_const_distType,
                                           func_growthRate = shape_const_growthRate,
                                           func_popSize = shape_const_focal_popValue,
                                           func_focalSize = shape_const_focal_popValue,
                                           func_manualGenerations = if(shape_const_growthForm == "constant"){
																			shape_numGenerations
																	   } else {
																			shape_const_growthGenerations
																	   },
                                           func_stepDivs = shape_size_timeStep)
# If this is logistic, then we now update the amount of loss, since only now have we calculated what the factor is
### NOTE: unless loss is random, this could have been done earlier
if(shape_const_growthForm == "logistic"){ shape_init_distSteps["popLost"] <- round(shape_const_focal_popValue/shape_init_distSteps["factor"],0) }

###################################################################################################
######################################## FITNESS LANDSCAPE ########################################
###################################################################################################

# This is the type of landscape model to be used in the simulation, we should describe of types: HoC, Additive, NK, or RMF
# See the fitnessLandscape function for more information
shape_simModel <- "RMF"
# This is the maximum number of mutations that a single new mutant can carry, minimum of 1 for obvious reasons.
## NOTE: This scipt does not support a maxHammign distance of a single mutant being greater than 1, this is implemented for future proffing.
shape_max_numMutations <- 1 #max(1, round(shape_genomeLength * shape_const_mutProb,0))
# This determine if we want to use relative fitness for our initial genotype
# NOTE: If not using relative fitness be certain that the distribution is one centered about zero
#		and such that fitness values drawn are providing values with biologically meaningful selective coefficient (s) values.
shape_const_relativeFitness <- TRUE
# This is a logical of whether or not we allow backmutations to occur in our simulation
shape_allow_backMutations <- TRUE
# This distribution and parameters should represent realistic biological expectations congruent with the mutation supply rate settings
# the following values work well with a skewNorm distribution: c(0.98,0.375,-8,0.01) and gives ~ 2.5% values > 1, with mean ~ 0.686, but mode ~ 0.9
# I also like the values c(0.92,0.32,-10,0.5) for skewNorm, they give a fatter beneficial tail with a good peak near 1, but more total beneficial
# Currently accepted values: Fixed,Gamma,Uniform,Normal,Chi2,beta,exp,evd,rweibull,frechet,skewNorm
shape_constDist <- "exp"
shape_const_distParameters <- c(100)
# This is a logical toggle the user can define so that when calculating fitness, in the fitnessLandscape function
# draws from shape_constDist get -1 (if this is TRUE) so that the raw value is treated as relative fitness which is converter to selection coefficients.
# This is implemented to differentiate between distributions like skewNorm which can simulate the distribution of
# relative fitnesses, from use of exp() which is likely going to be used as selection coefficients.
shape_const_distAsS <- TRUE

# This is the initial Hamming distance that our ancester should be from the theoretical global phenotype
# NOTE: This value only has meaning if the model is RMF
shape_const_RMF_initiDistance <- 5
# From Neidhart 2014, theta ~ 0.25 was on the high order of ruggedness according to empirical fitness landscapes.
# Thus we calculate our "c" value (or independent weight) by using our defined theta, and the distribution from which we'll be drawing our random component
shape_const_RMF_theta <- 0.35

# This is the K value for an NK landscape, it has no meaning otherwise.
shape_const_numInteractions <- 4

# This object has no meaning unless the fitness landscape model is "Fixed", in which case this is the file name for something in the workDir
# for a .csv file which can be read as a data.frame containing two columns named "binaryString", "fitness".
shape_const_fixedFrame <- NULL

###################################################################################################
################################## END OF CONSTANT DEFINITIONS ####################################
###################################################################################################

# So, we now set each of the shape_<var> as options in the global environment
tmpOptions <- ls(,envir = globalenv())
tmpOptions <- tmpOptions[which(grepl("shape_(.)+",tmpOptions))]
if(length(tmpOptions) > 0){
        options(sapply(tmpOptions,function(x){ eval(as.name(x)) },simplify=FALSE))
        # We then clear the globalenv of the shape objects
        rm(list = ls()[which(grepl("^shape_",ls()))])
} else {
        stop("No options to be set for SHAPE, something is not right in the Parameters source file")
}
