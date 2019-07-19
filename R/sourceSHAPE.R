# This is the main source file for the SHAPE package.
# It includes all the information from my Body, Functions and Analysis
##### RECALL:  # To compile this package don't forget to start with the following
#library(devtools)
#document()
#### We can check via
# devtools::check()
## This is for building
#build()
## Installing the package
#install()
## As well writting a pdf manual
#setwd("C:/Users/Jonathan/Documents/Programming/MyScripts/SHAPE/SHAPE")
#system("R CMD Rd2pdf . --title=rSHAPE --output=./SHAPE_manual.pdf --force --no-clean --internals")

# But, I've found from reading that when building our library we actually use the import functions
# and with Roxygen2 I use the following nomenclature that gets set to NAMESPACE
#' @importFrom abind abind
#' @importFrom graphics hist
#' @importFrom sn rsn
#' @importFrom VGAM rfrechet
#' @importFrom evd rgev rrweibull
#' @importFrom stats sd rpois rbinom rgamma runif rnorm rchisq rbeta rexp rweibull weighted.mean setNames var
#' @importFrom utils read.csv write.table
#' @import RSQLite
#' @import DBI
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#'
NULL

# Here are the library dependencies
library(abind) # This allows me to bind array objects along any dimension, including creating new ones.
library(RSQLite) # This allows SQL database calls
library(DBI) # This allows database interface
library(sn) # This allows the skewwed normal distribution
library(VGAM) # This includes the Fretchet distribution to be called
library(evd) # This allows the extreme value distributions to be used.
library(foreach) # This is for parallel processing
library(parallel) # This is for parallel backend robust to OS.
library(doParallel)


###################################################################################################
#################################### BEGIN OF PREDEFINITIONS ######################################
###################################################################################################
#' These are some global reference options that SHAPE will use and I consider the defaults.
#' SHAPE parameters can be changed by calling this function and changing values OR
#' by using the accessory SHAPE_parameters script, called in the SHAPE_runBody script.
#' This second approach is considered more practical for building and running experiments.
#'
#' @param shape_allow_backMutations This is a logical toggle controlling if revertant mutants are allowed.
#' @param shape_collapseString This is a string to collapse the progenitor and number of mutants pieces of information.
#' @param shape_constDist This is a character string to control the distribution used for drawing fitness value random components.
#' @param shape_const_relativeFitness This is a logical toggle which controls if the absolute fitness values calculated should be reinterpreted as relative fitness values.
#' @param shape_const_hoodDepth shape_const_hoodDepth This is an object to control which strains we get deep neighbourhood information for
#' It should be one of \strong{"none","limited","priority","full"}
#' setting this higher will cost more and more in post analysis runtime.
#' @param shape_const_focal_popValue This is the focal population value which has different meanings based on the growth model implemented.
#' @param shape_const_mutProb This is the probability of a mutation event - occuring relative to the number of mutable events - in a standard biological generation.
#' @param shape_const_distParameters This allows a single parameter to be passed for use in the distribution of fitness fitness effects.  NOTE: you are likely going to
#' want to pass multiple values in which case simply set this value prior to a run's start but after loading the library.
#' @param shape_const_distAsS This is a logical toggle controlling if fitness landscape values calculated should be interpreted as selection coefficients rather than relative
#' fitness values.
#' @param shape_const_RMF_initiDistance This is the distance of the independent global fitness optima away from the WT genotype.  It matters for the Rough Mount Fuji landscapes.
#' @param shape_const_RMF_theta This is the Rough Mount Fuji value that controls the scalar of the independent fitness contribution.
#' @param shape_const_numInteractions This is the number of sites which interact with respect to fitness calculations in models such as the NK.
#' @param shape_const_fixedFrame This defines the fitness landscape when our model is "Fixed", it must be user defined and be explicit to all genotypes possible.
#' @param shape_const_deathProb This is the proportion of individuals having a death event in a standard biological generation.
#' @param shape_const_birthProb This is the proportion of individuals with fitness == 1 having births events in a standard biological generation.
#' @param shape_const_ancestFitness This is the fitness value of the ancestral genotype.
#' @param shape_const_estProp This is the value controlling when SHAPE considers a population to be established.
#' @param shape_const_hoodThresh This is the numeric value controlling when a population is of sufficient size for SHAPE to consider
#' it worth having the genotype's mutational neighbourhood to be stored in a convenience DB for easier access - ie: this can save
#' computational time but will cost diskspace during the run.
#' @param shape_const_distType This is the type of stochastic disturbance events to be simulated.
#' @param shape_const_growthForm This is the growth form model to be simulated
#' @param shape_const_growthRate This is the number of offspring from every division event where 1 would mean replacement,
#' 2 is normal binary fission, etc....
#' @param shape_const_growthGenerations This is an optional integer value controlling if you want a standard number of time steps between
#' each stochastic disturbance function call.  Not defining this means it will be calculated based on other paratmerts defined.
#' @param shape_db_splitTables This is a logical toggle as to whether or not fitness landscape tables - for genotypes with the same number
#' of mutations - are allowed to be split into sub-tables.
#' @param shape_death_byDensity This is the logical toggle controlling if deaths are density dependent.
#' @param shape_death_densityCorrelation This is a positive numeric controlling the rate at which density dependent deaths increase from
#' minimal to maximal effect.  Where 1 is linear, > 1 creates an exponential form of curve and values < 1 will create a root function curve.
#' @param shape_death_densityCap If deaths are density dependent this is the maximal community size for when deaths are 100\% expected.
#' @param shape_envString This is a string used for programatically creating workspace environments for rSHAPE
#' @param shape_externalSelfing This is the logical toggle controlling if replicates are to be handled as individual external calls
#' rather than through the normal internal for loop.  It has limited value and was desgined for when you work on compute nodes with
#' limited wall time.
#' @param shape_external_stopFile This is the filename for a file which is used to control self-replciation of SHAPE when selfing is external.
#' @param shape_finalDir This is the directory where file from a remote server's compute node are to be back ported regularly.  Only matters
#' under the correct conditions.
#' @param shape_genomeLength This is the length of a simulant's genome, or in other words the number of sites where mutations can occur.
#' @param shape_includeDrift This is a logical toggle as to whether or not we should add stochasticity to the growth function
#' calculations.  It is meant to simulate drift in calculations that would otherwise be deterministic.
#' @param shape_init_distPars This is the vector of initial values of the dilution factor and random component of the stochastic disturbance function.
#' It needs to be set with a number and range of values approriate to the distribution to be simulated.
#' @param shape_maxReplicates This is the number of replicates to be run.
#' @param shape_maxRows This is the integer number of rows stored in a single table of the fitness landscape DB.  Only
#' matters is tables are aplit/
#' @param shape_muts_onlyBirths This is a logical flag to control if mutants only appear as a result of birth events.
#' @param shape_nextID This is the next genotype ID to be assigned for a genotype that get's created.
#' @param shape_numGenerations This is the number of generations to be simulated in the run.
#' @param shape_objectStrings This is a named character vector which are the string prefixes used when programatically naming objects.
#' @param shape_postDir This is the filepath to the directory where post-analysis results will be stored.
#' @param shape_recycle_repStart This is the first replicate being simulated once a SHAPE call is made.
#' @param shape_results_removeSteps This is a logical flag controlling if the steps log is removed after being processed.
#' @param shape_run_isRecycling This is a named vector of four logicals which control which parts of a run is meant to be recycled between replicates.
#' @param shape_save_batchBase This is a character string for naming your experiment.
#' @param shape_save_batchSet This is an integer value for the set of this experiment associated to this job.
#' @param shape_save_batchJob This is an integer value for the batch of this experiment associated to this job.
#' @param shape_scaleGrowth_byDeaths This is a logical flag that controls if growth is scaled by deaths so that the growth
#' form follows standard expectations.
#' @param shape_sepLines This is a string character that is used in collapsing multiple elements into a single character string
#' though namely employed in the summariseExperiment function.
#' @param shape_sepString This is a string character that is used for collpasing vectors of information into a single
#' character string, and subsequently splitting that information back out.
#' @param shape_serverFarm This is a logical flag of whether or not your simulations are going to be run on a remote server
#' or other situation with compute and host nodes where you might want to handle particularities I experienced and thus accounted for.
#' @param shape_simModel This is the fitness landscape model to be simulated.
#' @param shape_size_timeStep This is the proportion of a standard biological generation to be simulated in a single time step
#' of a SHAPE run.  Values greater than 1 are not guaranteed to work as expected.  Negative numbers will cause errors.
#' @param shape_stringsAsFactors I don't like strings to be factors and so SHAPE will avoid treating them as so.
#' @param shape_string_lineDescent This is a string that will be used to collapse vectors of character strings into a single string
#' It get's used when we are tracking sequential genotypes through the line of descent
#' @param shape_string_tableNames This is a string value used as the prefix when naming table in the fitness landscape DB.
#' @param shape_thisRep This is the replicate number of the first replicate processed in the called run.
#' @param shape_tmpGenoTable This is a temporary object of a table of genotype information that is to be passed along
#' different functions of SHAPE.  It's stored as an option since it can be build within a function where it is not returned
#' as an object but then used later.  There is little value in setting this manually.
#' @param shape_tmp_selfScript This is an optionally defined filepath location for a file that will exist to signal that
#' an externally replicating SHAPE run can stop.  This only matters if selfing is external.
#' @param shape_use_sigFig This is the number of significant figures that will be kept for processed output.
#' @param shape_toggle_forceCompletion This is a logical toggle controlling if a run crashes when it is ended prior to the
#' maximum number of replicates being completed.
#' @param shape_track_asWhole This is a logical toggle controlling if population sizes must be tracked as integer values
#' @param shape_track_distSize This is a numeric, the size of a disturbance caused by stochastic events.  It is the dilution factor
#' or the divisor of the community size.  It must be > 1 or is forced to that value.
#' @param shape_workDir This is the main working directory relative to which your SHAPE experiment will be built and run.  It defaults
#' to the -- tempdir -- of R when this value is NULL, I strongly recommend
#'
#' @section Warning:
#' Please pass a directory filepath to the argument of shape_workDir, rSHAPE will create this so it needn't exist yet.
#' If you leave it as the default -- ie NULL -- whatever is created will simply be lost in the temporary folder
#' of this R sessions' workspace.
#'
#' @examples
#' # This function builds the basic parameters for a run of SHAPE and I recommend as
#' # the most convenient wayfor setting your own parameters since this function will
#' # make appropriate derived settings based on values passed.
#' # You must at least call it before using runSHAPE() or shapeExperiment().
#'
#' # You can see there are a lot of parameters for SHAPE
#' args(defineSHAPE)
#' # Here are some default values that were just loaded as options
#' sapply(c("shape_workDir","shape_save_batchJob","shape_save_batchBase", "shape_simModel"),getOption)
#' # As an exmaple we change your working directory, the ID of the job and the fitness landscape model
#' options(list("shape_workDir" = paste(tempdir(),"~/alternativeFolder/",sep=""),
#'                 "shape_save_batchJob" = 3, "shape_save_batchBase" = "non_default_Experiment",
#'                 "shape_simModel" = "NK"))
#' sapply(c("shape_workDir","shape_save_batchJob","shape_save_batchBase", "shape_simModel"),getOption)
#' # NOTE: that manually setting the options will not create a new working directory for rSHAPE,
#' # you would need to do this yourself or could simply pass these arguments through a call
#' # to defineSHAPE().
#'
#' @export
defineSHAPE <- function(shape_allow_backMutations = TRUE,
                                 shape_collapseString = "__:__",
                                 shape_constDist = "exp",
                                 shape_const_relativeFitness = TRUE,
                                 shape_const_hoodDepth = "limited",
                                 shape_const_focal_popValue = 1e5,
                                 shape_const_mutProb = 1e-3,
                                 shape_const_distParameters = 20,
                                 shape_const_distAsS = FALSE,
                                 shape_const_RMF_initiDistance = 5,
                                 shape_const_RMF_theta = 0.35,
                                 shape_const_numInteractions = 4,
                                 shape_const_fixedFrame = NULL,
                                 shape_const_birthProb = 1,
                                 shape_const_deathProb = 1,
                                 shape_const_ancestFitness = 0,
                                 shape_const_estProp = 1e-3,
                                 shape_const_hoodThresh = 1e3,
                                 shape_const_distType = "bottleneck",
                                 shape_const_growthForm = "logistic",
                                 shape_const_growthRate = 2,
                                 shape_const_growthGenerations = NULL,
                                 shape_db_splitTables = TRUE,
                                 shape_death_byDensity = TRUE,
                                 shape_death_densityCorrelation = 4,
                                 shape_death_densityCap = NULL,
                                 shape_envString = "shapeEnvir",
                                 shape_externalSelfing = FALSE,
                                 shape_external_stopFile = "someNamed.file",
                                 shape_finalDir = NULL,
                                 shape_genomeLength = 1e2,
                                 shape_includeDrift = TRUE,
                                 shape_init_distPars = c("factor"= 100,
                                                          "random"= 1),
                                 shape_maxReplicates = 30,
                                 shape_maxRows = 2.5e7,
                                 shape_muts_onlyBirths = FALSE,
                                 shape_nextID = 0,
                                 shape_numGenerations = 100,
                                 shape_objectStrings = c("popDemographics" = "popDemo",
                                                         "repeatability" = "evoRepeat"),
                                 shape_postDir = NULL,
                                 shape_recycle_repStart = 1,
                                 shape_results_removeSteps = TRUE,
                                 shape_run_isRecycling = c("Landscape" = TRUE,
                                                          "Steps" = FALSE,
                                                          "Parameters"=TRUE,
                                                          "Neighbourhood"=FALSE),
                                 shape_save_batchBase = "yourJob",
                                 shape_save_batchSet = 1,
                                 shape_save_batchJob = 1,
                                 shape_scaleGrowth_byDeaths = TRUE,
                                 shape_sepString = "_",
                                 shape_sepLines = "__and__",
                                 shape_serverFarm = FALSE,
                                 shape_simModel = "HoC",
                                 shape_size_timeStep = 1,
                                 shape_stringsAsFactors = FALSE,
                                 shape_string_lineDescent = "_->_",
                                 shape_string_tableNames = "numMutations",
                                 shape_thisRep = 1,
                                 shape_tmpGenoTable = NULL,
                                 shape_tmp_selfScript = "~/random_nullFile.txt",
                                 shape_use_sigFig = 4,
                                 shape_toggle_forceCompletion = FALSE,
                                 shape_track_asWhole = FALSE,
                                 shape_track_distSize = NULL,
                                 shape_workDir = NULL){
    # This sets the global options
    options( sapply(ls(),function(envirParms){ eval(as.name(envirParms))  },simplify=FALSE) )
    # If the user has not defined a working directory, ie: shape_workDir == NULL then we set so here
    if(is.null(getOption("shape_workDir"))){
      options("shape_workDir" = gsub('\\','/',tempdir(),fixed=TRUE))
    }
    # Also, if the workDir does not end in a trailing '/', we add one
    if(!grepl("/$",getOption("shape_workDir"))){
      options("shape_workDir" = paste(getOption("shape_workDir"),'/',sep=""))
    }

    # If the user has not define an explicit shape_death_densityCap, we set it to the focal population
    if(is.null(getOption("shape_death_densityCap"))){
     options("shape_death_densityCap" = getOption("shape_const_focal_popValue"))
    }

    # You are welcome to change these logical toggles, but they're use is meaningful in that
    # replicates of a run should be based off the same parameters and fitness landscape, hence these
    # being the standard recycling of DB file settings.

    # If all the parameter recycling values have been set to FALSE, then we're not recycling and so
    # we force replicates to end after this instance.
    if(!any(getOption("shape_run_isRecycling"))){
      options("shape_maxReplicates" = getOption("shape_thisRep"))
    }

    # These are file path objects that are relative and hence calculated based on purely other parameters.
    options("shape_save_batchString" = name_batchString(funcBase = getOption("shape_save_batchBase"),
                                                        func_setID = getOption("shape_save_batchSet"),
                                                        func_jobID = getOption("shape_save_batchJob"),
                                                        func_sepString = getOption("shape_sepString")))
    # We calculate what the outDir should be
    options("shape_outDir" = paste(getOption("shape_workDir"),
                                   getOption("shape_save_batchString"),
                                   "/",sep=""))
    # Similarly, if the post-analysis directory is not explicitly stated, we define it here
    if(is.null(getOption("shape_postDir"))){ options("shape_postDir" = paste(getOption("shape_workDir"), "postAnal/",sep="")) }
    # Now, for the outDir and postDir, if they've been defined without a trailing "/", AND/OR
    # do not exist, then we update the values and create as required.
    for(thisArgument in c("shape_outDir","shape_postDir")){
      func_tmpValue <- getOption(thisArgument)
      # We check the the lagging "/"
      if(!grepl("/$",func_tmpValue)){
        options(setNames(list(paste(gsub("[\\]+$","",func_tmpValue),"/",sep="")),
                         thisArgument))
      }
      # Now we check that the directory exists
      if(!dir.exists(getOption(thisArgument))){
        dir.create(getOption(thisArgument),recursive = TRUE)
      }
    }
    # Now, if the user has set a non-NULL finalDir, we calulate it to be the outDir here
    if(is.null(getOption("shape_finalDir"))){
      options("shape_finalDir" = getOption("shape_outDir"))
    }
    # You could change this name, but I don't see why you'd need to and so for convenience I define it here.
    options("shape_processedData_filePattern" = paste("processed_runData_from_",
                                                      getOption("shape_save_batchBase"),
                                                      sep=""))
    options("shape_processedData_fileName" = paste(getOption("shape_outDir"),
                                                   "processed_runData_from_",
                                                   getOption("shape_save_batchString"),
                                                   "_",
                                                   getOption("shape_thisRep"),
                                                   ".RData",sep=""))
    # These are more optional filenames for an experiments secondary processing
    # These are filenames for information of collected batches of data
    options("shape_procExp_filenames" = c("fileList"= paste(getOption("shape_postDir"),
                                                                        "allFiles_from_",
                                                                        getOption("shape_save_batchBase"),
                                                                        ".RData",sep=""),
                                                       "parameters"= paste(getOption("shape_postDir"),
                                                                           "jobParameters_from_",
                                                                           getOption("shape_save_batchBase"),
                                                                           ".RData",sep=""),
                                                       "popDemographics"= paste(getOption("shape_postDir"),
                                                                                "popDemographics_from_",
                                                                                getOption("shape_save_batchBase"),
                                                                                ".RData",sep=""),
                                                       "repeatability"= paste(getOption("shape_postDir"),
                                                                              "repeatabilityData_from_",
                                                                              getOption("shape_save_batchBase"),
                                                                              ".RData",sep="")))

    # These are options that are not the be changed unless you're commited to recoding some-to-all of SHAPE's body.
    options( list("shape_max_numMutations" = 1,
                    "shape_popMat_colnames" = c("numMuts","genotypeID","popSize","fitness"),
                    "shape_reportMat_colnames" = c("numMuts","genotypeID","popSize","fitness","births","deaths","mutants","progenitor"),
                    "shape_processedObjects" = c("runDemographics","info_estLines"),
                    "shape_const_siteStates" = c(0,1),
                    "shape_init_distSteps" = compute_distGrowth(func_distFactor = if(getOption("shape_const_growthForm") == "logistic"){
                                                                                    getOption("shape_init_distPars")
                                                                                  }else{
                                                                                    0
                                                                                  },
                                                                func_growthType = getOption("shape_const_growthForm"),
                                                                func_distType = getOption("shape_const_distType"),
                                                                func_growthRate = getOption("shape_const_growthRate"),
                                                                func_popSize = getOption("shape_const_focal_popValue"),
                                                                func_focalSize = getOption("shape_const_focal_popValue"),
                                                                func_manualGenerations = if(getOption("shape_const_growthForm") == "constant"){
                                                                                            getOption("shape_numGenerations")
                                                                                          } else {
                                                                                            getOption("shape_const_growthGenerations")
                                                                                          },
                                                                func_stepDivs = getOption("shape_size_timeStep")),
                    "shape_saveParameters" = list("Population"= c("shape_const_focal_popValue",
                                                               "shape_genomeLength",
                                                               "shape_const_mutProb",
                                                               "shape_muts_onlyBirths",
                                                               "shape_numGenerations",
                                                               "shape_size_timeStep",
                                                               "shape_const_deathProb",
                                                               "shape_death_byDensity",
                                                               "shape_death_densityCorrelation",
                                                               "shape_death_densityCap",
                                                               "shape_const_ancestFitness",
                                                               "shape_const_estProp"),
                                                   "Growth_Disturbance"= c("shape_const_distType",
                                                                           "shape_const_growthGenerations",
                                                                           "shape_init_distPars",
                                                                           "shape_track_distSize",
                                                                           "shape_const_growthForm",
                                                                           "shape_const_growthRate",
                                                                           "shape_scaleGrowth_byDeaths",
                                                                           "shape_includeDrift"),
                                                   "FitnessLandscape"= c("shape_simModel",
                                                                         "shape_max_numMutations",
                                                                         "shape_const_relativeFitness",
                                                                         "shape_allow_backMutations",
                                                                         "shape_const_distAsS"),
                                                   "NK_modelElements" = c("shape_const_numInteractions",
                                                                          "shape_const_NK_interactionMat",
                                                                          "shape_const_siteBystate_fitnessMat"),
                                                   "RMF_modelElements" = c("shape_const_RMF_theta",
                                                                           "shape_const_RMF_indWeight",
                                                                           "shape_const_RMF_initiDistance",
                                                                           "shape_const_RMF_globalOptima"),
                                                   "Fixed_Landscape"=c("shape_const_fixedFrame"),
                                                   "DFE"= c("shape_constDist",
                                                            "shape_const_distParameters"),
                                                   "DataManagement" = c("shape_sepString",
                                                                        "shape_collapseString",
                                                                        "shape_string_tableNames",
                                                                        "shape_db_splitTables",
                                                                        "shape_maxRows",
                                                                        "shape_thisRep",
                                                                        "shape_maxReplicates",
                                                                        "shape_save_batchBase",
                                                                        "shape_save_batchSet",
                                                                        "shape_save_batchJob",
                                                                        "shape_save_batchString",
                                                                        "shape_save_batchIndex",
                                                                        "shape_fileName_dataBase"))) )
    # This suppresses writting out the NULL return
    invisible( NULL )
}

###################################################################################################
##################################### BEGIN OF FUNCTIONS ##########################################
###################################################################################################
#' This is a convenience script to build an named list of empty lists, where the names are based
#' on the genotype IDs being passed.
#'
#' @param func_focalID This should be any vector, that can be interpreted as character, and faithfully represent
#' the genotype IDs of interest for your pedigree.
#' @return a named list of empty lists.
#'
#' @examples
#' # this creates a named list, this trivial function exists for future flexibility and method design.
#' buildPedigree(c(1,"zebra","walrus",4))
#'
#' @export
buildPedigree <- function(func_focalID){
  # The pedigree frame is simply a list of lists
  return( sapply(as.character(func_focalID),function(thisID){
               list()
             },simplify=FALSE) )
}


# looks through the func_pedigreeFrame (passed in full as func_pedigreeAll), and func_lineageDemographics (passed in full as func_demoArray)
# which should include all parents that have ever been a parent to a lineage which needs parenthood, tracked with values >= 0 for as long as it's important
# for focalID lineage definition!
#' This function will look through a pedigree data.frame and recursively continue building that back through the history of the SHAPE run being processed.
#'
#' @param func_focalGenotype a vector of genotype IDs whose lineage you wish to identify.
#' @param func_startStep this is the first step in the SHAPE run from which you wish to consider re-tracing the lineage.
#' @param func_stepMatrix this is the matrix that represent what happened at each step in the SHAPE run.
#' @param func_progenitorList this is a list of the known progenitor(s) for our func_focalGenotypes
#' @param func_demoArray this is the whole array of step-wise SHAPE records for population demographics and feeds func_stepMatrix.
#' @param func_pedigreeAll this is a data.frame which contains all currently known pedigree information and informs our step-wise focus.
#' @param func_lineString this is the string that will be used to collapse the vector of progenitor genotype's into a single
#' charater string.  This collapse is done as a convenience for storage and retrieval.
#'
#' @return a vector of character strings, each of which is the found lineage of the func_focalGenotypes
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
findParent <- function(func_focalGenotype, func_startStep, func_stepMatrix, func_progenitorList,
                       func_demoArray, func_pedigreeAll, func_lineString = getOption("shape_string_lineDescent")){
  # We look for the first step at which this lineage did not exist, but earlier than the starting step
  func_missingSteps <- rownames(func_stepMatrix)[which(func_stepMatrix[,"popSize"] == 0)]
  func_missingSteps <- as.numeric(nameTable_step(func_missingSteps,funcSplit = TRUE))
  # We keep only the steps which are after our start step
  func_missingSteps <- func_missingSteps[which(func_missingSteps < func_startStep)]
  # If there are no missing steps then this genotype has always existed, so we return it as the parent
  if(length(func_missingSteps) == 0){
    return( func_focalGenotype )
    # this means there is something to be done
  } else {
    # We look for the parent(s) of our genotype on the max step
    func_fociParents <- unlist(lapply(names(func_progenitorList),function(thisProgenitor){
      # If this porgenitor gave a mutant at this step we return it, else NULL
      if(is.element(max(func_missingSteps+1), func_progenitorList[[thisProgenitor]])){
        return( thisProgenitor )
      } else {
        return( NULL )
      }
    }))
    # Then we pass the max step which is missing the foci along with the parent(s) down the recursion
    # This assumes our func_stepMatrix and func_progenitorList will always have the information for all
    # parental types
    return( paste(unlist(lapply(func_fociParents,function(thisProgenitor){
      findParent(func_focalGenotype = as.numeric(thisProgenitor),
                 func_startStep = max(func_missingSteps),
                 func_stepMatrix = func_demoArray[,, thisProgenitor],
                 func_progenitorList = func_pedigreeAll[[thisProgenitor]],
                 func_demoArray = func_demoArray,
                 func_pedigreeAll = func_pedigreeAll)
    })),
    func_focalGenotype,
    sep= func_lineString) )
  }

}


#' This is a function that steps forward through time steps of a SHAPE run and extracts population demographic
#' information.  This includes Fitness, Number of Lineages, and Transitions between dominant genotypes.
#' Most important it will also return the information related to which lineages will eventually establish in the population,
#' a piece of information that will be critical for downstream lineage specific information extraction.
#'
#' @param func_stepsCon This is the filepath to an SQLite database storing information for the stepwise changes of a SHAPE run.
#' @param func_estValue This value is used to define the threshold size required for a population before it is considered established.
#' @param func_landscapeCon This is the filepath to an SQLite database storing information for the complete explored and neighbouring fitness landscape of a SHAPE run.
#' @param func_hoodCon This is the filepath to an SQLite database storing information for high priority mutational neighbourhood information
#' (which is simply a subset of the full mutational landscape).
#' @param func_size_timeStep This is the proportion of a standard biological generation which is to be simulated in a single time step.
#'
#' @return This return a list object that contains various pieces of usefull summary demographic information.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
extract_popDemographics <- function(func_stepsCon, func_estValue, func_landscapeCon, func_hoodCon, func_size_timeStep){
  # We find the ordered series of step tables that we'll be calling in this extaction process
  func_allTables <- unlist(RSQLite::dbListTables(func_stepsCon))
  # We now order that set to an ascending series
  func_allTables <- func_allTables[order(as.numeric(nameTable_step(func_allTables,funcSplit = TRUE)), decreasing = FALSE)]
  # Now we establish some reporting lists the first is the general fitness and the next the number of lineages information
  # the next being for the information regarding transition timing, the last for tracking any lineage that is every established
  tmpReturn <- list("demoMat"= matrix(-1,nrow = length(func_allTables), ncol = 8,
                                      dimnames = list(func_allTables,c("minFit","meanFit","maxFit","sdFit","numLines","numEstablished","numMutants","popSize"))),
                    "transitionMat" = matrix(-1,ncol=5,dimnames=list(NULL,c("Step","genotypeID","numMuts","fitness","transitionGens"))),
                    "vec_estLineages" = 0,
                    "vec_final_estLineages" = 0,
                    "Hists" = vector(mode="list",length=length(func_allTables)))
  # We set the names of the histograms
  names(tmpReturn[["Hists"]]) <- func_allTables

  # We now go forward through each of the time steps and we extract information
  for(thisStep in func_allTables){
    # We go through each of the steps, extracting the pertinent information
    tmpData <- dbGetQuery(func_stepsCon, paste("SELECT numMuts, genotypeID, popSize, fitness, mutants FROM ",thisStep,sep=""))
    # Any genotype with fewer than 0 individuals will be recorded as having effectively zero
    tmpData[which(tmpData[,"popSize"] <= 0),"popSize"] <- 1e-1
    # We extract all established lineages
    tmp_estLineages <- querryEstablished(func_inMatrix = tmpData, func_estProp = func_estValue)
    # If there are any estblished line(s) we will be updating this information
    if(nrow(tmp_estLineages) > 0){
      # We now keep only those which are no already in the vec_estLineages
      tmpReturn[["vec_estLineages"]] <- unique(c(tmpReturn[["vec_estLineages"]], tmp_estLineages$genotypeID))
      # If we're on the final step, we'll record which lineages were established at this point.
      if(thisStep == func_allTables[length(func_allTables)]){
        tmpReturn$vec_final_estLineages <- tmp_estLineages$genotypeID
      }
    }

    # Now we check if there has been a transition between dominant genotypes, this is recorded by looking
    # if there are any dominant lineages, and thereafter by looking at the last row
    tmp_maxLine <- tmpData[which(tmpData$popSize == max(tmpData$popSize)),]
    # We find if there are any established lineage(s) not previously recorded at the same step, meaning we need to find which
    # rows of our transition matrix relate to the last time step(s)
    tmp_transRow <- which(tmpReturn[["transitionMat"]][,"Step"] == max(tmpReturn[["transitionMat"]][,"Step"]))
    # If any of the max sized lineage(s) in this step is not the same as the current max then we update our information about the max lineages in thisStep
    if(any(!is.element(tmp_maxLine$genotypeID,tmpReturn[["transitionMat"]][tmp_transRow,"genotypeID"]))){
      tmpReturn[["transitionMat"]] <- rbind(tmpReturn[["transitionMat"]],
                                            matrix(c(rep(as.numeric(nameTable_step(thisStep,funcSplit = TRUE)),nrow(tmp_maxLine)),
                                                     tmp_maxLine$genotypeID,
                                                     tmp_maxLine$numMuts,
                                                     tmp_maxLine$fitness,
                                                     rep(as.numeric(nameTable_step(thisStep,funcSplit = TRUE)) * func_size_timeStep -
                                                           ifelse(tmpReturn[["transitionMat"]][nrow(tmpReturn[["transitionMat"]]),"Step"] != -1,tmpReturn[["transitionMat"]][nrow(tmpReturn[["transitionMat"]]),"Step"],0)* func_size_timeStep,
                                                         nrow(tmp_maxLine)) ),
                                                   nrow=nrow(tmp_maxLine)))
    }

    # We insert the information into our demoFrame
    tmpReturn[["demoMat"]][thisStep,] <- c(min(tmpData$fitness),
                                           weighted.mean(tmpData$fitness,tmpData$popSize),
                                           max(tmpData$fitness),
                                           if(nrow(tmpData) == 1){ 0 } else { sd(tmpData[,"fitness"]* tmpData[,"popSize"]/sum(tmpData[,"popSize"])) },
                                           nrow(tmpData),
                                           nrow(tmp_estLineages),
                                           sum(tmpData$mutants),
                                           sum(tmpData$popSize))

    # We now create a histogram of the population fitness, with counts considering popSize
    # I actually do the fitness histogram in a "cheaty" fashion
    tmp_fitHist <- hist(tmpData$fitness,plot=FALSE)
    # Now using the counts and breaks that exist (where breaks is always +1 longer...) we sum the popSize
    # Noting that the lower bound of the first instance of breaks is included but is exclusive for the others.
    tmp_fitHist$counts <- log10(if(length(tmp_fitHist$counts) > 1){
      sapply(1:length(tmp_fitHist$counts),function(x){
        if(x == 1){
          sum(tmpData$popSize[intersect(which(tmpData$fitness >= tmp_fitHist$breaks[x]),
                                        which(tmpData$fitness <= tmp_fitHist$breaks[x+1]))])
        } else {
          sum(tmpData$popSize[intersect(which(tmpData$fitness > tmp_fitHist$breaks[x]),
                                        which(tmpData$fitness <= tmp_fitHist$breaks[x+1]))])
        }
      })
      # This means there is one fitness bin so all popSize are the count value...
    } else {
      sum(tmpData$popSize)
    })
    # Now we build the histogram of population sizes
    tmp_sizeHist <- hist(log10(tmpData$popSize),plot = FALSE)
    # then update the counts to be a matter of log2
    tmp_sizeHist$counts <- log2(tmp_sizeHist$counts)
    # We remove any infinite values from the counts of either histogram
    tmp_fitHist$counts[which(tmp_sizeHist$counts == "-Inf")] <- -1
    tmp_sizeHist$counts[which(tmp_sizeHist$counts == "-Inf")] <- -1

    tmpReturn[["Hists"]][[thisStep]] <- list("fitness"= tmp_fitHist,
                                             "lines"= tmp_sizeHist)

  }
  # We just drop that initialising row for transitionMat
  tmpReturn[["transitionMat"]] <- matrix(tmpReturn[["transitionMat"]][-1,],ncol=ncol(tmpReturn[["transitionMat"]]), dimnames= dimnames(tmpReturn[["transitionMat"]]))

  # We now return this information
  return( tmpReturn )
}


#' This is a function to extract genotype/lineage specific information.  This info will be mostly through time style of information
#' but will also include information about it's line of descent, growth pressures pre-establishment, and population size.
#' @param func_focalID This is the vector of genotype ID(s) of the focal lineage(s) for which information is to be extracted.
#' @param func_estValue This value is used to define the threshold size required for a population before it is considered established.
#' @param func_stepsCon This is the filepath to an SQLite database storing information for the stepwise changes of a SHAPE run.
#' @param func_landscapeCon This is the filepath to an SQLite database storing information for the complete explored and neighbouring fitness landscape of a SHAPE run.
#' @param func_hoodCon This is the filepath to an SQLite database storing information for high priority mutational neighbourhood information
#' @param func_refMatrix Is a matrix of a SHAPE run's population demographics at a step in time.  I will be querried for
#' information regarding a genotype's number of mutations and fitness value.
#' of genotypes, but is not required but is also required
#' @param func_subNaming This is a logical which controls if the tables which report on all genotypes with X mutations should be
#' forced into a single table or it SHAPE is allowed to split them into multiple tables.
#' @param func_genomeLength The number of positions simulated within the individual's genomes.
#' @param func_max_numMutations The maximum number of mutations that could occur in a single mutation event -- CAUTION: This should never
#' be anything other than 1 as per how SHAPE is currently implemented.
#' @param func_allow_backMutations This is a logical toggle controlling if reversions are allowed -- meaning loss of mutations.
#' @param func_descentSep This is the standard string used to collapse line of descent information.
#' @param func_hoodExplore This is an object to control which strains we get deep neighbourhood information for
#' It should be one of \strong{"none","limited","priority","full"}
#' setting this higher will cost more and more in post analysis runtime.
#' NOTE: That use of \strong{limited} requires that you pass a func_refMatrix of expected shape (has a "genotypeID" column)!
#' @param func_stringSep A common string separator used to merge information.
#'
#' @return This returns a list object with several pieces of summary information for the focal genotype ID.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
extractInfo_focalID <- function(func_focalID, func_estValue, func_stepsCon, func_landscapeCon,func_hoodCon, func_refMatrix, func_subNaming,
                                func_genomeLength = getOption("shape_genomeLength"), func_max_numMutations = getOption("shape_max_numMutations"),
                                func_allow_backMutations = getOption("shape_allow_backMutations"), func_descentSep = getOption("shape_string_lineDescent"),
                                func_hoodExplore = getOption("shape_const_hoodDepth"), func_stringSep = getOption("shape_sepString")){
  # We find the ordered series of step tables that we'll be calling in this extaction process
  func_allTables <- unlist(RSQLite::dbListTables(func_stepsCon))
  # We now order that set to a decreasing series
  func_allTables <- func_allTables[order(as.numeric(nameTable_step(func_allTables,funcSplit = TRUE)), decreasing = TRUE)]
  # I also order the focalID's
  func_focalID <- func_focalID[order(func_focalID)]

  # We establish the reporting information objects
  func_lineageDemographics <- array(-1,dim=c(length(func_allTables),5,length(func_focalID)),
                                    dimnames=list(func_allTables,c("Step","popSize","isEstablished","births","mutsIn"),as.character(func_focalID)))
  func_pedigreeFrame <- buildPedigree(func_focalID[order(func_focalID)])

  # Now we simply step back through each of the time steps and grab information pertinent to our focal lineages
  for(thisTable in func_allTables){
    # I start by getting the table index, this will be used when we need to pass information forward or backward
    tmpIndex <- which(func_allTables == thisTable)
    # We load the step information table
    tmpData <- dbGetQuery(func_stepsCon, paste("SELECT * FROM ",thisTable," WHERE genotypeID IN (",
                                               paste(unique(names(func_pedigreeFrame)),collapse=","),
                                               ")",sep=""))
    # I also querry which line(s) are established
    tmp_estLines <- querryEstablished(func_inMatrix = tmpData, func_estProp = func_estValue)

    # My next step should be to record the popSize of that lineage in this step, is it established, any births,
    # but also, and this is the trickier part, the information concerning number of mutants which comes from progenitor info
    for(func_thisMatrix in unique(names(func_pedigreeFrame))){
      func_dataRow <- which(tmpData[,"genotypeID"] == func_thisMatrix)
      # If there is some information but this genotype is not yet in our matrix then we need update our array
      if(!is.element(func_thisMatrix, dimnames(func_lineageDemographics)[[3]])){
        func_lineageDemographics <- abind(func_lineageDemographics,
                                          array(-1,dim=c(length(func_allTables),5,length(func_thisMatrix)),
                                                dimnames=list(func_allTables,c("Step","popSize","isEstablished","births","mutsIn"), func_thisMatrix)) )

      }
      # Now we actually update our demographics but if the lineage does not exist in this step we return zero values.
      if(length(func_dataRow) > 0){
        # I pre-sum the number of mutants which were fed into this lineage as I'll be subtracting those from births so these values are separated
        func_tmpSum <- sum(unlist(lapply(strsplit(tmpData[func_dataRow,"progenitor"], getOption("shape_collapseString"))[[1]],function(x){
          as.numeric(strsplit(x,func_stringSep)[[1]][2])
        })),na.rm=TRUE)
        func_lineageDemographics[thisTable,,as.character(func_thisMatrix)] <- c(as.numeric(nameTable_step(thisTable,funcSplit = TRUE)),
                                                                                tmpData[func_dataRow,"popSize"],
                                                                                is.element(func_thisMatrix, tmp_estLines[,"genotypeID"]),
                                                                                tmpData[func_dataRow,"births"] - func_tmpSum,
                                                                                func_tmpSum)
      } else {
        func_lineageDemographics[thisTable,,as.character(func_thisMatrix)] <- c(as.numeric(nameTable_step(thisTable,funcSplit = TRUE)),
                                                                                rep(0,4))
      }
      # Now we can work on the pedigree for this lineage, we start by checking if there are any parents for the main focal lineage
      # But only if there are rows to work with
      if(length(func_dataRow) > 0){
        # We have parents if the progenitor column does not contain a blank or WT value
        if(!is.element(tmpData[func_dataRow,"progenitor"],c("","WT"))){
          # These are the genotypeIDs which contributed to our focalID
          func_tmpParents <- unlist(lapply(strsplit(tmpData[func_dataRow,"progenitor"], getOption("shape_collapseString"))[[1]],function(x){
            tmpReturn <- as.numeric(strsplit(x,func_stringSep)[[1]][1])
            # Now because a lineage may exist and receive new mutants this will create
            # situations where we get NA values returned so we return null for those parents
            if(is.na(tmpReturn) || is.null(tmpReturn)){
              return( NULL )
            } else {
              return( tmpReturn )
            }
          }))
          # Lastly we see if we need to update the pedigree data.frame for this focalID by querrying which of it's parents are not yet represented
          func_missingParents <- func_tmpParents[which(!is.element(func_tmpParents, names(func_pedigreeFrame)))]
          if(length(func_missingParents) > 0){
            func_pedigreeFrame <- c(func_pedigreeFrame,
                                    buildPedigree(func_focalID = func_missingParents))
          }
        }
      }
    }
    # Now we update our tmpData records
    tmpData <- dbGetQuery(func_stepsCon, paste("SELECT * FROM ",thisTable," WHERE genotypeID IN (",
                                               paste(unique(names(func_pedigreeFrame)),collapse=","),
                                               ")",sep=""))
    # Now for each of the focal lineages' pedigreeID's we ask ...
    for(func_thisID in unique(names(func_pedigreeFrame))){
      func_dataRow <- which(tmpData[,"genotypeID"] == as.numeric(func_thisID))
      # If we can find the row for thisID, then we want to update it provided we find some "parents"
      if (length(func_dataRow) > 0 && !is.element(tmpData[func_dataRow,"progenitor"],c("","WT"))){
        # Then we want to update the lineage of this genotype, this requires us to first define parents
        func_tmpParents <- unlist(lapply(strsplit(tmpData[func_dataRow,"progenitor"], getOption("shape_collapseString"))[[1]],function(x){
          tmpReturn <- as.numeric(strsplit(x,func_stringSep)[[1]][1])
          # Now because a lineage may exist and receive new mutants this will create
          # situations where we get NA values returned so we return null for those parents
          if(is.na(tmpReturn) || is.null(tmpReturn)){
            return( NULL )
          } else {
            return( tmpReturn )
          }
        }))
        # Right, then we're going to store parents is all the steps at which they provide mutant(s),
        # this is a convenience to reduce the length of lists, but makes later extractions perhaps more involved...
        # For each of the parents, we check if they're in the the list of thisID, and if not add them
        for(func_thisParent in as.character(func_tmpParents)){
          func_recordStep <- as.numeric(nameTable_step(thisTable,funcSplit = TRUE))
          if(is.element(func_thisParent,names(func_pedigreeFrame[[func_thisID]]))){
            func_pedigreeFrame[[func_thisID]][[func_thisParent]] <- c(func_pedigreeFrame[[func_thisID]][[func_thisParent]], func_recordStep)
          } else {
            func_pedigreeFrame[[func_thisID]][[func_thisParent]] <- func_recordStep
          }
        }
      } # This is the conditional that there is something to update and that a progenitor is one of those things...
    } # This closes the loop of looking for focalID genotypes
  } # This closes out the for loop for the time step tables

  # I define which genotypes existed on the first step, this is so that they can be ignored when they appear as root genotypes.
  func_startGenotypes <- unname(unlist(dbGetQuery(func_stepsCon, paste("SELECT genotypeID FROM ",nameTable_step(0),sep=""))))
  # Right, now what we want is the lineage for all focalID genotypes which exist on the last step AND
  # any lineages which were focal but not a first genotype.  These last would have been
  # focalID's and so must be in the pedigree information already.
  func_endFoci <- intersect(func_focalID, unname(unlist(dbGetQuery(func_stepsCon, paste("SELECT genotypeID FROM ", func_allTables[1],sep="")))))

  func_nonEndFoci <- setdiff(func_focalID, c(func_startGenotypes, func_endFoci))
  # This calls our function of finding a parent, is calls itself recursively to work through the matrix.
  func_endLineages <- list()
  for(thisFoci in func_endFoci){
    func_endLineages[[as.character(thisFoci)]] <- findParent(func_focalGenotype = thisFoci,
                                                             func_startStep = as.numeric(nameTable_step(rownames(func_lineageDemographics)[which(func_lineageDemographics[,"popSize",as.character(thisFoci)] > 0)[1]],funcSplit = TRUE)),
                                                             func_stepMatrix = func_lineageDemographics[,,as.character(thisFoci)],
                                                             func_progenitorList = func_pedigreeFrame[[as.character(thisFoci)]],
                                                             func_demoArray = func_lineageDemographics,
                                                             func_pedigreeAll = func_pedigreeFrame)
  }
  func_nonendLineages <- list()
  if(length(func_nonEndFoci) > 0){
    for(thisFoci in func_nonEndFoci){
      func_nonendLineages[[as.character(thisFoci)]] <- findParent(func_focalGenotype = thisFoci,
                                                                  func_startStep = as.numeric(nameTable_step(rownames(func_lineageDemographics)[which(func_lineageDemographics[,"popSize",as.character(thisFoci)] > 0)[1]],funcSplit = TRUE)),
                                                                  func_stepMatrix = func_lineageDemographics[,,as.character(thisFoci)],
                                                                  func_progenitorList = func_pedigreeFrame[[as.character(thisFoci)]],
                                                                  func_demoArray = func_lineageDemographics,
                                                                  func_pedigreeAll = func_pedigreeFrame)
    }
  }

  # We need to know which fitness landscape tables exist
  func_landTables <- RSQLite::dbListTables(func_landscapeCon)
  # Great, now we start building a matrix which displays neighbourhood information for all unique transitions.
  # We find a unique transition by looking at the lineages at all intermediate positions between two values.
  func_uniqueTransitions <- unique(unlist(lapply(strsplit(as.character(unique(unlist(c(func_endLineages, func_nonendLineages)))), func_descentSep),function(thisLineage){
                                            # If there is only a single element, there has been no transition...
                                            if(length(thisLineage) > 1){
                                              return( paste(thisLineage[-length(thisLineage)], thisLineage[-1],sep= func_descentSep) )
                                            } else {
                                              NULL
                                            }
                                          })))
  # Before going to get he landscape topology information I'd like to extract the numMuts and fitness value(s) for unique
  # genotypes that are involved in transitions.
  func_uniquetransitionIDs <- unique(unlist(strsplit(as.character(unique(unlist(c(func_endLineages, func_nonendLineages)))), func_descentSep)))
  func_refInfo <- matrix(c(as.numeric(func_uniquetransitionIDs),rep(-1,length(func_uniquetransitionIDs) *2)),
                         ncol=3,dimnames=list(func_uniquetransitionIDs,c("genotypeID","numMuts","fitness")))
  # Ok, now if we've got a func_refMatrix object lets populate what we can into our func_refInfo
  if(!is.null(func_refMatrix)){
    # we define the unique genotypeID in the ref object
    func_tmpIDs <- intersect(unique(func_refMatrix[,"genotypeID"]), as.numeric(func_uniquetransitionIDs))
    func_refInfo[as.character(func_tmpIDs),c("numMuts","fitness")] <- func_refMatrix[unlist(lapply(func_tmpIDs,function(x){
      which(func_refMatrix[,"genotypeID"] == x)[1]
    })),
    c("numMuts","fitness")]
  }
  # Next if we're missing the fitness OR numMuts for any of these genotypeIDs we'll go and grab that from the fitness landscape
  func_missIDs <- unname(unlist(apply(func_refInfo,MARGIN=1,function(thisLine){
    if(any(thisLine[c("genotypeID","numMuts","fitness")] == -1)){
      return( thisLine["genotypeID"] )
    } else {
      return( NULL )
    }
  })))
  # If we're not missing the information for any IDs then we can skip the following set of calls
  if(length(func_missIDs) > 0){
    # So for each of these ID's we go and query for the fitness and or numMuts information as required
    # NOTE: numMuts will be informed by binaryString length after strsplit
    # I make a single database query which may mean I'm replicating information but such is life...
    func_updateInformation <- dbGetQuery(func_landscapeCon, paste("SELECT genotypeID,binaryString,fitness FROM ",
                                                                  func_landTables,
                                                                  ' WHERE genotypeID IN(',
                                                                  paste(func_missIDs,collapse=","),
                                                                  ')',
                                                                  collapse=" UNION "))
    # Now we update the binaryString information to be genome length
    func_updateInformation[,"binaryString"] <- unlist(lapply(strsplit(func_updateInformation[,"binaryString"],func_stringSep),length))
    colnames(func_updateInformation)[which(colnames(func_updateInformation) == "binaryString")] <- "numMuts"
    # We now update the information for any genotypes that were missing information
    for(thisCol in c("numMuts","fitness")){
      func_refInfo[as.character(func_missIDs[order(func_missIDs)]),thisCol] <- func_updateInformation[order(func_updateInformation[,"genotypeID"]), thisCol]

    }
  }

  # So we now look at the neighbourhood of the parental types and find where the fitness of the offspring fits within this
  # That will tell us what was the rank of the transition.  I won't gather information about the extent of neighbourhood
  # exploration, I have no distinct use for this information and its calculation is costly with respect to time.
  func_hoodTables <- RSQLite::dbListTables(func_hoodCon)
  # If there have been no transitions we return a NULL value, else the data.frame expected
  func_rankMat <- if(is.null(func_uniqueTransitions)){
    NULL
  } else {
    t(sapply(func_uniqueTransitions,function(thisTransition){
      # The first element is the progenitor while the second element is the offspring.
      # We grab the numMuts and fitness information for each from our func_refInfo
      tmp_transitionInfo = func_refInfo[strsplit(thisTransition, func_descentSep)[[1]],]
      # We gather the mutatation range of the parental type
      tmp_mutsRange <-  max(0,(tmp_transitionInfo[1,"numMuts"] - if(func_allow_backMutations){ func_max_numMutations }else{ 0 })):
        min(func_genomeLength ,(tmp_transitionInfo[1,"numMuts"] + (func_max_numMutations)))
      # We record the binary stirng information of the parental genotype so that it is included within the search
      tmp_parent_binaryString <- retrieve_binaryString(func_genotypeID = tmp_transitionInfo[1,"genotypeID"],
                                                       func_numMuts = tmp_transitionInfo[1,"numMuts"],
                                                       func_subNaming = func_subNaming,
                                                       func_landscapeCon = func_landscapeCon)[1,"binaryString"]
      # We now find the neighbours of the parental type
      tmp_parentalHood <- unique(c(tmp_parent_binaryString,
                                   if(is.element(nameTable_neighbourhood(tmp_transitionInfo[1,"genotypeID"]), func_hoodTables)){
                                     dbGetQuery(func_hoodCon,paste("SELECT * FROM ",nameTable_neighbourhood(tmp_transitionInfo[1,"genotypeID"]),sep=""))$neighbours
                                   } else {
                                     # We define all the possible nearestNeighbours for this lineage since they've not been stored.
                                     # This includes calling for the binary string of this genotype
                                     defineNeighbours(func_tmpGenotype = tmp_parent_binaryString,
                                                           func_tmpDirection = func_allow_backMutations)
                                   }))
      # Knowing the neighbours we can gather the fitness of each, this is done as we have their binary strings
      tmp_tmpStrings <- gsub("[[:space:]]","",paste("\'", tmp_parentalHood,"\'",collapse=','))
      tmp_querryTables <- func_landTables[unique(unlist(lapply(nameTable(tmp_mutsRange),function(tmpTable){
        which(grepl(tmpTable, func_landTables))
      })))]
      tmp_hoodFit <- unlist(lapply(tmp_querryTables,function(thisTable){
        unlist(dbGetQuery(func_landscapeCon, paste('SELECT fitness FROM ',
                                                   thisTable,
                                                   ' WHERE binaryString IN (',
                                                   tmp_tmpStrings,
                                                   ')',
                                                   sep="")))
      }), use.names=FALSE)
      tmp_hoodFit <- tmp_hoodFit[order(tmp_hoodFit, decreasing = TRUE)]
      tmp_childRank <- which(tmp_hoodFit == tmp_transitionInfo[2,"fitness"])[1]
      tmp_altPaths <- if(is.element(tmp_transitionInfo[1,"fitness"], tmp_hoodFit)){
                        which(tmp_hoodFit == tmp_transitionInfo[1,"fitness"])[1] - 1
                      } else {
                        0
                      }
      ### REMOVED: That way if there are problems they ought to be highlighted....
      # This is an unnecessary sanity check in-case somehow the offspring was not found in the parental hood
      # Which ought not to be possible but it doesn't affect the interpretation of functional transitions.
      #tmp_hoodFit <- c(tmp_transitionInfo[2,"fitness"], tmp_hoodFit)

      return(	c("absRank"= tmp_childRank,
                "hoodSize"=length(tmp_parentalHood),
                "hoodMin"=min(tmp_hoodFit),
                "hoodMax"=max(tmp_hoodFit),
                "progenitor_numMuts"= tmp_transitionInfo[1,"numMuts"],
                "progenitor_fitness"= tmp_transitionInfo[1,"fitness"],
                "offspring_numMuts"= tmp_transitionInfo[2,"numMuts"],
                "offspring_fitness"= tmp_transitionInfo[2,"fitness"],
                "progenitorID"=tmp_transitionInfo[1,"genotypeID"],
                "offspringID"=tmp_transitionInfo[2,"genotypeID"],
                "num_altPaths"=  tmp_altPaths,
                "relFit_altPaths" = if(tmp_altPaths > 0 && tmp_altPaths >= tmp_childRank){
                                      calc_relativeFitness(tmp_hoodFit[1:tmp_altPaths])[tmp_childRank]
                                    } else {
                                      0
                                    },
                "prop_maxFit" = if(tmp_altPaths > 0 && tmp_altPaths >= tmp_childRank){
                                  (tmp_altPaths - (tmp_childRank -1)) / tmp_altPaths
                                } else {
                                  0
                                }) )


    }))
  }
  # We can now return our information
  return( list("lineDemo"=func_lineageDemographics,
               "linePedigree"=func_pedigreeFrame,
               "landscapeTopology"=func_rankMat,
               "end_Lines_of_Descent"=func_endLineages,
               "transition_Lines_of_Descent"=func_nonendLineages) )
}


#' This is a wrapper function to process a SHAPE run and extract meaningful summary information.
#'
#' @param func_saveFile This is the filepath where the SHAPE run processed objects are to be saved.
#' @param func_subNaming This is a logical which controls if the tables which report on all genotypes with X mutations should be
#' forced into a single table or it SHAPE is allowed to split them into multiple tables.
#' @param func_stepsCon This is the filepath to an SQLite database storing information for the stepwise changes of a SHAPE run.
#' @param func_landscapeCon This is the filepath to an SQLite database storing information for the complete explored and neighbouring fitness landscape of a SHAPE run.
#' @param func_hoodCon This is the filepath to an SQLite database storing information for high priority mutational neighbourhood information
#' @param func_estProp This value is used to define the threshold size required for a population before it is considered established.
#' @param func_size_timeStep This is the proportion of a standard biological generation being considered to be within a single time step.
#' @param func_processObjects This is a vector of character strings which define the names of what objects will be produced and creates a global objects.  DO NOT CHANGE THESE VALUES.
#' @param func_hoodPriority This is an object to control which strains we get deep neighbourhood information for
#' It should be one of \strong{"none","limited","priority","full"}
#' setting this higher will cost more and more in post analysis runtime.
#'
#' @return This returns a string vector stating the result of trying to process for the specified filepath.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
runProcessing <- function(func_saveFile, func_subNaming, func_stepsCon,
                          func_landscapeCon, func_hoodCon, func_estProp, func_size_timeStep,
                          func_processObjects = getOption("shape_processedObjects"),
                          func_hoodPriority = getOption("shape_const_hoodDepth")){
  # We check to see if there already exists a file which would have the same name as the output from this file, if not then we process.
  if(file.exists(func_saveFile)){
    return( paste("Found a pre-processed file of same expected name: ", func_saveFile," did not perform processing of run",sep="") )
  } else {
    # For each of these sapply calls, we'll nest a sapply call if a particular run was replicated

    # Next is to extract demographic data from each of our sets
    runDemographics <- extract_popDemographics(func_stepsCon = func_stepsCon, func_estValue = func_estProp,
                                               func_landscapeCon = func_landscapeCon, func_hoodCon = func_hoodCon,
                                               func_size_timeStep = func_size_timeStep)
    # Now we're extracting information regarding lines which have established in our population
    info_estLines <- extractInfo_focalID(func_focalID = runDemographics[["vec_estLineages"]],
                                         func_estValue = func_estProp,
                                         func_stepsCon = func_stepsCon,
                                         func_landscapeCon = func_landscapeCon,
                                         func_hoodCon = func_hoodCon,
                                         func_refMatrix = runDemographics[["transitionMat"]],
                                         func_subNaming = func_subNaming,
                                         func_descentSep = getOption("shape_string_lineDescent"),
                                         func_hoodExplore = func_hoodPriority)
    # We now save our objects
    save(list= func_processObjects, file= func_saveFile)
    # We report having completed
    if(file.exists(func_saveFile)){
      return( paste("Processing completed for :", func_saveFile,sep="") )
    } else {
      return( paste("Save of :", func_saveFile,", failed as file does not exist.",sep="") )
    }

  }
}


#' This is a wrapper function where the birth and death related parameters of a SHAPE run are passed before the
#' appropriate functions (and their associated methods) are called.  This function will be called once per time step
#' of a SHAPE run.
#'
#' @param func_inSize This is the vector of population sizes within the community
#' @param func_inFitness This is the vector of fitness value for the community
#' @param func_bProb This is the general bith probability defined for this run of SHAPE
#' @param func_dProb This is the general death probability defined for this run of SHAPE
#' @param func_deathDen_logical This is a logical toggle to define if deaths are calculated
#' in a density dependent manner.
#' @param func_deathDen_max This is the community size at which maximum density dependent deaths
#' (ie: 100\% of func_inSize) occur.
#' @param func_deathDen_power This is a scaling factor that controls the rate of transition between
#' minimal and maximal values of the density dependent deaths.  Higher values mean a steeper transition
#' such that there are fewer deaths until higher densities are reached.
#' @param func_sizeStep This is a proportional scalar that will control what proportion of a standard
#' "generation" is simulated for each step within a SHAPE run.  NOTE: This parameter is not perfectly
#' validated to run as may be expected with all models.  For now, it should be left as a value of "1",
#' but exists for future implementation and testing.
#' @param func_growthForm This is the implemeted growth model to be simulated in this run.  Currently
#' this can be one of \strong{"logistic","exponential","constant","poisson"}.
#' @param func_carryingCapacity This is the maximum community size supported by tge simulated environment.
#' @param func_basalRate This is the basal growth rate, otherwise definable as the number of offspring
#' an individual will produce from a single birth event.
#' @param func_deathScale This is a logical toggle to define if the number of births should be scaled
#' by the number of deaths.  The exact interpretation of this varies by growth model, but in general
#' it forces growth to follow rates expected by standard pure birth models while still simulating
#' deaths within the community.
#' @param func_drift This is a logical toggle as to whether or not stochasticity is introduced into
#' the deterministic calculations that may be encountered within the growth function.  Its exact
#' implementation varies based on the growth model being simulated.
#' @param func_roundValues This is a logical toggle to define if the number of births and deaths are
#' forced to be tracked as integer values.  If TRUE, then any fractional amounts will be stochastically
#' rounded to the nearest integer with a probability of being rounded up equal to the decimal value
#' -- ie: 0.32 means 32\% chance of being rounded up --
#' @param func_inIDs This is a vector of the genotype IDs passed to this function, its order should
#' be representative of the ordered genotypeIDs passed for func_inSize and func_inFitness.
#'
#' @return A 2 column matrix of numeric values with columns "births" and "deaths", and rownames equal
#' to func_inIDs (as.character).
#'
#' @examples
#' # Imagine you've got an evolving community of three populations where
#' # in each time step 100% of individuals die and individuals with relateive
#' # fitness of 1 produce 2 offspring.  This growth function calculates the births
#' # and deaths of that community.
#' # First I show you when births are deterministic (proof of implementation):
#' growthFunction(func_inSize = c(100,100,100), func_inFitness = c(1,2,1.05),
#'                   func_bProb = 1, func_dProb = 1,
#'                   func_sizeStep = 1, func_growthForm = "exponential",
#'                   func_drift = FALSE, func_deathScale = TRUE)
#' # Now same things but with evolutionary drift thrown in
#' growthFunction(func_inSize = c(100,100,100), func_inFitness = c(1,2,1.05),
#'                func_bProb = 1, func_dProb = 1, func_sizeStep = 1,
#'                func_growthForm = "exponential", func_drift = TRUE,
#'                func_deathScale = TRUE)
#' # Now technically the values in the birth column is really the net population
#' # size and I'd previously set the births to be scaled by deaths but if this were
#' # not the case you'd get final population sizes of:
#' growthFunction(func_inSize = c(100,100,100), func_inFitness = c(1,2,1.05),
#'                    func_bProb = 1, func_dProb = 1, func_sizeStep = 1,
#'                    func_growthForm = "exponential", func_drift = TRUE,
#'                    func_deathScale = FALSE)
#'
#' @export
growthFunction <- function(func_inSize, func_inFitness, func_bProb, func_dProb, func_deathDen_logical = FALSE, func_deathDen_max = NULL,
                           func_deathDen_power = 4,func_sizeStep, func_growthForm = c("logistic","exponential","constant","poisson"),
                           func_carryingCapacity = NULL, func_basalRate = NULL, func_deathScale = FALSE,
                           func_drift = TRUE, func_roundValues = FALSE, func_inIDs = NULL){
  # This function will call the birth and death functions in proper order based on the parameterisation of the
  # growth model and how density dependence is being considered.

  # We start by triming down to the first growth form values passed, or this means that it will default to logistic
  # then we check that what is being used is something that can be recognised, else we complain.
  func_growthForm <- func_growthForm[1]
  if(!is.element(func_growthForm[1],
                  eval(as.list(args(growthFunction))$func_growthForm))){
    stop(paste("Incorrect growth form was passed as being ", func_growthForm," please review",sep=""))
  }
  # Other sanity check, if the vector of populations is all zero, we throw this as well
  if(all(func_inSize == 0)){
    # No individual means no births.... end of story.
    return( func_inSize )
  } else if( any(is.element(func_inSize,c("NaN","-Inf","Inf"))) ){
    # We complain about this, though it shouldn't be a problem...
    stop("One of the population elements had NaN, or Inf size, please review")
  }
  # We calculate the deaths first
  func_tmpDeaths <- deathFunction(func_inSize = func_inSize,
                                  func_inProb = func_dProb * func_sizeStep,
                                  func_roundValues = func_roundValues,
                                  func_depDensity = func_deathDen_logical,
                                  func_densityMax = func_deathDen_max,
                                  func_densityPower = func_deathDen_power)
  # then we calculate the births
  func_tmpBirths <- birthFunction(func_inSize = func_inSize,
                                  func_inFitness = func_inFitness,
                                  func_bProb = func_bProb,
                                  func_sizeStep = func_sizeStep,
                                  func_growthForm = func_growthForm,
                                  func_deaths = func_tmpDeaths,
                                  func_carryingCapacity = func_carryingCapacity,
                                  func_basalRate = func_basalRate,
                                  func_deathScale = func_deathScale,
                                  func_drift = func_drift,
                                  func_roundValues = func_roundValues)
  # We return the matrix of births and deaths with dimension names
  return( matrix(c(func_tmpBirths,func_tmpDeaths),ncol=2,
                 dimnames=list(func_inIDs,c("births","deaths"))) )
}



# This function takes the simulation run's parameters for type of growth, as well as the population's data, and returns the amount of growth
# for each lineage.  As values of growth calculated may include decimals, and this simulation tool deals in whole individuals, we round
# the values in a stochastic fashion with probability == size of the remainder.  Rounding is done in this way as to not prejudice small populations.
# This mehtod is prepared to handle values such as exponential, logistic, or constant where logistic is the default.
# the second line of variables are optinal and depend on which growth type is being called.  The func_deathScale argument is to control if
# we scale the number of births by the number of deaths, if this is done then the growth between time steps should match up to the analytical solutions
# of the respective growth forms; NOTE: That constant is not affected since it already considers deaths, nor is Poisson since it's growth form
# is simply a statement that the number of offspring is poisson distributed, and is not in fact explicitly a function about growth.
# The func_drift argument is a logical toggle for selecting if births are deterministic or simply an expectation that will be passed through a Poisson call.
# Last part is the control of rounding values, this affects whether or not values will be integer or not.
#' This function calculates the number of births for the vector of populations which are expected to be passed.
#' The number of parameters which can be passed may be more than the number required to use one of the growth forms.
#'
#' @param func_inSize This is the vector of population sizes within the community
#' @param func_inFitness This is the vector of fitness value for the community
#' @param func_bProb This is the general bith probability defined for this run of SHAPE
#' @param func_sizeStep This is a proportional scalar that will control what proportion of a standard
#' "generation" is simulated for each step within a SHAPE run.  NOTE: This parameter is not perfectly
#' validated to run as may be expected with all models.  For now, it should be left as a value of "1",
#' but exists for future implementation and testing.
#' @param func_growthForm This is the implemeted growth model to be simulated in this run.  Currently
#' this can be one of \strong{"logistic","exponential","constant","poisson"}.
#' @param func_deaths This is the vector of deaths for the genotypes within the community
#' @param func_carryingCapacity This is the maximum community size supported by tge simulated environment.
#' @param func_basalRate This is the basal growth rate, otherwise definable as the number of offspring
#' an individual will produce from a single birth event.
#' @param func_deathScale This is a logical toggle to define if the number of births should be scaled
#' by the number of deaths.  The exact interpretation of this varies by growth model, but in general
#' it forces growth to follow rates expected by standard pure birth models while still simulating
#' deaths within the community.
#' @param func_drift This is a logical toggle as to whether or not stochasticity is introduced into
#' the deterministic calculations that may be encountered within the growth function.  Its exact
#' implementation varies based on the growth model being simulated.
#' @param func_roundValues This is a logical toggle to define if the number of births and deaths are
#' forced to be tracked as integer values.  If TRUE, then any fractional amounts will be stochastically
#' rounded to the nearest integer with a probability of being rounded up equal to the decimal value
#' -- ie: 0.32 means 32\% chance of being rounded up --
#'
#' @return A vector of births with the same length as the vector of population sizes passed.
#'
#' # Imagine you've got an evolving community of three populations where in each time step individuals with
#' # relateive fitness of 1 produce 2 offspring.
#' birthFunction(func_inSize = c(100,100,100), func_inFitness = c(1,2,1.05), func_bProb = 1,
#' func_sizeStep = 1, func_growthForm = "exponential", func_drift = FALSE)
#' # Now with evolutionary drift
#' birthFunction(func_inSize = c(100,100,100), func_inFitness = c(1,2,1.05), func_bProb = 1,
#' func_sizeStep = 1, func_growthForm = "exponential", func_drift = TRUE)
#'
#' @export
birthFunction <- function(func_inSize, func_inFitness, func_bProb, func_sizeStep, func_growthForm = c("logistic","exponential","constant","poisson"),
                          func_deaths = NULL, func_carryingCapacity = NULL, func_basalRate = NULL, func_deathScale = FALSE, func_drift = TRUE,
                          func_roundValues = TRUE){
  # Great now we handle growth differently based on which form was called
  func_numBirths <- if(func_growthForm == "exponential"){
    # this is simply exponential growth for our lineage(s) which can be computed by the increase in individuals
    # And we remove the original amount of each population since this value will be added to the existing numbers.
    tmpBirths <- expGrowth(func_rate= func_inFitness * func_bProb *
                             if(!is.null(func_basalRate)){func_basalRate}else{2},
                           func_step= func_sizeStep,
                           func_startPop= func_inSize) -
      (func_inSize - if(func_deathScale){func_deaths}else{0})
    if(func_drift){
      tmpUpdate <- which(tmpBirths %%1 != 0)
      tmpBirths[tmpUpdate] <- addDrift(tmpBirths[tmpUpdate], func_integerValues = func_roundValues)
    }
    tmpBirths
  } else if (func_growthForm == "logistic"){
    ## Now we calculate the growth rate, per lineage, which is modified by the density dependent logisticMap term
    func_densityTerm <- (logisticMap(func_rate= func_inFitness * func_bProb * func_sizeStep *
                                       if(!is.null(func_basalRate)){func_basalRate}else{2},
                                     func_startPop= sum(func_inSize),
                                     func_maxPop= func_carryingCapacity)/sum(func_inSize)) - 1
    # Now we computed the actual number of newborn, if asked to compensate for the dead then we add back in a number of newborn
    # given the number of dead as determined by a call to growthFunction with constant growth form
    tmpBirths <- (func_inSize * func_densityTerm) +
      if(func_deathScale && !all(func_deaths == 0)){
        birthFunction(func_inSize = func_inSize,
                      func_inFitness = func_inFitness,
                      func_bProb = func_bProb,
                      func_sizeStep = func_sizeStep,
                      func_growthForm = "constant",
                      func_deaths = func_deaths,
                      func_carryingCapacity = func_carryingCapacity,
                      func_basalRate = func_basalRate,
                      func_deathScale = FALSE,
                      func_drift = FALSE,
                      func_roundValues = func_roundValues)
      }else{
        0
      }
    if(func_drift){
      tmpUpdate <- which(tmpBirths %%1 != 0)
      tmpBirths[tmpUpdate] <- addDrift(tmpBirths[tmpUpdate], func_integerValues = func_roundValues)
    }
    tmpBirths
  } else if (func_growthForm == "poisson"){
    # This is Poisson growth, where an individual has a number of offspring drawn from the poisson distribution
    # this is meant to keep a population approximately equal so long as prob_b == prob_d.
    unlist(lapply(func_inSize * func_inFitness * func_bProb * func_sizeStep,function(thisLineage){ rpois(1,thisLineage) }))
  } else if(func_growthForm == "constant"){
    # To consider constant growth, I calculate the deterministic growth, then if there is drift this value gets
    # passed to the drift function with basis on rounding.  These values are then scaled back to the number of deaths
    tmpBirths <- unlist(lapply(func_inSize * calc_relativeFitness(func_fitVector = func_inFitness) * func_bProb * func_sizeStep,
                               function(thisLineage){
                                 if(func_drift){
                                   addDrift(thisLineage, func_integerValues = func_roundValues)
                                 } else {
                                   thisLineage
                                 } }))
    # If there were no births we can simply return that vector of zeroes...
    if(sum(tmpBirths) != 0){
      tmpBirths <- tmpBirths/sum(tmpBirths) * sum(func_deaths)
    }
    tmpBirths

  }

  # We round our population values stochastically, if we're told to round values to integers, and if the value is not already an integer
  # - found to save time, and becuase above when we add drift, the values may already be integers if we're roundign values (due to Poisson stochasticity)
  if(func_roundValues){
    tmpUpdate <- which(func_numBirths %%1 != 0)
    if(length(tmpUpdate) > 0){
      func_numBirths[tmpUpdate] <- unlist(lapply(func_numBirths[tmpUpdate], function(thisValue){
        # We handle that the number of births may be a negative value
        if(sign(thisValue) == 1){
          floor(thisValue) + rbinom(1,1, thisValue %%1)
        } else {
          ceiling(thisValue) - rbinom(1,1, thisValue %%1)
        }
      }))
    }
  }
  # Now we adjust the number of births calculated, this is likely not actually required except that there is some stochastic rounding
  # which occurs meaning that we may overshoot population size.
  ### NOTE: this obviously only matters for constant growth
  if(func_growthForm == "logistic"){
    # We only need to adjust the size if the proposed growth is not ba;anced AND there is at least some proposed growth
    if((sum(func_inSize, func_numBirths) - if(func_deathScale){sum(func_deaths)}else{0}) > func_carryingCapacity && !all(func_numBirths == 0)){
      func_numBirths <- adjustBirths(func_adjVector= func_numBirths,
                                     func_sumTotal= func_carryingCapacity - sum(func_inSize) + if(func_deathScale){sum(func_deaths)}else{0},
                                     func_roundValues = func_roundValues)
    }
  } else if(func_growthForm == "constant"){
    if(sum(func_numBirths) != sum(func_deaths) && !all(func_numBirths == 0)){
      func_numBirths <- adjustBirths(func_adjVector= func_numBirths,
                                     func_sumTotal= sum(func_deaths),
                                     func_roundValues = func_roundValues)
    }
  }
  # We now return the new lineage sizes given the births
  return( func_numBirths )
}


# If the death rate is density dependent then a value must be passed to func_depDensity or we simply use the full func_inProb value
#' This allows SHAPE to simulate the death process as a deterministic value, and may be density dependent.
#' @param func_inSize This is the vector of population sizes within the community
#' @param func_inProb This is the general death probability defined for this run of SHAPE
#' @param func_roundValues This is a logical toggle to define if the number of births and deaths are
#' forced to be tracked as integer values.  If TRUE, then any fractional amounts will be stochastically
#' rounded to the nearest integer with a probability of being rounded up equal to the decimal value
#' -- ie: 0.32 means 32\% chance of being rounded up --
#' @param func_depDensity This is a logical toggle as to whether or not the calculation is density dependent.
#' If TRUE, then func_densityMax reuqires a value.
#' @param func_densityMax This is the community size at which maximum density dependent deaths
#' (ie: 100\% of func_inSize) occur.
#' @param func_densityPower This is a scaling factor that controls the rate of transition between
#' minimal and maximal values of the density dependent deaths.  Higher values mean a steeper transition
#' such that there are fewer deaths until higher densities are reached.
#'
#' @return A vector of the number of deaths caluclated for each of the populations represented by the func_inSize vector
#'
#' @examples
#' # Imagine you've got an evolving community of three populations where in each time step
#' # 100% of individuals die.
#' deathFunction(func_inSize = c(100,50,200), func_inProb = 1)
#' # What if their deaths were scaled based on population density,
#' # or an environmental carrying capacity?
#' deathFunction(func_inSize = c(100,50,200), func_inProb = 1,
#'               func_depDensity = TRUE, func_densityMax = 400)
#' deathFunction(func_inSize = c(100,50,200), func_inProb = 1,
#'               func_depDensity = TRUE, func_densityMax = 500)
#' deathFunction(func_inSize = c(100,50,200), func_inProb = 1,
#'               func_depDensity = TRUE, func_densityMax = 350)
#'
#' @export
deathFunction <- function(func_inSize, func_inProb = 0, func_roundValues = TRUE,
                          func_depDensity = FALSE, func_densityMax = NULL, func_densityPower = 4){
  # The number of death events defined deterministically using the probability of death and the popualion size, with rounding being
  # the only stochastic portion.  However, we can ignore any zero value populations, they won't have further deaths

  # Here we check if the calculation is density dependent or not
  func_tmpReturn <- if(func_depDensity && !is.null(func_densityMax)){
    func_inSize * func_inProb * (sum(func_inSize)/func_densityMax)^ func_densityPower
  } else {
    func_inSize * func_inProb
  }

  # this is where the stochastic rounding is performed.
  if(func_roundValues && any(func_tmpReturn %% 1 != 0)){
    func_tmpUpdates <- which(func_tmpReturn %% 1 != 0)
    func_tmpReturn[func_tmpUpdates] <- unlist(lapply(func_tmpReturn[func_tmpUpdates],function(thisPop){
      floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1)
    }))
  }
  return( func_tmpReturn )
}

#' This allows SHAPE to simulate the mutation process as a deterministic value.  At present, values must be tracked as integer results
#' for reasons of how I am passing to functions which identify what mutant genotype(s) are created.
#'
#' @param func_inSize This is the vector of the population sizes, or perhaps number of births, or sum of both, within the community.  Which vector
#' gets passed will depend on which growth form and other parameters are being implemented by SHAPE.
#' @param func_inProb This is the general mutation rate (probability) defined for this run of SHAPE. It is a per individual considered value, by
#' which I mean that each mutant will have a single new mutation (or reversion if allowed - handled elsewhere) and so this probability is based
#' on the vector of individuals passed and any context of if it is a "per generation" value relates to how time steps and birth probabilities are handled in the run.
#'
#' @return A vector of the number of mutants produced by each of the populations represented by the func_inSize vector
#'
#' @examples
#' # The number of mutants generated is forcibly integer but is based
#' # on the stochastic rounding of the product of the number of potentially
#' # mutable individuals and their probability of mutation.
#' mutationFunction(c(10,50,100),func_inProb = 0.3)
#' replicate(5,mutationFunction(c(10,50,100),func_inProb = 0.35))
#'
#' @export
mutationFunction <- function(func_inSize, func_inProb = 0){
  # Calcualte the deterministic value and then update if we are rounding.
  func_tmpReturn <- func_inSize * func_inProb
  # Here is where I enforce updating of value to be integer.
  func_tmpUpdates <- which(func_tmpReturn %% 1 != 0)
  func_tmpReturn[func_tmpUpdates] <- unlist(lapply(func_tmpReturn[func_tmpUpdates],function(thisPop){
    floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1)
  }))

  return( func_tmpReturn )
}

#' This is a simple little function used to represent drift by introducing stochasticity to the vector passed
#' by making poisson distribution calls.  At present it forces values to integers because I've not been able to
#' implement an appropriate continuous distribution for such calls that works with tested models and expected outcome.
#' @param func_inVector A vector of value to which stochasticity is to be added, integer values will be returned.
#' @param func_integerValues Logical toggle if a discrete or continous distribution is to be used for draws. DISABLED - as
#' testing could not identify a continuous distribution which works for obtaining expected results from established models.
#'
#' @return A vector of values, with same length as func_inVector
#'
#' @examples
#' # This adds drift by making draws from the Poisson distribution with a location parameter based on
#' # the elements to which drift is to be added.
#' replicate(10,addDrift(c(0.5,1,5,10,14.1)))
#'
#' @export
addDrift <- function(func_inVector, func_integerValues = TRUE){
  # At present I don't have a proper continuous Poisson like distribution from which to sample
  # Thus I suppress the ability to not round values.
  func_integerValues = TRUE ### THIS EXISTS TO FORCE RPOIS AS A COMPARISON DUE TO LACK OF CONTINUOUS POISSON LIKE DISTRIBUTION TO SAMPLE

  # If we're rounding to integer values then we use the poisson distribution, else we use gamma with default rate parameter of 1
  sign(func_inVector) * unlist(lapply(abs(func_inVector),function(thisValue){
    if(func_integerValues){
      rpois(1,thisValue)
    } else {
      ### TEST HAVE SHOWN THIS TO NOT WORK WELL
      func_tmpShape <- thisValue
      func_tmpRate <- 1 + (1.6/func_tmpShape)
      rgamma(1, func_tmpShape * func_tmpRate, func_tmpRate)
    }
  }))
}


#' This function is simply an implementation of the logistic growth equation where:
#' f(x) = K / (1 + ((K - N_0)/N_0) *exp-k(x-x_0))  ; Where x_0 is an adjustment to the position of the midpoint of the curve's maximum value
#' K = the curves maximum value, k = the steepness of the curve (growth rate), and N_0 is the starting population
#' it includes parameters to change the midpoint as well as change the natural exponent (ie: exp) to some other value.
#' NOTE: This is for continuous growth, and since SHAPE is discrete at present this is an unused function.
#'
#' @param func_rate The basal growth rate of individuals in the SHAPE run.
#' @param func_step This is the number of steps forward for which you wish to calculate the growth expected.
#' @param func_startPop The sum of the populations in the evolving community.
#' @param func_maxPop The carrying capacity of the enviromment being simulated.
#' @param func_midAdjust The midpoint which controls the point of inflection for the logistic equation.  Beware, change this at your own risk
#' as its impact will varrying based on the population sizes being simulated.  Ideally, don't change this value from its default.
#' @param func_basalExponent This defaults as the natural exponent "e" / "exp".  Change it at your own risk.
#'
#' @return Returns a single value representing the amount of logistic growth expected by the community
#'
#' @examples
#' # This calculates logistic growth based on the mathematical continuous time algorithm
#' logisticGrowth(func_rate = 2, func_step = 1, func_startPop = 1e2, func_maxPop = 1e4)
#' # It normally takes log2(D) steps for a binary fission population to reach carrying capacity,
#' # where D is max/start, in this case D = 100 and so it should take ~ 6.64 turns
#' logisticGrowth(func_rate = 2, func_step = c(1,2,3,6,6.64,7), func_startPop = 1e2, func_maxPop = 1e4)
#'
#' @export
logisticGrowth <- function(func_rate, func_step, func_startPop = NULL, func_maxPop = NULL, func_midAdjust = 0, func_basalExponent = exp(1)){
  return( func_maxPop / (1+ (((func_maxPop - func_startPop)/(func_startPop)) * func_basalExponent^(-func_rate * (func_step - func_midAdjust))) ) )
}


#' This is the discrete time logistic growth function known as the logistic map.  It calculates the amount of growth expected
#' in a step of time given by:  N_t+1 = N_t + r * (N_t (K - N_t)/K);
#' where N_t is community size at a time point, r is the per step growth rate, and K is the environmental carrying capacity.
#'
#' @param func_rate Per time step intrinsic growth rate of individuals
#' @param func_startPop The initial summed size of the evolving community
#' @param func_maxPop The carrying capacity of the simulated environment
#'
#' @return A single value as to the expected summed size of evolving populations in the considered environment.
#'
#' @examples
#' # This is the discrete time step form of the logistic equation, known as the logistic map.
#' # It takes a growth rate starting and max possible community size.
#' stepwise_Size <- 100
#' for(thisStep in 1:7){
#'   stepwise_Size <- c(stepwise_Size,
#'                      logisticMap(2,stepwise_Size[length(stepwise_Size)],1e4))
#' }
#' stepwise_Size
#' # When a population overshoots, it will loose members.
#' @export
logisticMap <- function(func_rate, func_startPop, func_maxPop){
  func_startPop + func_rate * (func_startPop * ((func_maxPop-func_startPop)/func_maxPop))
}



# This function will allow the starting or ending population size to be determined
# based on a given rate per step, step number, either  the start or stop population size,
# and with the option to set the basal value to exponentiate (it defaults to e making this an exponential growth function)
#' This function uses the exponential growth model and can either calculated the expected growth for a single time step
#' OR it can work backwards to calculated what was the expected starting population size prior to a step of exponential growth.
#'
#' @param func_rate This is the number of offpsring expected to be produced by an individual.  When calculating the expected
#' population size after a time step, we force this rate to be no less than 1 since this function has meaning only in the birth
#' function and so we do not want to calculate negative births (which would mean deaths).
#' @param func_step This is a proportional scalar that will control what proportion of a standard
#' "generation" is simulated for each step within a SHAPE run.  NOTE: This parameter is not perfectly
#' validated to run as may be expected with all models.  For now, it should be left as a value of "1",
#' but exists for future implementation and testing.
#' @param func_startPop This is the initial population size(s) for which you want to calculate a final size.  Leave NULL
#' if trying to calculated the expected initial size from a final population.
#' @param func_endPop This is the final population size(s) for which you want to calculate a initial size.  Leave NULL
#' if trying to calculated the expected final size from an initial population.
#'
#' @return numeric value
#'
#' @examples
#' # Exponential growth equation implemented but allowing either the final or initial population
#' # to be calculated based on whethere the initial or final community size is input.
#' expGrowth(func_rate = 2, func_step = 1,func_startPop = 100)
#' expGrowth(func_rate = 2, func_step = 1,func_endPop = 200)
#' expGrowth(func_rate = 2, func_step = 7,func_startPop = 100)
#' # You cannot set a growth rate less than 1 as this would then simulate deaths which is not
#' # allowed in this calculation.
#' expGrowth(func_rate = c(0.9,1,1.1), func_step = 1,func_startPop = 100)
#' @export
expGrowth <- function(func_rate, func_step,func_startPop = NULL, func_endPop = NULL){
  func_tmpRate <- (func_rate ^ func_step)
  # If the user has defined a starting population value we assume we're calculating the final value
  if(is.null(func_endPop) && !is.null(func_startPop)){
    if(any(func_tmpRate < 1)){ func_tmpRate[which(func_tmpRate < 1)] <- 1 }
    return( func_startPop * func_tmpRate )

    # If the user has entered an end value we assume we're calculating the starting value
  } else if(!is.null(func_endPop) && is.null(func_startPop)){
    return( func_endPop/func_tmpRate )
  } else {
    # otherwise we complain that this function is not meant to perform calculation of other terms.
    stop(paste("Not certain which value you wanted calculated since you defined a starting pop of ",func_startPop," and final pop of ",func_endPop," please review",sep=""))
  }
}


#' This function ensures that a vector of values will sum to a given number.  It's implemented in certain growth forms
#' (curently: \strong{constant} and \strong{logistic})
#'
#' @param func_adjVector Vector of values which must sum to the func_sumTotal.
#' @param func_sumTotal A single integer value which is to be the target summed value.
#' @param func_roundValues Logical toggle to control in values must be rounded to integers.
#'
#' @return A vector of values adjusted to sum to a single value.  These may have been forced to be rounded or could still
#' contain decimals.
#'
#' @examples
#' # In the event we're enforcing a vector to sum to a particular value, this function will
#' # force that vector to the sum and adjust proportionally to elements.  You can force values
#' # to become integers.
#' adjustBirths(func_adjVector = c(9,70,20), func_sumTotal = 100, func_roundValues = FALSE)
#' # When rounding, this is stochastic
#' replicate(10,adjustBirths(func_adjVector = c(9,70,20),
#'                           func_sumTotal = 100,
#'                           func_roundValues = TRUE))
#' # Same idea, different input vectors
#' adjustBirths(func_adjVector = c(10,75,20), func_sumTotal = 100, func_roundValues = FALSE)
#' replicate(10,adjustBirths(func_adjVector = c(10,75,20),
#'                              func_sumTotal = 100,
#'                              func_roundValues = TRUE))
#'
#' @export
adjustBirths <- function(func_adjVector, func_sumTotal, func_roundValues = getOption("shape_track_asWhole")){
  # This is a likely redundant validation as we ought not to be in this function unless this is true, but....
  if(sum(func_adjVector) != func_sumTotal){
    func_adjVector <- func_adjVector/sum(func_adjVector) * sum(func_sumTotal)
  }
  # If we're meant to round our values then we do so
  if(func_roundValues && any(func_adjVector %% 1 != 0)){
    func_tmpUpdates <- which(func_adjVector %% 1 != 0)
    func_adjVector[func_tmpUpdates] <- unlist(lapply(func_adjVector[func_tmpUpdates],function(thisPop){
      floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1)
    }))
  }
  # Due to rounding the value may anew not sum properly to the permissible value, but the problem ought not to be possible unless rounding is considered!
  # Hence our resolution to this using integer adjustments, but as a safety we won't bother if the difference is less than 1
  if(abs(diff(c(sum(func_adjVector), func_sumTotal))) >= 1){
    # The adjustment will be the difference between our permissible and proposed values
    func_adjustments <- func_sumTotal - sum(func_adjVector)
    # we distribute our adjustment randomly among the vector indexes with probability proportional to their relative to the total
    # Also, because this is an adjustment to the number of births, and we're doing integer adjustments, we won't adjust any values that are less than 1
    func_tmpAdjustable <- which(abs(func_adjVector) >= 1)
    func_proposedAdjustments <- NULL
    func_tmpCounter <- 0
    # Now, for as long as the porposed number of adjustments does not equal the change required, we recalculate, unless we've tried this 100 times in
    # which case we complain and break.
    while(sum(func_proposedAdjustments) != abs(func_adjustments) && func_adjustments >= 1 && func_tmpCounter < 100){
      # We propose adjustments
      func_proposedAdjustments <- table(sample(func_tmpAdjustable,
                                               abs(func_adjustments),
                                               replace=TRUE,
                                               prob = abs(func_adjVector[func_tmpAdjustable]/sum(func_adjVector[func_tmpAdjustable]))))
      # We reduce adjustments to be be no greater than the value passed
      func_proposedAdjustments <- sapply(names(func_proposedAdjustments) ,function(thisNamed){
        return( min(func_adjVector[as.numeric(thisNamed)], func_proposedAdjustments[thisNamed]) )
      })
      # This is an object format work-around for times when there are no proposed adjustments
      if(length(func_proposedAdjustments) == 0){ func_proposedAdjustments <- NULL }

      # We update our counter
      func_tmpCounter <- func_tmpCounter + 1
    }
    # We report if the counter was reached to break this chain
    if(func_tmpCounter >= 100){
      stop("While adjustingBirths and trying to aportion the adjustments got stuck in a loop for 100 replicates, please review...")
    }
    if(!is.null(func_proposedAdjustments)){
      # Otherwise we apply this proposed adjustment by either adding or subtracting the sampled numbers based on the sign
      func_adjVector[as.numeric(names(func_proposedAdjustments))] <- func_adjVector[as.numeric(names(func_proposedAdjustments))] +
        if(func_adjustments > 0){
          as.vector(func_proposedAdjustments)
        } else {
          - as.vector(func_proposedAdjustments)
        }
    }
  }
  # We now return our vector which should now be adjusted.
  return( func_adjVector )
}




# Now we need a function to generate the fitness values for a vector of genotypes that are passed along.
# We also provide a means to send forward the fitness of focal genotype(s) which may be the immedtae ancestor (Additive requires),
# an optimal genotype (RMF required), the independent fitness of all sites (NK required)
# (some models care about this, e.g. Fisher's Geometric....)
### NOTE: As we do not define all genotypes in the space (due to the genotype space of large genomes being impractical)
###			we can only reasonably implement landscape models which are not based on the fitness of anything other than the local space
#' This function will calculate the fitness values for genotypes being newly recorded to the fitness landscape.
#'
#' @param tmpGenotypes This is a vector of the binaryString values that represent the genotype(s) for which you want
#' to calculate new fitness values.
#' @param tmp_focalFitness This argument has different meaning depending upon the fitness landscape model being simulated.
#' It can be a vector of fitness values, a matrix, a single value, etc...
#' @param landscapeModel This is the character string that defines the fitness landscape model being simulated in this SHAPE run.
#' At present it can be one of: \strong{Additive, Fixed, HoC, NK, RMF}
#' @param tmp_ancestralFitness This is the fitness value of the pure WT genotype, it does not always have meaning.
#' @param tmp_weightsRMF This is the weighting of the constant/deterministic term calculated in the RMF fitness landscape equation.
#' @param tmp_optimaRMF This is the binary string genotype of the optimal genotype in the current RMF fitness landscape.  It needn't yet have
#' been yet explored, it is simply the genotype that will be the deterministic global optimum.
#' @param tmp_correlationsNK This is the matrix of fitness values and interactions between mutational states for the NK fitness lanscape model
#' @param tmp_const_numInteractionsNK This is the "K" value of the NK fitness landscape value and represents the number of other sites
#' correlated to the fitness of a focal site.
#' @param tmp_NK_ancestDep This is the fitness value of the WT mutant for an NK fitness landscape, it is passed as a computational
#' ease so that it needn't be calculated each time this function is called.
#' @param relativeFitness This is a logical toggle controlling if the fitness values returned should be relative fitness values
#' @param func_genomeLength This is the genome length of individuals.
#' @param func_distribution This is a character string representing which of the allowed distribution functions can be called
#' for draws of stochastic values when calculating fitness values.  See fitnessDist for those implemented.
#' @param func_distParameters This is a vector of the ordered distribution parameters expected by the distribution
#' referenced by func_distribution
#' @param func_distAsS This is a logical toggle to control in the final returned values should be considered as selection coefficients,
#' which is achieved by subtracting the calculated value by 1.
#' @param func_sepString This is a character string used for collapsing vectors of information, and expanding the collpased information back into
#' a vector of values.
#'
#' @return A vector of fitness values to be assgined for each of the newly explored genotypes defined in the vector tmpGenotypes
#'
#' @section Note:
#' There is no example as this does not have meaning outisde of a runSHAPE call.
#'
#' @export
fitnessLandscape <- function(tmpGenotypes, tmp_focalFitness, landscapeModel = "HoC",
                             tmp_ancestralFitness = getOption("shape_const_ancestFitness"), tmp_weightsRMF = getOption("shape_const_RMF_theta"),
                             tmp_optimaRMF = getOption("shape_const_RMF_globalOptima"), tmp_correlationsNK = getOption("shape_const_NK_interactionMat"),
                             tmp_const_numInteractionsNK = getOption("shape_const_numInteractions"), tmp_NK_ancestDep = getOption("shape_const_DepbySite_ancestFitness"),
                             relativeFitness = TRUE, func_genomeLength = getOption("shape_genomeLength"),
                             func_distribution = getOption("shape_constDist"), func_distParameters = getOption("shape_const_distParameters"),
                             func_distAsS = getOption("shape_const_distAsS"), func_sepString = getOption("shape_sepString")){
  #startTime <- proc.time() # This is for reporting on run times but is not needed....
  # We create a list of the mutations in each of the genotypes having been passed, only if not the "Fixed" model...
  if(landscapeModel != "Fixed"){
    tmp_genotypeList <- lapply(strsplit(tmpGenotypes,func_sepString),as.numeric)
  }
  # This is a function of small convenience to convert the state of the tmp_optimaRMF object
  if(!is.na(tmp_optimaRMF) && is.character(tmp_optimaRMF)){ tmp_optimaRMF <- as.numeric(strsplit(tmp_optimaRMF,func_sepString)[[1]]) }

  # This creates the overall fitness associated to mutations
  fitnessVec <- if(landscapeModel == "HoC"){
    # For the HoC since each site is genotype is uncorrelated we will simply generate a random number for each.
    tmpReturn <- fitnessDist(length(tmp_genotypeList), tmpDistribution = func_distribution, tmpParameters = func_distParameters)
    # If the distribution draws relative fitnesses but needs to be used as selection coefficients...
    if(func_distAsS){
      tmpReturn <- tmpReturn - 1
    }
    tmpReturn
    # Now we will deal with the instance of a simple additive model
  } else if(landscapeModel == "Additive") {
    # In the additive model there is no interaction, so we find the number and position of mutations in the genotype(s) - in tmp_genotypeList -
    # then review the constant independent fitness per site matrix which should be passed similar to the NK model.
    as.vector(sapply(tmp_genotypeList,function(thisGenotype){
      # So each genotype is already split into it's constituent mutations, so we simply add the values using the
      # siteBystate matrix which should be represented in tmp_focalFitness
      sum(c(tmp_ancestralFitness,  tmp_focalFitness[setdiff(thisGenotype,1:nrow(tmp_focalFitness)),"0"], tmp_focalFitness[thisGenotype,"1"]))
    }))
  } else if(landscapeModel == "NK") {
    ###### This is as per the model information in Kauffman 1989
    # In this case we simply need to know what is the independent fitness of all sites, in each state, thus tmp_focalFitness should be a
    # matrix with as many columns as there are states each site can assume, as many rows as genome length, and all states possible are included.
    if(!is.matrix(tmp_focalFitness) ||
        ncol(tmp_focalFitness) != length(getOption("shape_const_siteStates")) ||
        nrow(tmp_focalFitness) != func_genomeLength ||
        !all(is.element(getOption("shape_const_siteStates"),colnames(tmp_focalFitness))) ){
      stop("The number of states, and/or the shape of the tmp_focalFitness passed to a call for fitnessLandscape function, is not correct please review")
    }
    # AND we need to know which others sites each focal site depends upon.  BUT if the interactions == 0 then we can't ahve been passed a matrix
    # and thus tmp_correlationsNK should remain as a NULL value which can be handled downstream
    if(tmp_const_numInteractionsNK != 0){
      if(!is.matrix(tmp_correlationsNK) || ncol(tmp_correlationsNK) != tmp_const_numInteractionsNK || nrow(tmp_correlationsNK) != func_genomeLength){
        stop("The shape of the tmp_correlationsNK passed to a call for fitnessLandscape function, is not a matrix and or does not have as many columns as suggested by tmp_const_numInteractionsNK, please review")
      }
    } else if(tmp_const_numInteractionsNK == 0 && !is.null(tmp_correlationsNK)){
      # If there are no interactions then the object should be null
      stop("The tmp_correlationsNK object should have been passed as null since there are no reported NK correlations, please review")
    }

    # Ok, then we can calculate what the fitness is for a genotype by working out mean the fitness effect of each site (w_i).
    # The fitness value of each site is a function of the state of the site and all other sites upon which it depends.
    # We're defining this function as simply the mean of the independent fitness values
    ## So for each genotype, we pass the states of the sites, which in the initial implementation was recorded by the ID string being the mutant sites
    #### NOTE: As a computational shortcut, I've built the "shape_const_DepbySite_ancestFitness" object so that we know the bySite dependent fitness
    ####	for sites which do not hiold mutations, and thus we update the per site dependent fitness only for those sites which have dependence
    ####	on a mutant site.
    as.vector(sapply(tmp_genotypeList,function(thisGenotype){
      # This finds the independent siteWise fitness value of each site, based on their states as per thisGenotype
      tmp_indFitness <- tmp_focalFitness[,"0"]
      tmp_indFitness[thisGenotype] <- tmp_focalFitness[thisGenotype,"1"]
      # We simply return the mean of the ancestral dependent fitness values, having updated those sites which
      # are mutants.
      tmp_depFitness <- tmp_NK_ancestDep
      tmp_mutCorrelates <- unique(which(apply(tmp_correlationsNK,MARGIN=1,function(theseSites){ any(is.element(thisGenotype,theseSites)) })))
      for(thisSite in union(thisGenotype,tmp_mutCorrelates)){
        tmp_depFitness[thisSite] <- mean(tmp_indFitness[c(thisSite, tmp_correlationsNK[thisSite,])])
      }
      return( mean(tmp_depFitness) )
    }))
  } else if(landscapeModel == "RMF"){
    # in this case we need to know the optimal type (the string passed by tmp_focalFitness) and the weightings of independent
    # and interaction components of the genotype (assumed to be the first and second component of the tmp_wightsRMF  length 2 vector)
    if(length(tmp_weightsRMF) != 1){ stop("The vector of weights passed to RMF fitnessLandscape model is not length 1, please review.") }
    # From Neidhart 2014, we find that the fitness value for a genotype, in their modified RMF landscape model, is given by:
    # F(genotype) = -cD(genotype,optima_genotype) + random_component ; where c is the scaling of the independent fitness contribution
    # and the D() function is a measure of the Hamming distance, between a genotype with the independent fitness optima and the focal genotype.
    # The random component is drawn from a distribution and each genotype has it's own value, meaninging in a binary state genotype there are 2^(func_genomeLength) elements
    # hence we'll be drawing only when we create new genotype's and assign them fitness.
    tmp_randomTerm <- fitnessDist(length(tmp_genotypeList), tmpDistribution = func_distribution, tmpParameters = func_distParameters)
    # If the distribution draws relative fitnesses but needs to be used as selection coefficients...
    if(func_distAsS){
      tmp_randomTerm <- tmp_randomTerm - 1
    }
    ## The math is just the independent term which needs to know the distance to optima, and then a random component
    ## The random component should be of the order of selection coefficients
    as.vector( -tmp_weightsRMF *
                 unlist(lapply(tmp_genotypeList,function(thisGenotype){
                   # The distance between two genotypes can be calculated as the length of one - the number of sites
                   # the second shares in common, + the number of sites the second has in excess
                   tmpOverlap <-is.element(thisGenotype, tmp_optimaRMF)
                   length(tmp_optimaRMF) - length(which(tmpOverlap)) + length(which(! tmpOverlap))
                 })) +
                 tmp_randomTerm )

  }	else if (landscapeModel == "Fixed"){
    # This assumes the user has passed a named vector as the tmp_focalFitness with a complete description
    # of the fitness value for each genotype that could arise in simulations.  The names should be "binaryString" and values the "fitness".
    # The binary string values can be used to return the fitness values
    if(all(is.element(tmpGenotypes,names(tmp_focalFitness)))){
      tmp_focalFitness[tmpGenotypes]
    } else {
      stop("There was a problem trying to use the Fixed fitness landscape model, review input make certain all possible genotypes are declared.")
      Sys.sleep(20)
    }
  }# This is the fitnessVec creation closing

  # Now, if we're using relativeFitness, then we calculate the relative fitness value using our function, which already handles instances
  # of the ancestral fitness being passed as a negative value or of zero or NULL.
  if(relativeFitness){ fitnessVec <- calc_relativeFitness(fitnessVec, func_ancestFit = tmp_ancestralFitness) }
  # We round the fitness vector off to the 4th decimal place ought of simple convenience and aesthetics.
  return( round(fitnessVec,4) )
}



### NOTE: The user must pass at least the minimum correct number of parameters for the distribution chosen
###			and they must be in the same order as the R function expects them
#' This is the function that will call for draws from distributions.
#'
#' @param tmpDraws This is the number of draws sought from the distribution being called
#' @param tmpDistribution This is the character string that represents the implemented distribution you want called.
#' It must be one of: \strong{Fixed, Gamma, Uniform, Normal, Chi2, beta, exp, evd, rweibull, frechet, skewNorm}
#' @param tmpParameters This is the ordered vector of parameters to be passed in order to parameterise the distribution from which
#' you want to draw
#'
#' @return A vector of values with length equal to tmpDraws
#'
#' @examples
#' # This draws from distributions
#' fitnessDist(10, "Uniform", c(0,1))
#' fitnessDist(10, "Normal", c(0,1))
#' fitnessDist(10, "exp", 1)
#'
#' @export
fitnessDist <- function(tmpDraws, tmpDistribution, tmpParameters){
  if(tmpDistribution == "Fixed"){
    # It has been found that if the number of values from which to be sampled is == 1 (meaning we passed two values), then it samples integers...
    # so we force a value to be replciated twice if this is the case
    return( sample(if(length(tmpParameters) == 2){ rep(tmpParameters[2],2) } else { tmpParameters[2: length(tmpParameters)]},
                   tmpDraws, replace = as.logical(tmpParameters[1])) )
  }else if(tmpDistribution == "Gamma"){
    return( rgamma(tmpDraws, shape = tmpParameters[1], rate = tmpParameters[2]) )
  } else if(tmpDistribution == "Uniform") {
    return( runif(tmpDraws, min = tmpParameters[1], max = tmpParameters[2]) )
  } else if(tmpDistribution == "Normal") {
    return( rnorm(tmpDraws, mean = tmpParameters[1], sd = tmpParameters[2]) )
  } else if(tmpDistribution == "Chi2") {
    return( rchisq(tmpDraws, df = tmpParameters[1], ncp = tmpParameters[2]) )
  } else if(tmpDistribution == "beta") {
    return( rbeta(tmpDraws, shape1 = tmpParameters[1] , shape2 = tmpParameters[2]) )
  }  else if(tmpDistribution == "exp") {
    return( rexp(tmpDraws, rate = tmpParameters[1]) )
  }  else if(tmpDistribution == "evd") {
    return( evd::rgev(tmpDraws, loc = tmpParameters[1], scale = tmpParameters[2], shape = tmpParameters[3]) )
  }  else if(tmpDistribution == "rweibull") {
    return( evd::rrweibull(tmpDraws, loc = tmpParameters[1], scale = tmpParameters[2], shape = tmpParameters[3]) )
  }  else if(tmpDistribution == "frechet") {
    return( VGAM::rfrechet(tmpDraws, location = tmpParameters[1], scale = tmpParameters[2], shape = tmpParameters[3]) )
  } else if(tmpDistribution == "skewNorm") {
    # We set the location parameter to 1.05 so the mean is near 1 defined by the omega and alpha
    tmpReturn <- sn::rsn(tmpDraws, xi = tmpParameters[1], omega = tmpParameters[2], alpha = tmpParameters[3], tau = tmpParameters[4])
    # This distribution allows for negative value space, thus we adjust....
    tmpReturn[which(tmpReturn < 0)] <- 0
    return( as.vector(tmpReturn) )
  }
}



# This function will take a focal genotype and then create all the unique genotypes which have one mutation more or less
# which are already not within the reference list.
# NOTE: I'm having this function inherit some pre-definitions, these must be respected!
#' This function searches the nearby mutational space of a focal genotype, identifies which genotypes in that space have not
#' yet been identified, and create new database entries for any new genotypes.
#'
#' @param tmp_focalGenotype This is the focal genotype for which we want to create missing mutational neighbours.
#' @param tmp_focalFitness This is the fitness value of the tmp_focalGenotype.
#' @param maxHamming The maximum number of sites that could be changed by mutation of the tmp_focalGenotype.
#' NOTE: At present I've not made the code work for anything other than a value of 1.  So do not update
#' without updating associated code. where appropriate.
#' @param tmp_landModel This is the character string that defines the fitness landscape model being simulated in this SHAPE run.
#' At present it can be one of: \strong{Additive, Fixed, HoC, NK, RMF}
#' @param tmp_sepString This is a character string used to collapse vectors of characters.
#' @param tmpDirection This is a logical which controls if reversions are allowed (ie: if TRUE sites can revert from mutated to WT)
#' @param tmp_relativeFitness This is a logical which controls if fitness values are to be calculated as relative and no absolute values
#' that would otherwise be calculated via calls to the fitness landscape model.
#' @param tmp_currNeighbours This is an optinal vector that would define the genotype of all neighbours within the 1 step mutational
#' neighbourhood of the tmp_focalGenotype genotype.  If NULL then this vector will be calculated within the function.
#' @param tmp_genCon This is the filepath for the database file that contains the fitness landscape information.
#' @param tmp_tableSplit This is a logical which controls if the tables which report on all genotypes with X mutations should be
#' forced into a single table or it SHAPE is allowed to split them into multiple tables.
#' @param tmp_maxRows The maximum number of rows allowed in a database table before a new table is created.  This has no meaning
#' if tmp_tableSplit is FALSE.
#' @param tmp_genomeLength The length of the genomes, or number of mutable sites/positions, being simulated.
#' @param tmp_distAsS This arugment is passed through to downstream function, but will control if the stochastic portion of
#' fitness effect will be considered as selection coefficients (meaning subtracting 1 from the initially drawn value).
#' @param ... Additional arguments that may get passed to internal functions.
#'
#' @return This invisibly returns NULL, this function is to perform work on databases.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
createGenotypes <- function(tmp_focalGenotype, tmp_focalFitness, maxHamming, tmp_landModel = "HoC", tmp_sepString = getOption("shape_sepString"),
                            tmpDirection = getOption("shape_allow_backMutations"), tmp_relativeFitness = getOption("shape_const_relativeFitness"),
                            tmp_currNeighbours = NULL, tmp_genCon, tmp_tableSplit = getOption("shape_db_splitTables"), tmp_maxRows = getOption("shape_maxRows"),
                            tmp_genomeLength = getOption("shape_genomeLength"), tmp_distAsS = getOption("shape_const_distAsS"), ...){
  # If we haven't previously defined what are the possible neighbours for tmp_focalGenotype, we do so here
  if(is.null(tmp_currNeighbours)){ tmp_currNeighbours <- defineNeighbours(func_tmpGenotype = tmp_focalGenotype, func_tmpDirection = tmpDirection)}
  # We can define aspects of our focalGenotype from the binaryString which ought to have been passed
  tmp_numMuts <- length(strsplit(tmp_focalGenotype, tmp_sepString)[[1]])
  # We find the range of mutants that could be reached in the next step, and bound this by the no-mutant type and the shape_genomeLength
  tmp_mutsRange <-  max(0,(tmp_numMuts - if(tmpDirection){ maxHamming }else{ 0 })): min(tmp_genomeLength ,(tmp_numMuts + (maxHamming)))
  # Unless we allow back mutations AND the maxHamming is > 1, we don't consider the current numMuts as part of the neighbourhood.
  if(!tmpDirection && maxHamming <= 1){
    tmp_mutsRange <- tmp_mutsRange[-which(tmp_mutsRange == tmp_numMuts)]
  }

  # We'll need to query what tables are within the database, both within this function and the tmp_found_neededNeighbours call
  tmp_dbTables <- RSQLite::dbListTables(tmp_genCon)
  #startTime <- proc.time() # This is for reporting on run times but is not needed....
  # We now search for which neighbours already exist within our SQL database
  tmp_found_neededNeighbours <- find_neededNeighbours(tmp_possibleNeighbours = tmp_currNeighbours,
                                                           tmp_focal_numMuts = tmp_numMuts,
                                                           tmpRange_numMuts = tmp_mutsRange,
                                                           tmp_refTables = tmp_dbTables,
                                                           tmp_genCon = tmp_genCon)
  # Now so long as there are neighbours that need to be created, we do so.
  if(length(tmp_found_neededNeighbours) != 0){
    # Now we assign each of these neighbours to a numMuts category for table updating reasons
    # This is done by splitting the strings of all the new neighbours to be assigned, and finding how many mutations based
    # on the length of mutation positions recorded in the compressed nomenclature
    neighbourMuts <- lapply(tmp_mutsRange,function(x){
                         return( tmp_found_neededNeighbours[which(sapply(strsplit(tmp_found_neededNeighbours, tmp_sepString),length) == x)] )
                        })
    names(neighbourMuts) <- sapply(tmp_mutsRange,nameTable, "func_subNaming"=tmp_tableSplit)
    # Now the names of neighbourMuts will be the stringConstant table name for names in our database which reference to a particular numMuts value
    # So we look if there are any tables which contain these names, if not we create a table, if so see the  <else>  section.
    # Also, if there any neighbourMuts list elements without length... meaning there are not mutants needed to be created, we don't pass that instance.
    for(thisTable in names(neighbourMuts)[which(lapply(neighbourMuts,length) > 0)]){
      # We'll be assigning new fitnessLandscape space defined by thisTable's neighbours into some table..
      #startTime <- proc.time() # This is for reporting on run times but is not needed....
      ### NOTE: This call would have a problem generating our ancestral genotype while it's binary string is considered as "",
      ###			I consider this a trivial problem as I expect to always have started runs by defining the no-mutant genotype.
      options("shape_tmpGenoTable" = create_genotypeFrame(tmpID = getOption("shape_nextID"): (getOption("shape_nextID") + length(neighbourMuts[[thisTable]]) - 1),
                                                          tmpStrings = neighbourMuts[[thisTable]],
                                                          tmpFitnesses = fitnessLandscape(tmpGenotypes = neighbourMuts[[thisTable]],
                                                                                               tmp_focalFitness = tmp_focalFitness,
                                                                                               landscapeModel = tmp_landModel,
                                                                                               relativeFitness = tmp_relativeFitness,
                                                                                               func_distAsS = tmp_distAsS)))
      # We now update the value of the getOption("shape_nextID") value for the next iteration.
      options(shape_nextID = getOption("shape_nextID") + length(neighbourMuts[[thisTable]]))
      # This checks if we need to make a new table or are adding to existing tables, it's controlled by an external logical
      # Of whether or not we have a max size of DB table sizes
      if(!any(grepl(thisTable, tmp_dbTables)) || !tmp_tableSplit){
        # In this case we are simply taking the genotypeFrame for all neighbours in thisTable and copying it to the database
        dbWriteTable(tmp_genCon,
                     name = ifelse(tmp_tableSplit, paste(thisTable,1,sep=""), thisTable),
                     value = getOption("shape_tmpGenoTable"),
                     append=TRUE)
      } else {
        # This is trickier, this means we need to identify the number of tables which share this string, and for the last created one
        # we querry if it has the maximum number of rows, if so we create a new one, if not we add to it with these
        tmp_similarTables <- tmp_dbTables[which(grepl(thisTable, tmp_dbTables))]
        # Now to find which is the oldest we look for the highest value in the 3rd separated position
        tmpOldest <- tmp_similarTables[which.max(sapply(tmp_similarTables,function(x){ as.numeric(strsplit(x,tmp_sepString)[[1]][3]) }))]
        # We querry the number of rows in this table and if there are fewer than our shape_maxRows we'll add, otherwise we create a new one
        if(nrow(dbGetQuery(tmp_genCon,paste("SELECT genotypeID FROM ", tmpOldest,sep=""))) < tmp_maxRows){
          # This means we can simply add our new table to the existing one
          dbWriteTable(tmp_genCon,
                       name = tmpOldest,
                       value = getOption("shape_tmpGenoTable"),
                       append=TRUE)
        } else {
          # We write a new table which simply takes the index (third indexed sepString position) and advance to the next one....
          dbWriteTable(tmp_genCon,
                       name = paste(paste(strsplit(tmpOldest,tmp_sepString)[[1]][-3],collapse=tmp_sepString),as.numeric(strsplit(tmpOldest,tmp_sepString)[[1]][3])+1,sep=tmp_sepString),
                       value = getOption("shape_tmpGenoTable"))
        }

      }# This closes out  if we needed to make a new table or add to an existing one
    } # This closes out going through all the table types that might need to be created as a result of numMuts in genotypes

    # Now we go to the ancestral genotype's table and we update that it has been explored
    # This requires we find it's table, so we search for all tables in it's mut range
    tmp_neededTables <- tmp_dbTables[which(grepl(nameTable(tmp_numMuts, func_subNaming=tmp_tableSplit), tmp_dbTables))]
    # If the focalGenotype is anything other than our absolute ancestor we need to find the tables and use it's binaryString
    # However, for the ancestor we don't do this since it has no binaryString value.... we'd have needed it's genotypeID
    tmp_ancestFind <- 0
    tmp_ancestTable <- tmp_dbTables[which(grepl(nameTable(0, func_subNaming=tmp_tableSplit), tmp_dbTables))]
    names(tmp_ancestFind) <- tmp_ancestTable
    if(tmp_focalGenotype != ""){
      # We find the proper table and then the genotypeID of our ancestor
      tmp_ancestFind <- sapply(tmp_neededTables, function(thisTable){
        as.vector(unlist(dbGetQuery(tmp_genCon, paste("SELECT genotypeID FROM ",thisTable,' WHERE binaryString = "', tmp_focalGenotype,'"',sep=""))))
      },simplify=FALSE)
      # There should be a single returned non NULL value
      tmp_ancestTable <- names(tmp_ancestFind)[which(sapply(tmp_ancestFind,length) == 1)]
    }
    dbExecute(tmp_genCon,
              paste("UPDATE ", tmp_ancestTable ,' SET isExplored = 1 WHERE genotypeID = ', tmp_ancestFind[tmp_ancestTable],sep=""),
              synchronous = NULL)

  } # This closes out if we had any neighbours that needed to be created.

  # We don't return anything, we've performed work....
  invisible( NULL )
}



# This function will help us investigate the local neighbourhood of a focal genotype and see if we need to create new genotype information
# NOTE: It assumes the tmp_focalGenotype being passed to it is either a character string, or a numeric vector
#' This function querries if a suite of genotypes exist within the fitness landscape database.
#'
#' @param tmp_possibleNeighbours This is a vector of all possible mutants that we're trying to querry within the fitness landscape database.
#' @param tmp_focal_numMuts This is the number of mutations in the focal genotype, it controls - along with other parameters -
#' what tables of the fitness landscape database are querried.
#' @param tmp_refTables This is the a vector of named tables that exist within the fitness landscape.  It can not be passed
#' in which case the database at tmp_genCon is querried for this information.
#' @param maxHamming The maximum number of sites that could be changed by mutation of the tmp_focalGenotype.
#' @param tmp_tableSplit This is a logical which controls if the tables which report on all genotypes with X mutations should be
#' forced into a single table or it SHAPE is allowed to split them into multiple tables.
#' @param tmp_genomeLength The length of the genomes, or number of mutable sites/positions, being simulated.
#' @param tmpDirection This is a logical which controls if reversions are allowed (ie: if TRUE sites can revert from mutated to WT)
#' @param tmpRange_numMuts This is the range of number of mutations which a mutant neighbour may posses.  If not supplied that will
#' be calculated in line via other parameters passed to the function.
#' @param tmp_genCon This is the filepath for the database file that contains the fitness landscape information.
#'
#' @return A vector of the genotypes that need to be created as they've not yet been defined within the fitness landscape.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
find_neededNeighbours <- function(tmp_possibleNeighbours, tmp_focal_numMuts, tmp_refTables,
                                  maxHamming = getOption("shape_max_numMutations"),
                                  tmp_tableSplit = getOption("shape_db_splitTables"), tmp_genomeLength = getOption("shape_genomeLength"),
                                  tmpDirection = getOption("shape_allow_backMutations"), tmpRange_numMuts = NULL, tmp_genCon){
  #startTime <- proc.time() # This is for reporting on run times but is not needed....
  # If we haven't been passed the vector of all tables in our reference database, then we define it here
  if(is.null(tmp_refTables)){ tmp_refTables <- RSQLite::dbListTables(tmp_genCon) }
  # Now we find which of these are not elements of binaryStrings already within our database, we look for possible neighbours that are not
  # found within the database (called using the select function of "dplyr"
  # This is started by querrying for binaryStrings from among meaningfull tables in our database , we define which meaningful tables exist:
  # If we've passed a range for the number of mutants then we skip this, otherwise we redefine the NULL parameter
  tmp_neededTables <- NULL
  if(is.null(tmpRange_numMuts)){
    tmpRange_numMuts <- max(0,(tmp_focal_numMuts - if(tmpDirection){ maxHamming }else{ 0 })): min(tmp_genomeLength,(tmp_focal_numMuts + maxHamming))
    tmp_neededTables <- tmp_refTables[unlist(lapply(tmpRange_numMuts[-which(tmpRange_numMuts == tmp_focal_numMuts)],function(x){ which(grepl(nameTable(x, func_subNaming=tmp_tableSplit), tmp_refTables)) }))]
  } else {
    tmp_neededTables <- tmp_refTables[unlist(lapply(tmpRange_numMuts,function(x){ which(grepl(nameTable(x, func_subNaming=tmp_tableSplit), tmp_refTables)) }))]
  }
  # We now ask if there are any table from which we can compare for Neighbours, if there aren't (length == 0) we needn't do anything as all neighbours need to be made.
  if(length(tmp_neededTables) == 0){
    # Well in this case we can return all possible neighbours as ones that need to be created.
    return( unlist(tmp_possibleNeighbours) )
  } else {
    # Then we querry our database, among the meaningfull tables, and take the possible neighbours that are not in the database yet.
    # We grab as a vector all the binary strings which are within those meaningfull table
    tmp_binaryStrings <- gsub("[[:space:]]","",paste("\'", tmp_possibleNeighbours,"\'",collapse=','))

    # We now look for which of these possible neighbours exist within our database and return the needed ones
    return ( setdiff(tmp_possibleNeighbours,
                     unlist(dbGetQuery(tmp_genCon,
                                       paste('SELECT binaryString FROM ',
                                             tmp_neededTables,
                                             ' WHERE binaryString IN (',
                                             tmp_binaryStrings,')',
                                             collapse=" UNION "))
                     )) )
  }
}



# This is a function which informs us of which binary strings are the nearest neighbours of a focal genotype
# We expect that at least a binaryString single length character vector is passed, from this we can work out number of mutations,
# Otherwise if a vector is passed and/or number of mutations therein we use that directly,
# lastly we need to be informed if we're allowing back mutations, technically the user can define the maxHamming distance, but we've only coded it for 1...
#' The function will identify the binary string of all possible neighbours to a focal genotype.  It is important when querrying the fitness landscape.
#'
#' @param func_tmpGenotype This is the binary string of the focal genotype for which we want to define possible neighbours.
#' @param func_tmpDirection This is a logical which controls if reversions are allowed (ie: if TRUE sites can revert from mutated to WT)
#' @param func_maxHamming The maximum number of sites that could be changed by mutation of the tmp_focalGenotype.
#' NOTE: At present I've not made the code work for anything other than a value of 1.  So do not update
#' without updating associated code, where appropriate.
#' @param func_sepString This is a character string used to collapse vectors of characters.
#' @param func_genomeLength The length of the genomes, or number of mutable sites/positions, being simulated.
#'
#' @return Vector of all the genotypes in the neighbouring mutational space accessible within 1 mutation event
#'
#' @examples
#' # If you had some individuals with a genome length of 10 sites, and an
#' # individual with no mutations, as well as one with a single mutation at
#' # position 7, each had a mutant.  This would define the possible one step
#' # mutational neighbours.  I also allow back mutations
#' defineNeighbours(c(""), func_tmpDirection = FALSE, func_maxHamming = 1,
#'                  func_sepString = "_", func_genomeLength = 10)
#' defineNeighbours(c("7"), func_tmpDirection = FALSE, func_maxHamming = 1,
#'                  func_sepString = "_", func_genomeLength = 10)
#' #' # Same idea, but if we allow back-mutations (ie: reversions)
#' defineNeighbours(c("7"), func_tmpDirection = TRUE, func_maxHamming = 1,
#'                  func_sepString = "_", func_genomeLength = 10)
#'
#' @export
defineNeighbours <- function(func_tmpGenotype, func_tmpDirection, func_maxHamming = getOption("shape_max_numMutations"),
                             func_sepString = getOption("shape_sepString"), func_genomeLength = getOption("shape_genomeLength")){
  # We expect to be passed the func_tmpGenotype as the binaryString which is the position of each mutation with a separator
  func_tmpGenotype <- as.numeric(strsplit(func_tmpGenotype, func_sepString)[[1]])
  # Now, if there are no mutations in this genome, we can only add mutations, so we define each position in the genome as a possible mutant type
  return( if(length(func_tmpGenotype) == 0){
            # This instance should only occur when the func_tmpGenotype is the WT, and thus the single step mutants involve
            # a mutation at each possible position in the genome.
            as.character(seq(1, func_genomeLength))
            ################################################################################################################################
            #### THESE TWO "else if" CALLS are hard coded for a maxHamming of 1, these will require meaningful updates if this changes. ####
            ################################################################################################################################
            # Otherwise we need to consider adding mutations to each other position in the forward direction if one direction only
          } else if(length(func_tmpGenotype) > 0 && !func_tmpDirection){
            unlist(lapply(seq(1, func_genomeLength)[-func_tmpGenotype],function(thisPos){
              paste(c(func_tmpGenotype,thisPos)[order(c(func_tmpGenotype,thisPos))],collapse= func_sepString)
            }))
            # This is the case where we allow back mutations, meaning we permit the addition of mutations and their removal
          } else if(length(func_tmpGenotype) > 0 && func_tmpDirection){
            # So we'll return the concatenation of all those values where we add a mutation and remove one.
            c(unlist(lapply(seq(1, func_genomeLength)[-func_tmpGenotype],function(thisPos){
              paste(c(func_tmpGenotype,thisPos)[order(c(func_tmpGenotype,thisPos))],collapse= func_sepString)
            })),
            unlist(lapply(1:length(func_tmpGenotype),function(thisPos){
              paste(c(func_tmpGenotype[-thisPos])[order(c(func_tmpGenotype[-thisPos]))],collapse= func_sepString)
            })) )
          } ) # This closes the return call
}


# This is a function for creating genotype_refDatabse index data.frames
#' This is a convenience function to ensure that we have a standard shaped data.frame.  It is used to initiate
#' a new table for the fitness landscape.
#'
#' @param tmpID A numeric vector of the unqiue identifiers for genotypes
#' @param tmpStrings A vector of the character strings that represent the binary string of genotypes
#' @param tmpFitnesses A vector of the numeric fitness values to be input
#'
#' @return A 4 column data frame with column names of genotypeID, binaryString, fitness, isExplored
#'
#' @examples
#' # This is just a convenience function for outputting vectors in a data.frame with
#' # standard named columns.
#' create_genotypeFrame(c(1,10,50),c("1","1_7","6_12"),c(1,0.25,1.57))
#' @export
create_genotypeFrame <- function(tmpID, tmpStrings,tmpFitnesses){
  if(length(tmpStrings) != length(tmpFitnesses)){
    stop("There was a problem trying to create a genotype data frame, as a result of the number of genotypes and fitnesses passed do not match")
  }
  return( data.frame("genotypeID"=tmpID,"binaryString"= tmpStrings, "fitness"= tmpFitnesses, "isExplored"= 0, stringsAsFactors=FALSE) )
}



# This is a function for the naming of tables to be copied to the SQL database.  NOTE: The func_tmpIndex value will be reduced to a single string element
# We also incldue a logical to suppress the use of sub-inxeding when creating names
#' This is a standardising function which allows SHAPE to programatiically name tables for the fitness landscape OR split
#' a named table and extract the embedded information from its naming.
#'
#' @param func_tmpMutations Integer value(s) for the number of mutations to be expected in mutants stored within the named tables.
#' @param func_tmpIndex An optinal element that will be used to insert a unique vector ID
#' @param func_baseString This is the standard prefix character string used in table naming.
#' @param func_sepString This is a character string used to collapse vectors of characters.
#' @param func_splitName A logical toggle to control if this function is splitting a named table or not.  So, FALSE (default)
#' means we're creating a table name whereas TRUE is splitting a named table into it's parts.
#' @param func_subNaming This is a logical which controls if the tables which report on all genotypes with X mutations should be
#' forced into a single table or it SHAPE is allowed to split them into multiple tables.
#'
#' @return If func_splitName is TRUE, then a vector of table names is returned, it would be best practice to not assume recycling of passed
#' elements and so pass equally lengthed vectors as input.  If FALSE, we split the table and return the data detailing the number of mutations
#' which ought to be present for genotypes stored in the named table.
#'
#' @examples
#' # This creates a table name in a standard way, it can also split table names to extract info.
#' defineSHAPE()
#' nameTable(2,1,"myTest","_",FALSE,FALSE)
#' nameTable("myTest_2",func_splitName = TRUE)
#'
#' @export
nameTable <- function(func_tmpMutations, func_tmpIndex = NULL, func_baseString = getOption("shape_string_tableNames"),
                      func_sepString = getOption("shape_sepString"), func_splitName = FALSE, func_subNaming = getOption("shape_db_splitTables")){
  # As a afety I've included an option to split rather than build table names to recover the numMuts stored therein
  if(func_splitName){
    return( as.vector(unlist(lapply(func_tmpMutations, function(x) { strsplit(x, func_sepString)[[1]][2] }))) )
  }
  # Now if we're building names we need to know if we are indexing or not....
  if(func_subNaming){
    if(!is.null(func_tmpIndex)){
      return( paste(func_baseString,func_tmpMutations,paste(func_tmpIndex,collapse=func_sepString),sep=func_sepString) )
    } else {
      # We've added this trailing func_sepString element so that numeric indexces can later be better tracked within names using grepl calls with fixed = TRUE
      return( paste(func_baseString,func_tmpMutations,"",sep=func_sepString) )
    }
    # This ignores the presence of the func_tmpIndex value
  } else {
    return( paste(func_baseString,func_tmpMutations,sep=func_sepString) )
  }
}



# This is to be used for the names of our step report object tables
#' This is a standardising function which allows SHAPE to programatiically name tables for the step-wise record OR split
#' a named table and extract the embedded information from its naming.
#'
#' @param func_Index Integer value(s) for the step of a SHAPE run which will be recorded by this table
#' @param func_sepString This is a character string used to collapse vectors of characters.
#' @param funcSplit A logical toggle to control if this function is splitting a named table or not.  So, FALSE (default)
#' means we're creating a table name whereas TRUE is splitting a named table into it's parts.
#'
#' @return If funcSplit is TRUE, then a vector of table names is returned.  If FALSE, we split the table and return the
#' data detailing the step number being recorded on the named table.
#'
#' @examples
#' # This creates a table name in a standard way, it can also split table names to extract info.
#' defineSHAPE()
#' nameTable_step(2,FALSE)
#' nameTable_step("Step_2",TRUE)
#'
#' @export
nameTable_step <- function(func_Index, funcSplit = FALSE, func_sepString = getOption("shape_sepString")){
  # This checks if we're asking for the name to be split or not
  if(funcSplit){
    # We return only the second piece of information given the setup of this function's naming practice.
    return( unlist(lapply(strsplit(func_Index,func_sepString),function(x){ x[2] })) )
  } else {
    return( paste("Step",func_Index,sep=func_sepString) )
  }
}


# This is to be used for the names of our neighbourhood object tables
#' This is a standardising function which allows SHAPE to programatiically name tables for the neighbourhood record OR split
#' a named table and extract the embedded information from its naming.
#'
#' @param func_Index Integer value(s) for the unique genotype ID whose neighbourhood which will be recorded by the named table
#' @param func_sepString This is a character string used to collapse vectors of characters.
#' @param funcSplit A logical toggle to control if this function is splitting a named table or not.  So, FALSE (default)
#' means we're creating a table name whereas TRUE is splitting a named table into it's parts.
#'
#' @return If funcSplit is TRUE, then a vector of table names is returned.  If FALSE, we split the table and return the
#' data detailing the genotype ID whose neighbourhood is being recorded on the named table.
#'
#' @examples
#' # This creates a table name in a standard way, it can also split table names to extract info.
#' defineSHAPE()
#' nameTable_neighbourhood(2,FALSE)
#' nameTable_neighbourhood("Step_2",TRUE)
#'
#' @export
nameTable_neighbourhood <- function(func_Index, funcSplit = FALSE, func_sepString = getOption("shape_sepString")){
  # This checks if we're asking for the name to be split or not
  if(funcSplit){
    # We return only the second piece of information given the setup of this function's naming practice.
    return( unlist(lapply(strsplit(func_Index,func_sepString),function(x){ x[2] })) )
  } else {
    return( paste("genotypeID",func_Index,sep=func_sepString) )
  }
}


# This is a little cheater function to close and re-open our database connection after each write
# We pass a filename for the connection we want.
#' This is a convenience function to refresh connections to database files.
#'
#' @param func_conName The filepath to which an SQLite connection is sought.
#' @param func_existingCon If any value other than NULL, then any existing connection is first dropped
#' prior to attempting to form a connection to the func_conName filepath.
#' @param func_type This should be a character string of either \strong{connect}, in which case a
#' connection is made/refreshed to the filepath in func_conName", or any other value will cause disconnection
#'
#' @return An SQLite connection object to an SQLite database.
#'
#' @examples
#' # This function can be called to set, resset SQL connections
#' fileName_testCon <- paste(tempdir(),"/testCon.sqlite",sep="")
#' testCon <- reset_shapeDB(fileName_testCon)
#' reset_shapeDB(testCon, func_type = "disconnect")
#'
#' @export
reset_shapeDB <- function(func_conName, func_existingCon = NULL, func_type = "connect"){
  # Now we either open or close the connection for the conName passed
  if(func_type == "connect"){
    # If there is an existing connection we'll close it
    if(!is.null(func_existingCon)){ RSQLite::dbDisconnect(func_existingCon) }
    return( RSQLite::dbConnect(RSQLite::SQLite(), func_conName) )
  } else {
    return( RSQLite::dbDisconnect(func_conName) )
  }
}



# This is a function which will allow us to generate matrices easily for reporting on the population(s) in a step
#' This is a convenience function to ensure that our population demographics are stored in a data frame
#' and exists because R's standard functions can collapse single row frames to named vectors.  It requires that
#' all passed vectors be of the same length
#'
#' @param func_numMuts This is a vector of the number of mutations held within each tracked genotype.
#' @param func_genotypeID This is a vector of the unique genotype ID for each tracked population in the community.
#' @param func_popSizes This is a vector of the number of individuals for each population of genotypes in the community.
#' @param func_fitnesses This is a vector of the fitness for each genotpe being tracked.
#' @param func_births This is a vector of the number of births produced by each population in this time step.
#' @param func_deaths This is a vector of the number of deaths in each population in this time step.
#' @param func_mutants This is a vector of the number of mutants produced by each population in this time step.
#' @param func_progenitor This is a vector of character strings expressing any progenitor genotypes which generated a mutant
#' that fed into each genotype's population in this time step.
#' @param func_reportMat_colnames DO NOT MODIFY - This is the vector of character strings to be assigned as the column names.
#'
#' @return A data frame with columns named as per func_reportMat_colnames.
#'
#' @examples
#' # This returns a data.frame with a standard format
#' defineSHAPE()
#' reportPopulations(1:3,2:4,c(10,50,100),rep(1,3),
#'                  rep(0,3),c(10,10,10),c(1,2,0),c("","0_->_1","2"))
#' @export
reportPopulations <- function(func_numMuts, func_genotypeID, func_popSizes, func_fitnesses, func_births, func_deaths, func_mutants, func_progenitor,
                              func_reportMat_colnames = getOption("shape_reportMat_colnames")){
  # We require that the user send us each piece of information being the same length, otherwise we complain.  Similarly if there are less
  # elements than epected (as per func_reportMat_colnames) we complain as well
  length_allInputs <- c(length(func_numMuts), length(func_genotypeID), length(func_popSizes), length(func_fitnesses),
                        length(func_births), length(func_deaths), length(func_mutants), length(func_progenitor))
  if( all(sapply(length_allInputs[-1],function(tmpLengths){ length_allInputs[1] == tmpLengths  })) ){
    if(length(length_allInputs) != length(func_reportMat_colnames)){
      stop("Trying to create population matrix with insufficient number of elements passed, please review")
      return( NULL )
    } else {
      # If so then we'll use these to build a matrix, otherwise we'll complain about the information passed along
      tmpReturn <- data.frame(func_numMuts, func_genotypeID, func_popSizes, func_fitnesses, func_births, func_deaths, func_mutants, func_progenitor,stringsAsFactors=FALSE)
      dimnames(tmpReturn) <- list(rep(func_genotypeID,max(length_allInputs)/length(func_genotypeID)),
                                  func_reportMat_colnames)
      return( tmpReturn )
    }
  } else {
    stop("Trying to create population matrix with vectors of information that are not all the same size, please review")
    return( NULL )
  }
}



#' This is a function to search our mutational database and then find the binary string of the genotypeID passed.
#' This function is more efficient when the number of mutations for each genotypeID be passed as this helps reduce
#' the tables of the mutational space that are searched.  This matters when large genotypes are simulated.
#'
#' @param func_genotypeID This is a vector of the unique genotype ID for each tracked population in the community.
#' @param func_numMuts This is a vector of the number of mutations held within each tracked genotype.
#' @param func_subNaming This is a logical which controls if the tables which report on all genotypes with X mutations should be
#' forced into a single table or it SHAPE is allowed to split them into multiple tables.
#' @param func_landscapeCon This is the filepath to an SQLite database storing information for the complete explored and neighbouring fitness landscape of a SHAPE run.
#'
#' @return This returns a vector of character strings that represent the binary strings of the genotypes
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
retrieve_binaryString <- function(func_genotypeID, func_numMuts = NULL, func_subNaming,
                                  func_landscapeCon){
  # We find the needed tables by using the func_numMuts value passed
  tmp_refTables <- RSQLite::dbListTables(func_landscapeCon)
  # If we've been informed about the number of mutations this ID has, then we subset our tables
  if(!is.null(func_numMuts)){
    tmp_refTables <- tmp_refTables[which(grepl(nameTable(func_numMuts, func_subNaming = func_subNaming), tmp_refTables))]
  }
  # We now query all the table(s) for the information
  func_tmpReturn <- dbGetQuery(func_landscapeCon, paste("SELECT binaryString,isExplored FROM ",
                                                        tmp_refTables,
                                                        ' WHERE genotypeID IN(',
                                                        paste(func_genotypeID,collapse=","),
                                                        ')',
                                                        collapse=" UNION "))
  # If nothing was found then we'll have a nrow == 0 object, in which case we simply search all tables in the DB
  # I don't use recursion to avoid being caught in an infinite loop, and instead put a crash out marker
  if(nrow(func_tmpReturn) == 0){
    func_tmpReturn <- dbGetQuery(func_landscapeCon, paste("SELECT binaryString,isExplored FROM ",
                                                          RSQLite::dbListTables(func_landscapeCon),
                                                          ' WHERE genotypeID IN(',
                                                          paste(func_genotypeID,collapse=","),
                                                          ')',
                                                          collapse=" UNION "))
  }
  # Last check
  if(nrow(func_tmpReturn) == 0){
    stopError(paste("Problem finding the binaryString information for: ",func_genotypeID," please review",sep=""))
  } else {
    return( func_tmpReturn )
  }
}



# This function will compute the size of disturbance which occurs in a population, and thus the number of steps until its recovery
# The argument for maxSize exists but I am uncertain that the code's postAnalysis script will handle if this value changes.... it may, haven't thought about it very hard.
# The cappedGrowth logical is to define if the maxSize value is meant to be the carrying capacity or some initial value of the population.
#' This function is used to calculate the effect size and timing of the next stochastic population disturbance in a SHAPE run.
#'
#' @param func_distFactor This is the expected effect size of the disturbance, it should be a named vector with elements
#' \strong{factor, random} which are each used as per the func_distType
#' @param func_growthType This is the growth model of the SHAPE run
#' @param func_distType This is the type of disturbance to be simulated.  Currently I've implemented
#' \strong{bottleneck, random} options for constant bottlenecks or normally distributed random effect sizes
#' @param func_growthRate This is the basal growth rate of the SHAPE run
#' @param func_popSize This is a vector of the number of individuals in each of the populations
#' @param func_focalSize This only matters if the growth model is exponential in which case the disturbance
#' is always such that the community size is reduced to the func_focalSize value
#' @param func_manualGenerations If not NULL, it will be rounded to an integer value and taken as the manually
#' controlled number of generations between disturbances.  Otherwise, the disturbance factor and growth rate are
#' used to estimate the number of steps required for a community with relative fitness 1 to rebound.
#' @param func_stepDivs This is the value that controls what proportion of a standard biological "generation" is
#' simulated in each step of a SHAPE run.
#'
#' @return A named vector with three elements describing the simulated reduction factor of populations,
#' the number of individuals lost, and the number of steps estimated until the next disturbance.
#'
#' @examples
#' # This calculates the information for the next planned stochastic disturbance event.
#' # Consider a situation where there is a disturbance reducing populations 100 fold,
#' # and it occurs either in a proscriptive number of steps, or we calculate it based
#' # on recovery time as per the growth rate and growth model parameters.
#' compute_distGrowth("bottleneck","exponential","bottleneck",
#'                      2,1e4,1e2,5,1)
#' compute_distGrowth("bottleneck","exponential","bottleneck",
#'                    2,1e4,1e2,NULL,1)
#' # If growth is constant or Poisson, then disturbances are effectively supressed
#' compute_distGrowth("bottleneck","poisson","bottleneck",
#'                    2,1e4,1e2,NULL,1)
#'
#' @export
compute_distGrowth <- function(func_distFactor, func_growthType, func_distType, func_growthRate,
                               func_popSize, func_focalSize, func_manualGenerations = NULL, func_stepDivs){
  # If the growth form is constant growth them we simply return that there will never be a disturbance
  # If the user has not set a number of generations prior to disturbance we set this as a large number.
  if(is.element(func_growthType, c("poisson","constant"))){
    return( c("factor"=0,
              "popLost"= 0,
              "stepReq"= if(is.null(func_manualGenerations )){
                               1e12
                             } else {
                               ceiling(func_manualGenerations/func_stepDivs)
                             }) )
  }

  # We calculate the factor of population loss based on the growth type, for exponential or poisson this is simply a resseting of the population
  # to it's initial size
  func_tmpFactor <- as.vector(if(func_growthType == "logistic"){
    if(func_distType == "bottleneck"){
      func_distFactor["factor"]
    } else if(func_distType == "random") {
      rnorm(1,mean= func_distFactor["factor"],sd= func_distFactor["random"])
    } else {
      stop(paste('There was a problem with the definition of <func_distType> = ', func_distType, ' please review',sep=""))
    }
  } else if (func_growthType == "exponential") {
    (sum(func_popSize)/func_focalSize)
  })
  # Since, in theory, the dilution factor could be a value smaller than 1 (due to distribution size or growth form), this would mean that
  # a proposed dilution would be of a form that increases the population size.  Thus if func_tmpFactor < 1 we set it to 1
  if(func_tmpFactor < 1){ func_tmpFactor <- 1 }

  # Now knowing the factor of loss, we calculate the number of steps in which growth will occur based on the division type
  # We also adjust this by the number of divisions with each step of growth.  We take the ceiling since we cannot run partial growth steps....
  func_growthSteps <- if(is.null(func_manualGenerations)){
                        # If there aren't a mannual number of generations defined, then we take the log, base of growth rate,
                        # of the disturbance factor.  However, if there was no disturbance, because perhaps this is an initialisation
                        # step then we will set the this value to a default of 100 fold dilution to replicate serial passaging
                        log(if(func_tmpFactor <= 0){ 100 } else { func_tmpFactor }, base = func_growthRate)
                      } else {
                        round(func_manualGenerations,0)
                      }
  func_growthSteps <- ceiling(func_growthSteps/func_stepDivs)
  # We return what we;ve calculated
  return( c("factor" = func_tmpFactor,
            "popLost" = 0,
            "stepReq" = as.vector(func_growthSteps)) )
}



# This is a function which calculates the amount of loss suffered by a population based on a dilution factor suffered by the population
# We return the lesser of the initial size OR the drawn value - as we cannot lose more than the initial number
#' This function actually calculates the stochastic loss to populations.
#'
#' @param func_inPopulation This is a vector of the number of individuals in the populations within the community.
#' @param func_dilutionFactor This is expected proportion of the current population sizes that should remain.
#'
#' @return A vector of the resultant population sizes remaining.
#'
#' @examples
#' # A vector of population sizes is randomly sampled to be around the product of size and factor
#' replicate(5,lossSampling(c(1e4,2e4,3e4),0.01))
#'
#' @export
lossSampling <- function(func_inPopulation, func_dilutionFactor){
  # Here, the loss experienced is quite simply based on the the size of each lineage
  return( apply(cbind(func_inPopulation, func_dilutionFactor),MARGIN=1,function(thisPop){
               min(thisPop[1],rpois(1, lambda = thisPop[1] * thisPop[2]))
             }) )
}



# This function will allow me to dynamically build unique names for jobs performed with SHAPE, sepString is presumed to be in the evironment
# When splitting the jobID, setID, and repID values should be passed as logicals to control if we expect those to exist in what was passed.
#' This function is used to build or split character string to be used for naming batches of SHAPE runs.
#'
#' @param funcBase If building names this is the basal string element prefixing the name.  If splitting, it is the vector of names to be split.
#' @param func_setID If building names, a vector of the unique set IDs to be named, otherwise a logical of whether or not the batch naming structure includes sets
#' @param func_jobID If building names, a vector of the unique job IDs to be named, otherwise a logical of whether or not the batch naming structure includes jobs
#' @param func_repID If building names, a vector of the unique replicate IDs to be named, otherwise a logical of whether or not the batch naming structure includes replicates
#' @param funcSplit Logical toggle TRUE if splitting names, FALSE to build string characters
#' @param func_sepString This is the standard string separator for the SHAPE run
#'
#' @return Either a vector of character strings for the created batch names, or a matrix with the decomposed
#' elements of the split batch name strings
#'
#' @examples
#' # This simply produces or splits a standard named string.
#' name_batchString("myTest",1,9,3,FALSE,"_")
#' name_batchString("myTest_1_9_3",TRUE,TRUE,TRUE,TRUE,"_")
#'
#' @export
name_batchString <- function(funcBase, func_setID = NULL, func_jobID = NULL, func_repID = NULL,
                             funcSplit = FALSE, func_sepString = getOption("shape_sepString")){
  # Only the base need be defined, if the others are not supplied then the base is assumed to be the name
  # The split option allows a user to extract information from the name if passed as the base
  if(!funcSplit){
    # I return a vector of name(s) built from the information passed, warning recycling may happen I don't control it.
    return( unname(apply(cbind(funcBase, func_setID, func_jobID, func_repID),MARGIN=1,function(thisInfo){
                          paste(paste(thisInfo,collapse=func_sepString,sep=""),sep="")
                        })) )
  } else {
    # We check that all the c(func_jobID, func_setID, func_repID) are logicals
    if(!all(is.logical(c(func_setID, func_jobID, func_repID)))){
      stop(print("Problem splitting the batchString of ",paste(funcBase,collapse=" ")," missing logicals.",sep=""))
    }
    # I return a matrix, possibly with only a single row, of the split information.
    func_tmpReturn <-  matrix(sapply(strsplit(funcBase,func_sepString),function(thisSplit){
                               # The tmpReturn is built by building a vector based on the presence of ID logicals
                               tmp_returnString <- thisSplit[1]
                               thisSplit <- thisSplit[-1]
                               for(func_tmpIndex in c(func_setID, func_jobID, func_repID)[which(c(func_setID, func_jobID, func_repID))]){
                                 if(func_tmpIndex){
                                   tmp_returnString <- c(tmp_returnString,thisSplit[1])
                                   thisSplit <- thisSplit[-1]
                                 } else {
                                   tmp_returnString <- c(tmp_returnString,NA)
                                 }
                               }
                               return( tmp_returnString )
                             }),
                             ncol = length(funcBase))
    dimnames(func_tmpReturn) <- list(c("base",
                                       if(!is.null(func_setID)){if(as.logical(func_setID)){"setID"}else{NULL}}else{NULL},
                                       if(!is.null(func_jobID)){if(as.logical(func_jobID)){"jobID"}else{NULL}}else{NULL},
                                       if(!is.null(func_repID)){if(as.logical(func_repID)){"repID"}else{NULL}}else{NULL}),
                                     funcBase)
    return( func_tmpReturn )
  }
}



#' This is a function to calculate the relative fitness for a vector of fitnesses.  As a frame of reference
#' it can use either an ancestral fitness value or the mean fitness of the passed vector.  If the frame of
#' reference is a value of zero - OR - the func_absDistance is set to TRUE then instead the vector is centered
#' around a value of 1 where negative values will be set to zero.
#'
#' @param func_fitVector a numeric vector of values to be interpreted as fitnesses
#' @param func_ancestFit An optional single numeric value to be used as a frame of reference
#' for calculating relative fitness.
#' @param func_weights An optional vector of weights to be used for calculating relative fitness
#' as an absolute distance from the mean of the func_fitVector vector.
#' @param func_absDistance A logical toggle to override if relative fitnesses are to be calculated as the absolute
#' distance from 1.  Will be overrode if either the mean of func_fitVector or func_ancestFit are zero.
#'
#' @return A vector of relative fitness values of length equal to the input vector.
#'
#' @examples
#' # This calculates relative fitness values either based on the mean of the community or
#' # based on an ancestral fitness value.
#' defineSHAPE()
#' calc_relativeFitness(c(0.9,1,1.1))
#' calc_relativeFitness(c(0.9,1,1.1),func_ancestFit = 0)
#' calc_relativeFitness(c(0.9,1,1.1),func_ancestFit = 1)
#' calc_relativeFitness(c(0.95,1,1.1))
#'
#' @export
calc_relativeFitness <- function(func_fitVector, func_ancestFit = NULL, func_weights = NULL,
                                 func_absDistance = (getOption("shape_simModel") == "RMF")){
  # If we've passed an ancestral fitness and it is zero ,then we're calculating fitness values as selection coefficients
  # thus we want to use a distance measure for our relative fitness.
  if(!is.null(func_ancestFit)){
    if(func_ancestFit == 0){
      # This is irrespective of the shape_simModel passed as it's based on what the ancestral fitness value was defined to be
      func_absDistance <- TRUE
    }
  }
  # There are some special circumstances to be handled.  Simplest case, is that we've been passed along
  # and ancestral fitness value, that is not zero.  If the ancestral fitness value is zero, we simply ignore it
  # and use our other mehtods of centering values and calculating relative fitness.  Also if we pass an ancestral fitness
  # but it is zero, then we're calculating distributions to represent (s), so we use distance.
  if( (!is.null(func_ancestFit) && func_absDistance) ){
    # We return the distance from the ancestor to a minimum of zero
    return( unlist(lapply(1 + (func_fitVector - func_ancestFit),function(func_thisFitness){ max(0, func_thisFitness) })) )
  } else if(!is.null(func_ancestFit) && func_ancestFit != 0){
    # Simplest case is the ancestral fitness is not zero. This is simple, we just divide the fitness vector by the ancestral fitness value
    # But multiply by it's sign so that the "direction" of fitness increase is preserved in the event of a negative ancestor.  We also consider
    # that negative fitness values
    return( unlist(lapply(sign(func_ancestFit) *func_fitVector/func_ancestFit, function(func_thisFitness){ max(0, func_thisFitness) })) )
    # If the ancestral fitness is zero, then we won't be using it, thus we pass to our other methods.
  } else {
    # This means we use the fitness vector that was passed will be centered around 1, and negative numbers are treated as zero
    # While this tep would account for the magnitude of values, it's not appropriate in so much as the user ought to be
    # defining DFE parameters such that the values are reasonably selection coefficients or relatiev fitness values.
    #func_fitVector <- func_fitVector/(max(abs(func_fitVector)))
    return( unlist(lapply(1 + (func_fitVector - if(is.null(func_weights)){mean(func_fitVector)}else{weighted.mean(func_fitVector,func_weights)}), function(func_thisFitness){ max(0, func_thisFitness) })) )
  }
}



# If the func_estProp is numeric and between 0 and 1, then we calculate that a line is established
# if it's popSize > sum(popSize) * func_estProp.  Otherwise we try to evaluate the expression passed.  I have no real
# handlers for the validity of this expression. USER BEWARE!
# A common expression uses the suggestions of Desai 2007 where a lineage establishes when it has 1/s individuals
#' This function is used to find which elements of a population matrix are deemed as established.  Established is determined
#' by having a number of individuals greater than or equal to a definable proportion of the summed community size.
#'
#' @param func_inMatrix This is a matrix which must contain at least one column named as func_sizeCol which contains
#' the number of individuals in the communities' populations.  But it may also be required to include a column
#' func_fitCol if func_estProp is "Desai".
#' @param func_sizeCol DO NOT MODIFY - this is the column name that is querried to find population sizes
#' @param func_fitCol DO NOT MODIFY - this is the column name that is querried to find population fitness -
#' only important if func_estProp is set to "Desai"
#' @param func_estProp If this value is less than 1 - This is the proportion of the current community size which
#' is used to define a population as established it returns the rows of.  If this value is greater than 1, it is
#' the minimum number of individuals required before a population is considered as established.
#' Lastly, it can be the character string "Desai", at which point - as per Desai 2007 - a lineage is established
#' once it has 1/s individuals.
#'
#' @return A subset form of the input func_inMatrix matrix object containing the populations which are calculated
#' as established.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
querryEstablished <- function(func_inMatrix, func_sizeCol = "popSize", func_fitCol = "fitness", func_estProp = 0.01){
  if(is.numeric(func_estProp)){
    # If the value is less than or equal to 1 we take it to mean a proportion of the population
    if(func_estProp <= 1){
      return( func_inMatrix[which(func_inMatrix[, func_sizeCol] >= (func_estProp * sum(func_inMatrix[, func_sizeCol]))),] )
      # Otherwise we take it to mean an exact value, we use the floor call here so that a value of ex: 1.1
      # could be used in order to force return that any lineage which arises is established.
    } else {
      return( func_inMatrix[which(func_inMatrix[, func_sizeCol] >= floor(func_estProp)) ,] )
    }
    # NOTE: This will only work as long as all fitness values passed are greater than 1.
  } else if(func_estProp == "Desai") {
    return( func_inMatrix[which(apply(func_inMatrix[, c(func_sizeCol, func_fitCol)],MARGIN=1,function(thisRow){
      thisRow[func_sizeCol] >= 1/log(thisRow[func_fitCol])
    })),] )
  } else {
    return( func_inMatrix[eval(parse(text= func_estProp)),] )
  }
}



#' This is a convenience wrapper for sending an error and ending the SHAPE run as well as the R environment.
#' It will print a message and then traceback() report before pausing and quiting the R session.  This exists
#' to help debugging when SHAPE is run in batch-mode.
#' @param func_message The message to be sent to screen prior to ending the R session.
#'
#' @section Note:
#' There is no example as this functions role is to print a message
#' and then quit the R run.
#'
#' @export
stopError <- function(func_message){
  stop(func_message)
  traceback()
  Sys.sleep(5)
}


#' This is a function to trim a string by removing the first and last character, it's used to trim quotation marks
#' used in the parameter input
#'
#' @param funcIn a vector of character strings which you want trimmed
#'
#' @return character vector of length equal to the input
#'
#' @examples
#' # It removes leading and trailing string positions, use when quotations are known to exist.
#' trimQuotes(c('"someWords"','otherwords"',"is_changed"))
#'
#' @export
trimQuotes <- function(funcIn){ substr(funcIn,2,nchar(funcIn)-1) }

#' This is a function to add quotation marks around each element of a character string vector
#'
#' @param funcIn a vector of character strings which you want padded by quotation marks
#'
#' @return character vector of length equal to the input
#'
#' @export
addQuotes <- function(funcIn){ paste('"',funcIn,'"',sep="") }



#' This is a function to programatically create R batch submission script names
#'
#' @param inVar This is the vector of character string(s) to be used for naming
#'
#' @return A vector of character string of length equal to input.
#'
#' @examples
#' # Returns a standard named string
#' name_subScript(c("myJob","otherContent"))
#'
#' @export
name_subScript <- function(inVar){
  paste("SHAPE_submit_",inVar,".sh",sep="")
}

#' This is a function to programatically create R batch submission script names
#'
#' @param inVar This is the vector of character string(s) to be used for naming
#'
#' @return A vector of character string of length equal to input.
#'
#' @export
name_batchSubmit <- function(inVar){
  paste("submit_jobBatch_",inVar,".sh",sep="")
}

#' This is a function to programatically create R script names
#'
#' @param inVar This is the vector of character string(s) to be used for naming
#'
#' @return A vector of character string of length equal to input.
#'
#' @examples
#' # Returns a standard named string
#' name_bodyScript(c("myJob","otherContent"))
#'
#' @export
name_bodyScript <- function(inVar){
  paste("SHAPE_body_",inVar,".r",sep="")
}

#' This is a function to programatically create R script names
#'
#' @param inVar This is the vector of character string(s) to be used for naming
#'
#' @return A vector of character string of length equal to input.
#'
#' @examples
#' # Returns a standard named string
#' name_parameterScript(c("myJob","otherContent"))
#'
#' @export
name_parameterScript <- function(inVar){
  paste("SHAPE_parameters_",inVar,".r",sep="")
}


#' This quick little function is a means for me to create the strings of
#' environments and subsequently extract information back out.
#'
#' @param func_Index This is the vector of numeric, or otherwise unique ID values for the
#' environments to be created.  Or if funcSplit == TRUE, then these are the names to be split.
#' @param funcSplit A logical toggle of whether you are building or splitting the name
#' @param funcBase This is the character string used as a prefix to identify environment objects
#'
#' @return A vector of character string of length equal to input.
#'
#' @examples
#' # Returns a standard named string
#' test_envNames <- nameEnviron(1:10)
#' nameEnviron(test_envNames, funcSplit = TRUE)
#' @export
nameEnviron <- function(func_Index, funcSplit = FALSE, funcBase = getOption("shape_envString")){
  # This checks if we're asking for the name to be split or not
  if(funcSplit){
    # We return only the second piece of information given the setup of this function's naming practice.
    return( unlist(lapply(strsplit(func_Index,"_e_"),function(x){ x[2] })) )
  } else {
    return( paste(funcBase,func_Index,sep="_e_") )
  }
}



#' This quick little function is a means for me to create the strings of
#' environments and subsequently extract information back out.
#'
#' @param func_inString This is the vector of numeric, or otherwise unique ID values for the
#' environments to be created.  Or if funcSplit == TRUE, then these are the names to be split.
#' @param func_splitStr A logical toggle of whether you are building or splitting the name
#' @param func_inPrefix This is the character string used as a prefix to identify environment objects
#'
#' @return A vector of character string of length equal to input.
#'
#' @examples
#' # Returns a standard named string
#' test_objectNames <- nameObject(1:10, "testObject")
#' nameObject(test_objectNames, "testObject", func_splitStr = TRUE)
#' @export
nameObject <- function(func_inString, func_inPrefix, func_splitStr = FALSE){
  if(func_splitStr){
    unlist(lapply(strsplit(func_inString, func_inPrefix),function(x){ x[length(x)] }))
  } else {
    paste(func_inPrefix, func_inString,sep="")
  }
}





#' This function is designed to establish an initial object which maps the fitness values
#' of genome positions based on the state of that site.  At present, this has no meaning
#' if the model of simulation is no NK, Additive, or Fixed.  Where the first is Kauffman's NK
#' model and form of calculations, Additive is what that word would make you think for fitness
#' effects of mutations at sites, and Fixed is when user supplied a defined fitness matrix that
#' describes the entire fitness landscape.  NOTE: This function should likely be called without
#' supplying any non-default arguments as it will use the shape_ options defined.
#'
#' @param func_simModel This is the fitness landscape model being simulated
#' @param func_const_fixedFrame This is a contextual object that described constant fitness effects
#' @param func_const_siteStates These are the posibble states for genome sites, at present this
#' ought to be "0" and/or "1"
#'
#' @return A contextually meaningful matrix describing fitness effects of mutations/genotypes,
#' where based on the context NULL may be returned.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
set_siteByState_fitnessMat <- function(func_simModel = getOption("shape_simModel"),
                                       func_const_fixedFrame = getOption("shape_const_fixedFrame"),
                                       func_const_siteStates = getOption("shape_const_siteStates")){
  # This build the item that will be returned
  func_tmpReturn <- if(is.element(func_simModel,c("NK","Additive"))){
                      # This is for the NK or Additive forms of simulation modelling
                      sapply(func_const_siteStates,function(thisState){
                        # If we're using the NK model then the ancestral state is the ancestral fitness value, for additive
                        # this is the selection coefficient.
                        return(  if(thisState == "0" && getOption("shape_const_relativeFitness")){
                          rep(if(func_simModel == "NK"){
                                getOption("shape_const_ancestFitness")
                              } else if(func_simModel == "Additive"){
                                0
                              },
                              getOption("shape_genomeLength"))
                          # If we're considering the derived state, then we use the distribution information
                        } else {
                          tmpReturn <- fitnessDist(getOption("shape_genomeLength"),
                                                   tmpDistribution = getOption("shape_constDist"),
                                                   tmpParameters = getOption("shape_const_distParameters"))
                          # If the distribution draws relative fitnesses but this needs to be used as selection coefficients...
                          if(getOption("shape_const_distAsS")){
                            tmpReturn <- tmpReturn - 1
                          }
                          tmpReturn
                        }
                        ) })
                    } else if(func_simModel == "Fixed") {
                      # This assumes the func_const_fixedFrame object is a table created by the user and
                      # with the colnames binaryString and fitness
                      func_tmpFileName <- paste(getOption("shape_workDir"),func_const_fixedFrame,sep="")
                      if(file.exists(func_tmpFileName)){
                        # In this case we exchange what should be a filename to be the value loaded
                        func_tmpLandscape <- read.csv(func_tmpFileName,
                                                      header=TRUE,
                                                      stringsAsFactors = FALSE)
                        options("shape_const_fixedFrame"= func_tmpLandscape)
                      } else {
                        stopError(paste("Could not find the fixed fitness landscape file: ",func_tmpFileName," please review",sep=""))
                      }
                      if(all(is.element(colnames(func_tmpLandscape),c("binaryString","fitness")))){
                        # Ok this means we need to convert what was passed into a named vector.  It is assumed the information passed
                        # for binaryString is a sequence of 0 and 1:  ex:   000,   010,   100,  etc.....
                        tmpReturn <- func_tmpLandscape[,"fitness"]
                        names(tmpReturn) <- unlist(lapply(strsplit(func_tmpLandscape[,"binaryString"],""),function(x){
                                                        if(any(x == "1")){
                                                          paste(which(x == "1"),sep="",collapse=getOption("shape_sepString"))
                                                        } else {
                                                          ""
                                                        }
                                                      }))
                        tmpReturn
                      } else {
                        stopError("Fixed fitness landscape matrix information does not contain the two column names binaryString and fitness")
                      }
                    # This means there is no meaningful context to this value
                    } else {
                      NULL
                    }
  # If a siteBystate object was created then we need to give it some column names, also we ensure it is a matrix
  # It may not be a matrix if our genome length = 1.... which is used for trivial purposes but can happen.
  if(is.element(func_simModel,c("NK","Additive"))){
    if(!is.matrix(func_tmpReturn)){
      if(length(func_tmpReturn) != 2){
        stopError(paste("There was a problem trying to create the <shape_const_siteBystate_fitnessMat> matrix, it's length was ",
                        length(func_tmpReturn),
                       "please review",
                       sep=""))
      } else {
        func_tmpReturn <- matrix(func_tmpReturn,
                                nrow=1,
                                ncol= length(func_const_siteStates))
      }
    }
    # Now the func_tmpReturn should be a matrix with a number of  columns equal to the number of states
    colnames(func_tmpReturn) <- func_const_siteStates
  }
  # We return the value
  return( func_tmpReturn )
}


#' This is a convenience function for setting the dependent fitness values of sites in
#' an NK fitness landscape model.  This allows the dependent fitness of sites to be calculated
#' once and then referenced as mutations occur.  It makes exploring this style of fitness landscape
#' a bit more computationally friendly - as it generally isn't.
#'
#' @param func_simModel This is the fitness landscape model being simulated
#' @param func_const_siteBystate_fitnessMat This is the sitewise independent fitness contributions in the fitness landscape
#' @param func_const_NK_interactionMat This defines the sitewise dependencies based on the K interactions.
#'
#' @return Either the dependent sitewise fitness contributions in an NK fitness landscape, or NA.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
set_DepbySite_ancestFitness <- function(func_simModel = getOption("shape_simModel"),
                                        func_const_siteBystate_fitnessMat = getOption("shape_const_siteBystate_fitnessMat"),
                                        func_const_NK_interactionMat = getOption("shape_const_NK_interactionMat")){
  # We collect the value to be returned, this only has meaning if we're simulating NK fitness landscape
  func_tmpReturn <- if(func_simModel == "NK"){
                      sapply(1:nrow(func_const_siteBystate_fitnessMat),function(thisSite){
                          mean(func_const_siteBystate_fitnessMat[c(thisSite, func_const_NK_interactionMat[thisSite,]),"0"])
                        })
                    } else {
                      NA
                    }
  # We return the value
  return( func_tmpReturn )
}


#' In a RMF fitness landscape model, there is a weighting value applied to the independent
#' fitness contribution term.  This function calculates that value for the run
#'
#' @param func_simModel This is the model of fitness landscape being considered
#' @param func_numDraws This is the number of draws taken from the independent term's distribution
#' so that we can identify the amount of variance in that distribution.  It should be a large integer -- eg 5e7
#' @param func_distType This is the distribution string reference for this run
#' @param func_distParms These are the parameters for this runs distribution function
#' @param func_const_RMF_theta This is the theta value which is multiplied to the variance in the distribution.
#' The value returned will be a product of this numeric and the variance calulated.
#' From Neidhart 2014 theta is measured as:  theta = c / sqrt var random_component
#' and so if we want to calculate "c" we return the product of theta and sqrt of variance in the distribution
#'
#' @return A single numeric value, which may be NA if a non Rough Mount Fuji model is being simulated
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
set_RMF_indWeight <- function(func_simModel = getOption("shape_simModel"), func_numDraws = 1e8,
                              func_distType = getOption("shape_constDist"), func_distParms = getOption("shape_const_distParameters"),
                              func_const_RMF_theta = getOption("shape_const_RMF_theta")){
  # We define our value and then return it.
  func_tmpReturn <- if(func_simModel == "RMF"){
                      # We make a large call to our fitness function, find the variance, and then use this to help set our value
                      # Why?  Because of literature on RMF fitness landscapes suggesting a range as observed in nature
                      # being relative to some theta value and variance in the distribution of terms.
                      func_const_RMF_theta * sqrt(var(fitnessDist(func_numDraws,
                                                                   tmpDistribution = func_distType,
                                                                   tmpParameters = func_distParms)))
                    } else {
                      NA
                    }

  return( func_tmpReturn )
}



#' This is a function to just return a matrix that defines the sitewise dependencies
#' for an NK fitness landscape.  If K == 0 or, this is not an NK simulation, it return NULL
#'
#' @param func_simModel This is the fitness landscape model being simulated
#' @param func_genomeLength This is the number of sites in the genome being simulated
#' @param func_numInteractions An integer value defining the number of sites
#' that interact with each other site
#'
#' @return Either NULL, or a matrix with K + 1 columns, detailing the sites
#' interacting with a focal site - identified by the row number and the
#' cell values of the columns.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
set_const_NK_interactionsMat <- function(func_simModel = getOption("shape_simModel"), func_genomeLength = getOption("shape_genomeLength"),
                                         func_numInteractions = getOption("shape_const_numInteractions")){
  # Calculate the reutrn and then return it.
  func_tmpReturn <- if(func_simModel != "NK" || func_numInteractions == 0){
                      NULL
                    } else if(func_simModel == "NK" && func_numInteractions == 1){
                      as.matrix(sapply(1:func_genomeLength,function(thisSite){
                                    sample((1:func_genomeLength)[-thisSite],
                                               func_numInteractions,
                                               replace=FALSE)
                                  }),ncol=1)
                    } else {
                      t(sapply(1:func_genomeLength,function(thisSite){
                            sample((1:func_genomeLength)[-thisSite],
                                   func_numInteractions,
                                   replace=FALSE)
                          }))
                    }

  return( func_tmpReturn )
}



#' This function samples the space of all possible genotypes and then defines
#' one that will be considered as the independent fitness contribution global optima.
#'
#' @param func_simModel This is the fitness landscape model being simulated
#' @param func_genomeLength The number of sites in the genome being simulated
#' @param func_initDistance This is the number of mutations found in the global optimal genotype
#' @param func_sepString This is the string collapse separator used in the run
#'
#' @return A character string of genome positions at which there ought to be mutations to be optimal
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
set_const_RMF_globalOptima <- function(func_simModel = getOption("shape_simModel"), func_genomeLength = getOption("shape_genomeLength"),
                                       func_initDistance = getOption("shape_const_RMF_initiDistance"),
                                       func_sepString = getOption("shape_sepString")){
  # Step 1 - we sample the genome length to know what positions will be mutated at optima
  func_tmpReturn <- if(func_simModel == "RMF"){
                      # Now, if the wild-type - or zero mutations - is the optima this is a special case
                      # and since in SHAPE the zero mutant type is recorded as "", we return that here.
                      tmpReturn <- if(func_initDistance == 0){
                                      ""
                                    } else {
                                      sample(1: func_genomeLength, func_initDistance,replace=FALSE)
                                    }
                      paste(tmpReturn[order(tmpReturn,decreasing=FALSE)],collapse= func_sepString)
                    } else {
                      NA
                    }
  return( func_tmpReturn )
}





#' This is the function that runs the main body, or meaningful execution, of SHAPE experiments.
#' In other words this is the main work-horse function that calls all the other parts and will
#' execute you simulation run.  It has the main parts of:
#' 1. Stochastic Events;
#' 2. Deaths;
#' 3. Births;
#' 4. Mutations;
#' and during mutations this is where the mutational landscape is queried and updated as required.
#' NOTE: Many of its internal operations are controlled by options with the suffix "shape_" and are
#' not explicitly passed as arguments at call to this function.
#' @param func_inputFrames This is a list of data.frames, either 1 or 2 elements, reporting on the last
#' one or two steps in the simulation.
#' @param func_currStep This is an integer value counting the absolute step in the simulation, its value is never reset.
#' @param func_stepCounter This is an integer value which is a counter in the most tradititional sense.
#' It's job is to track if it's time for a Stochastic event to trigger and its value is reset at that point.
#' @param func_growthModel This is the growth model of the SHAPE run, it is passed here as a computational
#' convenience since it is used numerous times in the function
#' @param func_growthRate This is the growth rate of the SHAPE run, it is passed here as a computational
#' convenience since it is used numerous times in the function
#' @param func_landscapeModel This is the fitness landscape model of the SHAPE run, it is passed here as a computational
#' convenience since it is used numerous times in the function
#' @param func_fileName_dataBase This is the filepaths of DBs of the SHAPE run, it is passed here as a computational
#' convenience since it is used numerous times in the function
#'
#' @return Returns a new list of 2 data.frames reporting on the state of SHAPE community for the last 2 time steps -
#' ie: the one just run, and the most prior step.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
runReplicate <- function(func_inputFrames, func_currStep, func_stepCounter,
                         func_growthModel = getOption("shape_const_growthForm"),
                         func_growthRate = getOption("shape_const_growthRate"),
                         func_landscapeModel = getOption("shape_simModel"),
                         func_fileName_dataBase = getOption("shape_fileName_dataBase")){
  # We define the working frame as the data.frame of the last step
  workFrame <- func_inputFrames[[nameTable_step(func_currStep-1)]]
  # This loop is where the disturbance tracking occurs, this could be made into a function given it writes to global.
  if(func_stepCounter > getOption("shape_track_distSize")[nrow(getOption("shape_track_distSize")),"stepReq"]){
    # We check which populations have at least a half individual otherwise they aren't considered alive.
    # This value is important as it is the least I will consider meaningful wether or not populations are
    # being tracked as either whole or fractional entities.
    currLineages <- unname(which(workFrame[,"popSize"] >= 0.5))
    if(length(currLineages) == 0){
      print(paste("All lineages have died as of timeStep ", func_currStep,sep=""))
      return( func_inputFrames )
    }
    # This is an unfortunate requirement to make certain that we've got a workFrame which is in the proper format
    workFrame <- reportPopulations(func_numMuts= workFrame[currLineages,"numMuts"],
                                   func_genotypeID= workFrame[currLineages,"genotypeID"],
                                   func_popSizes= workFrame[currLineages,"popSize"],
                                   func_fitnesses= workFrame[currLineages,"fitness"],
                                   func_births=workFrame[currLineages,"births"],
                                   func_deaths=workFrame[currLineages,"deaths"],
                                   func_mutants=workFrame[currLineages,"mutants"],
                                   func_progenitor=workFrame[currLineages,"progenitor"])
    # Now that we've trimmed our workFrame to the living lineage we can perform lossSampling.
    tmp_lossFactors <- compute_distGrowth(func_distFactor = if(func_growthModel == "logistic"){
                                                              getOption("shape_init_distPars")
                                                            }else{
                                                              0
                                                            },
                                          func_growthType = func_growthModel,
                                          func_distType = getOption("shape_const_distType"),
                                          func_growthRate = func_growthRate,
                                          func_popSize = workFrame[,"popSize"],
                                          func_focalSize = getOption("shape_const_focal_popValue"),
                                          func_manualGenerations = if(func_growthModel == "constant"){
                                                                      getOption("shape_numGenerations")
                                                                    } else {
                                                                      getOption("shape_const_growthGenerations")
                                                                    },
                                          func_stepDivs = getOption("shape_size_timeStep"))

    # If the dilution factor is not greater than 1 we won't be performing any dilution.
    if(tmp_lossFactors["factor"] > 1){
      # Now that we've calculated the dilution factor that will be used we call our lossSampling function to see which individuals will remain
      tmp_sampledPop <- lossSampling(func_inPopulation = workFrame[,"popSize"],
                                     func_dilutionFactor = 1/if(tmp_lossFactors["factor"] == 0){1}else{tmp_lossFactors["factor"]})
      # We update our loss tracker firstby the amount of loss realised, then by the actual factor of loss this represented
      tmp_lossFactors["popLost"] <- sum(workFrame[,"popSize"] - tmp_sampledPop)
      tmp_lossFactors["factor"] <- sum(workFrame[,"popSize"])/sum(tmp_sampledPop)
      workFrame[,"popSize"] <- tmp_sampledPop
      # We track the population change which occurs as a result of this change
      options("shape_track_distSize" = rbind(getOption("shape_track_distSize"), c("Step"= func_currStep, tmp_lossFactors)))
    }
    # We update our counter_logSteps object back to the initial value.
    func_stepCounter <- 1
    options("shape_counter_logSteps" = func_stepCounter)
  }
  # We check which populations have at least a half individual otherwise they aren't considered alive.
  currLineages <- as.vector(which(workFrame[,"popSize"] >= 0.5, useNames=FALSE))
  if(length(currLineages) == 0){
    print(paste("All lineages have died as of timeStep ", func_currStep,sep=""))
    return( func_inputFrames )
  }
  # This is an unfortunate requirement to make certain that we've got a workFrame which is in the proper format
  workFrame <- reportPopulations(func_numMuts= workFrame[currLineages,"numMuts"],
                                 func_genotypeID= workFrame[currLineages,"genotypeID"],
                                 func_popSizes= workFrame[currLineages,"popSize"],
                                 func_fitnesses= workFrame[currLineages,"fitness"],
                                 func_births=workFrame[currLineages,"births"],
                                 func_deaths=workFrame[currLineages,"deaths"],
                                 func_mutants=workFrame[currLineages,"mutants"],
                                 func_progenitor=workFrame[currLineages,"progenitor"])
  # We calculate the growth for each of our lineage(s)
  num_birthDeaths <- growthFunction(func_inSize = workFrame[,"popSize"],
                                    func_inFitness = workFrame[,"fitness"],
                                    func_bProb = getOption("shape_const_birthProb"),
                                    func_dProb = getOption("shape_const_deathProb"),
                                    func_deathDen_logical = getOption("shape_death_byDensity"),
                                    func_deathDen_max = getOption("shape_death_densityCap"),
                                    func_deathDen_power = getOption("shape_death_densityCorrelation"),
                                    func_sizeStep = getOption("shape_size_timeStep"),
                                    func_growthForm = func_growthModel,
                                    func_carryingCapacity = getOption("shape_const_focal_popValue"),
                                    func_basalRate = func_growthRate,
                                    func_deathScale = getOption("shape_scaleGrowth_byDeaths"),
                                    func_drift = getOption("shape_includeDrift"),
                                    func_roundValues = getOption("shape_track_asWhole"),
                                    func_inIDs = rownames(workFrame))
  #print("starting reportPopulations") # This is for reporting on run times but is not needed....
  # Now, because the number of births may be negative (for reasons of constant growth or growth function with carrying capacities adjusting sizes - its rare...)
  # we adjust the num_birthDeaths so that negative births are shuffled to deaths.
  tmpUpdates <- which(num_birthDeaths[,"births"] < 0)
  if(length(tmpUpdates) > 0){
    # We subtract the negative births from (basically add to...) the amount of deaths
    num_birthDeaths[tmpUpdates,"deaths"] <- num_birthDeaths[tmpUpdates,"deaths"] - num_birthDeaths[tmpUpdates,"births"]
    num_birthDeaths[tmpUpdates,"births"] <- 0
  }

  # Now is the section where we verify the number of mutants which may arise in each lineage
  # We check if there are any lineages which could be selected
  tmp_mutableLineages <- if(!getOption("shape_allow_backMutations") &&
                            is.element(getOption("shape_genomeLength"),workFrame[,"numMuts"])){
    # If we don't allow back mutations and there is a mutant which is fully derived we don't consider it
    rownames(num_birthDeaths)[which(workFrame[rownames(num_birthDeaths),"numMuts"] !=
                                      getOption("shape_genomeLength"))]
  } else {
    # Otherwise all lineages can be considered
    rownames(num_birthDeaths)
  }
  # A lineage cannot be mutable if, after births and deaths, the population size is not at least 1
  # This is because from a practical standpoint, we've balanced growth possibly by scaling deaths AND
  # since the generation of a mutant results in subtraction from the parental lineage, we can't allow
  # the maths to not add up...
  tmp_mutableLineages <- tmp_mutableLineages[which(workFrame[tmp_mutableLineages,"popSize"] +
                                                     num_birthDeaths[tmp_mutableLineages,"births"] -
                                                     num_birthDeaths[tmp_mutableLineages,"deaths"] >= 1)]

  # So for all lineages which could generate a mutant, we call the mutation process to define the number which arise
  # Since a birth events include not just the offpsring but also the parental type, we multiply the number which
  # may have gained a mutation during replication by the basic growth rate value.  If mutations occur not only at birth
  # then we adjust the pop-Size by the (births * (basic growth rate - 1))
  proposedMutants  <- mutationFunction(func_inSize = unlist(lapply(num_birthDeaths[tmp_mutableLineages,"births"] *
                                                                     func_growthRate/
                                                                     (func_growthRate - 1),function(x){ max(0,x) })),
                                       func_inProb = getOption("shape_const_mutProb")) +
    if(getOption("shape_muts_onlyBirths")){
      0
    } else {
      mutationFunction(func_inSize = unlist(lapply(workFrame[tmp_mutableLineages,"popSize"] -
                                                     num_birthDeaths[tmp_mutableLineages,"deaths"] -
                                                     num_birthDeaths[tmp_mutableLineages,"births"] *
                                                     1/(func_growthRate - 1), function(x){ max(0,x) })),
                       func_inProb = getOption("shape_const_mutProb") * getOption("shape_size_timeStep"))
    }
  names(proposedMutants) <- tmp_mutableLineages
  # We now udpate our num_birthDeaths object with the mutants generated but noting that a lineage
  # cannot generate more mutants than it has living memebers.
  num_birthDeaths <- cbind(num_birthDeaths,
                           "mutants"= as.vector(sapply(rownames(workFrame),function(x){
                             if(is.element(x,names(proposedMutants))){
                               return( min(proposedMutants[x],
                                           workFrame[x,"popSize"] +
                                             num_birthDeaths[x,"births"] -
                                             num_birthDeaths[x,"deaths"]) )
                             } else {
                               0
                             }
                           })))
  # We initialise the object for reporting on which mutants are created.  This chunk of code
  # is one of the most computationally expensive parts of the simulations, after the cost
  # of querrying the DB files once they get to be large.
  all_newMutants <- NULL
  if(sum(num_birthDeaths[,"mutants"]) > 0){
    for(thisLineage in rownames(num_birthDeaths)[which(num_birthDeaths[,"mutants"] > 0)] ){
      # This re-opens the connections so we can access the fitness landscape object.
      connections_dataBase <- sapply(names(func_fileName_dataBase),function(x){
        reset_shapeDB(paste(getOption("shape_outDir"), func_fileName_dataBase[[x]],sep=""),
                      func_type = "connect")
      })
      # We define the focal genotype we're working with
      this_focalGenotype <- retrieve_binaryString(func_genotypeID = as.numeric(thisLineage),
                                                  func_numMuts = workFrame[thisLineage,"numMuts"],
                                                  func_subNaming = getOption("shape_db_splitTables"),
                                                  func_landscapeCon = connections_dataBase$genotypeSpace)
      # This is an unfortunate work-around for an instance where I've found that a lineage had been miss-classified by numMuts
      # which caused there to be no value returned, I similarly check for multiple instances....
      # I believe to have fixed all situations that would have led to this
      if(length(this_focalGenotype) == 0 || nrow(this_focalGenotype) != 1){
        # I'm interested in how often this sanity check is used...
        print(paste("Safety net for mis-classified genotypeID has been used for ", this_focalGenotype," on step ",func_currStep,sep=""))
        # We look through all of the tables and return the one who has a search result for our genotypeID
        tmp_mutSearch <- sapply(nameTable(RSQLite::dbListTables(connections_dataBase$genotypeSpace), func_splitName = TRUE), function(this_numMuts){
          return( nrow(retrieve_binaryString(func_genotypeID = as.numeric(thisLineage),
                                             func_numMuts = as.numeric(this_numMuts),
                                             func_subNaming = getOption("shape_db_splitTables"),
                                             func_landscapeCon = connections_dataBase$genotypeSpace)) )
        })
        # Now if we've found only one table which held our value, we update the workFrame's numMuts and re-run our call.
        if(length(which(tmp_mutSearch == 1)) == 1){
          workFrame[thisLineage,"numMuts"] <- as.numeric(names(tmp_mutSearch)[which(tmp_mutSearch == 1)])
          # We now update the value for this_focalGenotype, as we should have been able to define a proper table.
          this_focalGenotype <- retrieve_binaryString(func_genotypeID = as.numeric(thisLineage),
                                                      func_numMuts = workFrame[thisLineage,"numMuts"],
                                                      func_subNaming = getOption("shape_db_splitTables"),
                                                      func_landscapeCon = connections_dataBase$genotypeSpace)
        } else {
          stop(paste("There was a problem while trying to define a unique table containing ",thisLineage," during ",func_currStep," please review",sep=""))
        }
      }
      # We check to see if thisLineage is already within the neighbourhood reference database, this is a time saving database.
      tmp_allNeighbours <- NULL
      if(is.element(nameTable_neighbourhood(thisLineage),RSQLite::dbListTables(connections_dataBase$nearestNeighbours))){
        tmp_allNeighbours <- dbGetQuery(connections_dataBase$nearestNeighbours,
                                        paste("SELECT * FROM ",nameTable_neighbourhood(thisLineage),sep=""))$neighbours
      } else {
        # We define all the possible nearestNeighbours for this lineage since they've not been stored.
        tmp_allNeighbours <- defineNeighbours(func_tmpGenotype = this_focalGenotype[1,"binaryString"],
                                              func_tmpDirection = getOption("shape_allow_backMutations"))
        # Now if the size of thisLineage is large enough we store this in the Neighbourhood shortcut database.
        if(workFrame[thisLineage,"popSize"] >= getOption("shape_const_hoodThresh")){
          dbWriteTable(connections_dataBase$nearestNeighbours,
                       nameTable_neighbourhood(thisLineage),
                       data.frame("neighbours"=tmp_allNeighbours))
        }
      }

      # We see if this mutant's neighbourhood has been explored, if not we need to add it to the fitness landscape
      if(!as.logical(this_focalGenotype[1,"isExplored"])){
        # This uses the fitness landscape models to calculate the fitness value for all mutants in this neighbourhood.
        createGenotypes(tmp_focalGenotype = this_focalGenotype[1,"binaryString"],
                        tmp_focalFitness = if(is.element(func_landscapeModel,c("NK","Additive"))){
                                              getOption("shape_const_siteBystate_fitnessMat")
                                            } else {
                                              workFrame[thisLineage,"fitness"]
                                            },
                        maxHamming = getOption("shape_max_numMutations"),
                        tmp_landModel = func_landscapeModel,
                        tmp_relativeFitness = getOption("shape_const_relativeFitness"),
                        tmpDistribution = getOption("shape_constDist"),
                        tmpParameters = getOption("shape_const_distParameters"),
                        tmp_currNeighbours = tmp_allNeighbours,
                        tmp_genCon = connections_dataBase$genotypeSpace,
                        tmp_tableSplit = getOption("shape_db_splitTables"),
                        tmp_maxRows = getOption("shape_maxRows"),
                        tmp_distAsS = getOption("shape_const_distAsS"))
      }
      # Now we generate our mutants by sampling all possible neighbours with replacement
      newMutants <- table(sample(tmp_allNeighbours,
                                 num_birthDeaths[thisLineage,"mutants"], replace = TRUE))
      # Now for the mutants that were drawn from the neighbourhood, we go and retrieve the information from the fitness landscape database.
      # To know where to find our mutant information we find the number of mutations for the mutant, then go search those for the genotype.
      tmp_numMuts <- unlist(lapply(strsplit(names(newMutants),getOption("shape_sepString")),length))
      tmp_dbTables <- RSQLite::dbListTables(connections_dataBase$genotypeSpace)
      tmp_dbTables <- tmp_dbTables[unique(unlist(lapply(nameTable(unique(tmp_numMuts)),function(thisString){
                                                    which(grepl(thisString,tmp_dbTables))
                                                  })))]
      # Here we get the information for the mutant by querrying our SQL database
      # This is the actual string for the call to information
      tmp_newStrings <- gsub("[[:space:]]","",paste("\'", names(newMutants),"\'",collapse=','))
      # This is the call to get our information about mutants.  It's followed by a check that all the relevant information was found
      # If it wasn't we try to create the missing pieces.
      tmp_infoAdd <- as.matrix(dbGetQuery(connections_dataBase$genotypeSpace,
                                          paste("SELECT genotypeID,fitness,binaryString FROM ",
                                                tmp_dbTables,
                                                ' WHERE binaryString IN (',
                                                tmp_newStrings,
                                                ')',
                                                collapse=" UNION ")))
      # If we didn't find some of the mutant information we try again after trying to re-fill the fitness landscape
      # This ought not to be an issue or be required.
      if(!all(is.element(names(newMutants),tmp_infoAdd[,"binaryString"]))){
        # We explore the neighbouring space of this focal lineage
        createGenotypes(tmp_focalGenotype = this_focalGenotype[1,"binaryString"],
                        tmp_focalFitness = if(is.element(func_landscapeModel,c("NK","Additive"))){
                                              getOption("shape_const_siteBystate_fitnessMat")
                                            } else {
                                              workFrame[thisLineage,"fitness"]
                                            },
                        maxHamming = getOption("shape_max_numMutations"),
                        tmp_landModel = func_landscapeModel,
                        tmp_relativeFitness = getOption("shape_const_relativeFitness"),
                        tmpDistribution = getOption("shape_constDist"),
                        tmpParameters = getOption("shape_const_distParameters"),
                        tmp_currNeighbours = tmp_allNeighbours,
                        tmp_genCon = connections_dataBase$genotypeSpace,
                        tmp_tableSplit = getOption("shape_db_splitTables"),
                        tmp_distAsS = getOption("shape_const_distAsS"))
        # We go to collect the information again... it should work now.
        tmp_infoAdd <- as.matrix(dbGetQuery(connections_dataBase$genotypeSpace,
                                            paste("SELECT genotypeID,fitness,binaryString FROM ",
                                                  tmp_dbTables,
                                                  ' WHERE binaryString IN (',
                                                  tmp_newStrings,
                                                  ')',
                                                  collapse=" UNION ")))
      }
      for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
      # Now we should, no matter what, have information for our mutant and add it to the mutant tracking matrix
      tmp_newOrder <- unlist(lapply(names(newMutants),function(x){ which(tmp_infoAdd[,"binaryString"] == x) }))
      newMutants <- cbind(tmp_numMuts, tmp_infoAdd[tmp_newOrder,"genotypeID"], newMutants, tmp_infoAdd[tmp_newOrder,"fitness"])
      colnames(newMutants) <- getOption("shape_popMat_colnames")
      rownames(newMutants) <- as.character(newMutants[,"genotypeID"])
      # This updates the mutant tracking matrix
      all_newMutants <- rbind(all_newMutants,
                              cbind(newMutants,
                                    "progenitor"= unlist(lapply(newMutants[,"popSize"],function(thisSize){
                                      paste(thisLineage,thisSize,sep=getOption("shape_sepString"))
                                    })) ) )
    } # This closes out the thisLineage split for each row of num_birthDeaths which has a mutant
  } # This closes out the conditional that there are mutants to have entered this loop.

  # We update our tmp_stepChanges with the number of births, deaths and mutations which occured for our existing (workFrame) lineages
  # We also report on new mutants and which lineages birthed them new mutants, this is achieved through our progenitorReport object update in the newMutants section
  # We use the progenitorReport object to easily reference which lineages

  # First we report on all the existing lineages, this is the easy part.
  tmp_stepChanges <- data.frame(cbind(workFrame[rownames(num_birthDeaths),
                                                getOption("shape_popMat_colnames")],
                                      num_birthDeaths,
                                      "progenitor"="",stringsAsFactors=FALSE)
                                ,stringsAsFactors=FALSE)
  # If there are no newMutants we can ignore all of this since the changes are just birth based...
  if(!is.null(all_newMutants)){
    # For reasons I can't track, if two separate progenitor lineages produce the same mutant genotypeID
    # I've found that the ID's are not always tracked properly.  This fixes that.
    all_newMutants[,"genotypeID"] <- gsub("[[:space:]]","",all_newMutants[,"genotypeID"])
    # Then we check if any of these newMutants come from different progenitors which would mean the same ID appears twice.
    tmp_repIDs <- table(all_newMutants[,"genotypeID"])
    tmp_repIDs <- names(which(tmp_repIDs > 1))
    # If there are any replicated genotypeID's then we'll want to merge those rows
    if(length(tmp_repIDs) > 0){
      tmp_removeRows <- NULL
      for(thisID in tmp_repIDs){
        tmpRows <- as.vector(which(all_newMutants[,"genotypeID"] == thisID))
        # We update the first row to be the combination of all rows
        all_newMutants[tmpRows[1],"popSize"] <- as.character(sum(as.numeric(all_newMutants[tmpRows,"popSize"])))
        all_newMutants[tmpRows[1],"progenitor"] <- paste(all_newMutants[tmpRows,"progenitor"],collapse= getOption("shape_collapseString"))
        # then flag the other rows to be removed
        tmp_removeRows <- c(tmp_removeRows, tmpRows[-1])
      }
      # Checking there are rows to remove is a sanity check but should always be true
      if(!is.null(tmp_removeRows)){
        all_newMutants <- matrix(all_newMutants[-tmp_removeRows,],ncol=ncol(all_newMutants),
                                 dimnames=list(rownames(all_newMutants[-tmp_removeRows]),colnames(all_newMutants)))
      }
    }
    # Now since we're about to add our all_newMutants information to the workFrame and tmp_stepChanges data.frames, we put it in the same format
    all_newMutants <- data.frame(all_newMutants,stringsAsFactors=FALSE)
    for(thisCol in getOption("shape_popMat_colnames")){
      all_newMutants[, thisCol] <- as.numeric(all_newMutants[, thisCol])
    }

    # Now we update the stepChanges information by looking at which mutants are new to the population or existing
    # For those which exist we simply update the popSize, for new ones we add rows
    tmpTypes <- is.element(all_newMutants[,"genotypeID"], tmp_stepChanges[,"genotypeID"])
    # For the TRUE responses we're updating our existing matrix so we find the similarly ordered rows
    tmpExisting <- if(any(tmpTypes)){
      unlist(lapply(all_newMutants[which(tmpTypes),"genotypeID"], function(x){
        which(tmp_stepChanges[,"genotypeID"] == x)
      }))
    } else {
      NULL
    }
    # If anything already exists then we go through each and simply update the tmp_stepChanges row(s)
    if(!is.null(tmpExisting)){
      # Since the tmpExisting is based on the order of the which(tmpTypes), we can call the additions, or combinations,
      # as vectors of the matrices

      # We update the births to be the sum of the existing births and the newly arived mutants
      tmp_stepChanges[tmpExisting,"births"] <- apply(cbind(tmp_stepChanges[tmpExisting,"births"],
                                                           all_newMutants[which(tmpTypes),"popSize"]),MARGIN=1,sum)
      # We update the progenitor information to be the existing information plus the newly arrived mutants.
      tmp_stepChanges[tmpExisting,"progenitor"] <- apply(cbind(tmp_stepChanges[tmpExisting,"progenitor"],
                                                               all_newMutants[which(tmpTypes),"progenitor"]), MARGIN = 1, function(x){
                                                                 paste(x,collapse=getOption("shape_collapseString"))
                                                               })
      # We now remove the which(tmpTypes) and ensure it is still a data.frame of proper typed data.
      all_newMutants <- all_newMutants[-which(tmpTypes),]
    }

    # For mutant genotypeID's that aren't in the population we are adding information to the population matrix.
    tmp_stepChanges <- rbind(tmp_stepChanges,
                             reportPopulations(func_numMuts= all_newMutants[,"numMuts"],
                                               func_genotypeID= all_newMutants[,"genotypeID"],
                                               func_popSizes= all_newMutants[,"popSize"],
                                               func_fitnesses= all_newMutants[,"fitness"],
                                               func_births= all_newMutants[,"popSize"],
                                               func_deaths= rep(0,nrow(all_newMutants)),
                                               func_mutants= rep(0,nrow(all_newMutants)),
                                               func_progenitor= all_newMutants[,"progenitor"]))
  }
  # Now we update the popSizes in func_currStep where the popSize for a lineage changes by births - deaths, mutants are already accounted for above.
  # We store our changes in a reporting matrix so that we adjust our large report_stepStates only once.
  # NOTE: We only need to update the births and deaths for those lineage which existed before mutations occured, as new mutants have popSize given by their generation
  tmp_stepChanges[rownames(num_birthDeaths),"popSize"] <- tmp_stepChanges[rownames(num_birthDeaths),"popSize"] +
    tmp_stepChanges[rownames(num_birthDeaths),"births"] -
    tmp_stepChanges[rownames(num_birthDeaths),"deaths"] -
    tmp_stepChanges[rownames(num_birthDeaths),"mutants"]

  connections_dataBase <- sapply(names(func_fileName_dataBase),function(x){
    reset_shapeDB(paste(getOption("shape_outDir"), func_fileName_dataBase[[x]],sep=""),
                  func_type = "connect")
  })
  # We ensure the population state, at the end of this step, is in proper format shape
  tmp_returnMatrix <- setNames(list(reportPopulations(func_numMuts= tmp_stepChanges[,"numMuts"],
                                                      func_genotypeID= tmp_stepChanges[,"genotypeID"],
                                                      func_popSizes= tmp_stepChanges[,"popSize"],
                                                      func_fitnesses= tmp_stepChanges[,"fitness"],
                                                      func_births= tmp_stepChanges[,"births"],
                                                      func_deaths= tmp_stepChanges[,"deaths"],
                                                      func_mutants= tmp_stepChanges[,"mutants"],
                                                      func_progenitor= tmp_stepChanges[,"progenitor"])),
                                    nameTable_step(func_currStep))
  # We write out the state of the population at the end of this step.
  dbWriteTable(connections_dataBase$timeStep_States,
               name = names(tmp_returnMatrix),
               value = tmp_returnMatrix[[1]])
  # We close the connection.
  for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
  # So, we return the current state and up to the last state
  return( c(func_inputFrames[[length(func_inputFrames)]],
            tmp_returnMatrix)  )
}




#' This is the actual running of shape, it will initialise objects and values which are calculated from
#' the parameters that have been set - see the options with the suffix 'shape_'.  It will establish the
#' database output files and other initial conditions and then perform replicate simulations as appropriately defined.
#' In essense this is the master wrapper function for all other functions. If you want to test/see SHAPE's default run
#' then simply call this function after loading the library you'll see an experiment built under your root directory.
#' It at least requires that defineSHAPE have been run, else this is going to fail.
#'
#' @param loop_thisRep This is the first replicate value to be simulated in this run, it is standard 1 but can be changed
#' to help with recovery in the middle of a series of replicates.
#' @param workingReplicates This is the maximum replicate number to to simulated in this call.  It is meaningfully different
#' from the number of replicates to be run only when loop_thisRep != 1.
#' @param tmpEnvir_recycleParms This is an environment used to temporarily store loaded RData file objects so that parameters
#' from previous runs, that were stored in RData, can be read back in as required.
#'
#' @examples
#' # First step is to set parameters for the run, this could be done manually but I
#' # recommend using the defineSHAPE function which has a default setting for all
#' # possible parameters and will calculate the value of derived/conditional parameters.
#' defineSHAPE()
#' # Now you can run the simulations, you should get printout to your stdout.
#' \donttest{ runSHAPE() }
#' # Now go and check the SHAPE working directory, which can be found at:
#' getOption("shape_workDir")
#' list.files(getOption("shape_workDir"))
#' # You'll have an experiment folder as well as post-analysis folder
#' # created each with appropriate output!
#' @export
runSHAPE <- function(loop_thisRep = getOption("shape_thisRep"),
                     workingReplicates = seq(getOption("shape_thisRep"),getOption("shape_maxReplicates"),by=1),
                     tmpEnvir_recycleParms = new.env()){
  # This checks if defineSHAPE is likely to have been run
  if(length(which(grepl("^shape_",names(options())))) <
     length(as.list(args(defineSHAPE))) - 1){
    print("Doesn't look like expected options have been set for SHAPE,\n either run defineSHAPE() before trying again or doing the same prior to setting options manually.")
    invisible( NULL )
  }

  # This controls how many, and the number for which, replicates are to be run in this instance.
  for(loop_thisRep in workingReplicates){
    # We set the loop's rep to be the shape_thisRep
    options("shape_thisRep"=loop_thisRep)
    # We now update the processed filename
    options("shape_processedData_fileName" = paste(getOption("shape_outDir"),
                                                   "processed_runData_from_",
                                                   getOption("shape_save_batchString"),
                                                   "_",
                                                   getOption("shape_thisRep"),
                                                   ".RData",sep=""))
    ####### THESE ARE FILENAME AND OBJECT NAME CONSTANTS USED IN THE SCRIPT #######
    # This defines a naming system that provides unique names for each experiment and parameter combination as controlled by the
    # batchString and replicate number information.
    ### CAUTION: NEVER include digits or underscores (ie: "_" ) in this string literal it may break post analysis.
    options("shape_save_batchIndex" = sapply(c("Landscape", "Steps", "Parameters", "Neighbourhood"),function(thisIndex){
                                            paste(getOption("shape_save_batchString"),
                                                  thisIndex,
                                                  if(getOption("shape_run_isRecycling")[thisIndex]){
                                                    getOption("shape_recycle_repStart")
                                                  } else {
                                                    getOption("shape_thisRep")
                                                  },sep=getOption("shape_sepString"))
                                          }))
    # This is the name of the file output for saving a run
    save_fileName <- paste(paste(getOption("shape_save_batchIndex")["Parameters"],sep=getOption("shape_sepString")),".RData",sep="")
    # This is a name used for our database files, they'll be unique to this run.
    options("shape_fileName_dataBase" = list("genotypeSpace"=paste(getOption("shape_save_batchIndex")["Landscape"],".sqlite",sep=""),
                                             "timeStep_States"=paste(getOption("shape_save_batchIndex")["Steps"],".sqlite",sep=""),
                                             "nearestNeighbours"=paste(getOption("shape_save_batchIndex")["Neighbourhood"],".sqlite",sep="")))

    #################################################################################################
    ############################### LANDSCAPE MODEL SPECIFIC INFORMATION ############################
    #################################################################################################
    # If we're recycling our landscape, we must recapture the NK and RMF parameters in order for
    # the runs to be properly interpretable, else we're redefining this.
    if(getOption("shape_run_isRecycling")["Landscape"] &&
       c(getOption("shape_thisRep") > getOption("shape_recycle_repStart") ||
         file.exists(paste(getOption("shape_outDir"),
                           paste(getOption("shape_save_batchString"),
                                 "Parameters",
                                 getOption("shape_recycle_repStart"),
                                 sep = getOption("shape_sepString")),
                           ".RData",sep=""))) ){
      # These store information that would be stored in the loaded file but must be maintained as is for this run
      # these objects will be used downstream to restore the values.
      tmp_realReplicate <- getOption("shape_thisRep")
      tmp_realDist <- getOption("shape_track_distSize")
      tmp_real_batchStrings <- getOption("shape_save_batchIndex")
      tmp_realNames <- getOption("shape_fileName_dataBase")

      # We create an environment , into which we can load the existing information so the needed parts can be obtained
      tmp_tryFile <- paste(getOption("shape_outDir"),
                           paste(getOption("shape_save_batchString"),
                                 "Parameters",
                                 getOption("shape_recycle_repStart"),
                                 sep=getOption("shape_sepString")),
                           ".RData",sep="")
      if(file.exists(tmp_tryFile)){
        load(tmp_tryFile, envir= tmpEnvir_recycleParms )
      } else {
        stop("There was a problem trying to load the seed runs' parameters during call to recycle Landscape.  Please review")
      }
      # Now dependeing on wheter or not we're recycling all parameters or just the fitness landscapes
      # which happens to be the DFE, FitnessLandscape, NK_modelElements, RMF_modelElements, sets....
      # we load different sets of parameters
      loadSet <- if(!getOption("shape_run_isRecycling")["Parameters"]){
        c("FitnessLandscape","DFE","NK_modelElements", "RMF_modelElements")
      } else {
        names(getOption("shape_saveParameters"))
      }

      # we now simply go through the sets to load an place them into the global environment
      for(thisSet in loadSet){
        for(thisNamed in names(tmpEnvir_recycleParms$runParameters[[thisSet]])){
          options(setNames(list(tmpEnvir_recycleParms$runParameters[[thisSet]][[thisNamed]]),thisNamed))
        }
      }

      # Here is where we restore the necessary information back to its initialised values.
      options(list("shape_thisRep" = tmp_realReplicate,
                   "shape_track_distSize" = tmp_realDist,
                   "shape_save_batchIndex" = tmp_real_batchStrings,
                   "shape_fileName_dataBase" = tmp_realNames))
      # We now remove our loaded parameters to save workspace.
      rm(list=ls(, envir = tmpEnvir_recycleParms), envir = tmpEnvir_recycleParms)

    } else {
      ### This section builds fitness landscape model objects for the NK, RMF, Additive and Fixed models ###

      # This is the number of interactions a site has in an NK model, this value is not used unless shape_simModel == "NK", but should be
      # some real positive integer.  The maximum value can be 1 less than the length of the genome.
      options("shape_const_numInteractions" = min(getOption("shape_const_numInteractions"),getOption("shape_genomeLength") - 1))

      # This is the matrix of interactions between sites for the NK model, if we have no interactions among sites, or are not using the NK model
      # we return this as NULL, otherwise we want a matrix with the number of columns equal to the number of interactions
      options("shape_const_NK_interactionMat" = set_const_NK_interactionsMat())

      # This builds the indepednent effect for each site, based on state, where the wild-type values are zeroes
      options("shape_const_siteBystate_fitnessMat" = set_siteByState_fitnessMat())

      # Lastly, as a computational convenience for calculating the per genotype fitenss values.  We calculated the per-site, dependent
      # fitness values of the ancestral genotype.  If this way, for each novel genotype we need only update the value of sites which depend
      # upon a mutant sites
      options("shape_const_DepbySite_ancestFitness" = set_DepbySite_ancestFitness())

      ### These following values are the are only meaningful when the RMF model is implemented

      # This value is the weighting of the independent fitness component > 0, but since in the RMF fitness model equation the distance to the
      # global optima impacts the dynamics so strongly (through c - the independent contribution), and we want to have our dynamics to have some meaningful
      # relationship, we'll base this value on some idealised theta value (which is a measure of the contribution of independent part and the random component).
      # From Neidhart 2014 theta is measured as:  theta = c/(sqrt(var(random_component)))
      options("shape_const_RMF_indWeight" = set_RMF_indWeight())

      # For RMF some optimal genotype exists and here we'll define what that genotype is.  If the model is not RMF, we return NA
      options("shape_const_RMF_globalOptima" = set_const_RMF_globalOptima())

    }

    ###################################################################################################
    ##################################### END OF PREDEFINITIONS #######################################
    ###################################################################################################


    ###################################################################################################
    ######################################### BEGIN OF BODY ###########################################
    ###################################################################################################
    # We want to report on the run time so initialise that reporter here where the simulations are near to starting.
    startTime <- proc.time()

    ######## HERE WE ESTABLISH DATABASE CONNECTIONS, THE FITNESS LANDSCAPES, AND REPORTING OBJECTS #######
    # This creates our SQL Lite database which will hold our growing landscape and reporting objects, this is done since our method
    # may (and is intended to) create large objects which are best not held in R directly.  If the connections exist, this resets them.
    connections_dataBase <- sapply(names(getOption("shape_fileName_dataBase")),function(x){
      reset_shapeDB(paste(getOption("shape_outDir"), getOption("shape_fileName_dataBase")[[x]],sep=""),
                    func_type = "connect")
    })

    ##################################################################################################################
    ########################## NOW WE ESTABLISH THE INITIAL POPULATION AND OTHER OBJECTS #############################
    ##################################################################################################################
    # We initialise the fitness landscape and then the step-wise population tracker database with the initial population.
    # Here is where if we allow a user to redefine the starting population we'd make changes to SHAPE.

    # This is the fitness landscape database
    if(!getOption("shape_run_isRecycling")["Landscape"]  ||
       c(getOption("shape_thisRep") == getOption("shape_recycle_repStart")  &&
         length(RSQLite::dbListTables(connections_dataBase$genotypeSpace)) == 0) ){
      options("shape_tmpGenoTable" = create_genotypeFrame(getOption("shape_nextID"),
                                                          "",
                                                          if(getOption("shape_const_relativeFitness")){
                                                            calc_relativeFitness(getOption("shape_const_ancestFitness"))
                                                          } else {
                                                            getOption("shape_const_ancestFitness")
                                                          }))
      dbWriteTable(connections_dataBase$genotypeSpace,
                   name = if(getOption("shape_db_splitTables")){ nameTable(0,1) }else{ nameTable(0) },
                   value = getOption("shape_tmpGenoTable"))
    }


    # This is the step-wise population tracking database, we initialise a version that will live
    # inside the scope of this replicate's loop, it is a named list object so that we can ensure
    # that we're loading the correct next step, else we'll have to look at the DB
    current_workMatrix <- NULL
    if(!getOption("shape_run_isRecycling")["Steps"]  || getOption("shape_thisRep") == getOption("shape_recycle_repStart")){
      current_workMatrix <- setNames(list(reportPopulations(func_numMuts = 0,
                                                            func_genotypeID = 0,
                                                            func_popSizes = getOption("shape_const_focal_popValue")/
                                                              if(getOption("shape_const_growthForm") == "logistic"){
                                                                getOption("shape_init_distSteps")["factor"]
                                                              }else{
                                                                1
                                                              },
                                                            func_fitnesses = if(getOption("shape_const_relativeFitness")){
                                                              calc_relativeFitness(getOption("shape_const_ancestFitness"))
                                                            } else {
                                                              getOption("shape_const_ancestFitness")
                                                            },
                                                            func_births = getOption("shape_const_focal_popValue")/
                                                              if(getOption("shape_const_growthForm") == "logistic"){
                                                                getOption("shape_init_distSteps")["factor"]
                                                              }else{
                                                                1
                                                              },
                                                            func_deaths = 0,
                                                            func_mutants = 0,
                                                            func_progenitor = "WT")),
                                     paste("Step",0,sep=getOption("shape_sepString")))
      # We name our steps using the step naming function, and setup our population matrix with the population reporting function.
      dbWriteTable(connections_dataBase$timeStep_States,
                   name = names(current_workMatrix),
                   value = current_workMatrix[[1]],
                   overwrite = TRUE)
    }
    # We initialise the object for tracking disturbance and timing.
    options("shape_track_distSize" = rbind(c("Step"=0,getOption("shape_init_distSteps"))))

    ##################################################################################################################
    ############################ RESSETTING ANCESTRAL FITNESS VALUE FOR RMF MODEL ####################################
    ##################################################################################################################
    # Now if using the RMF model, we reset the ancestral fitness using the fitness calculation
    # But we have to check if we're recycling our information here, otherwise we'll have loaded this
    if(!getOption("shape_run_isRecycling")["Landscape"]  ||
       c(getOption("shape_thisRep") == getOption("shape_recycle_repStart")  && !file.exists(paste(getOption("shape_outDir"), save_fileName,sep=""))) ){
      # For the RMF model we'll want the ancestral fitness to be the weighted distance component plus
      # what has been set as the shape_const_ancestFitness value.  We get the distance compoenent only by fudging the distribution
      # as being a null zero distribution of values.  The assumption here is that "shape_const_ancestFitness" was originally set to
      # reflect the upper end of the distribution passed for future use.
      if(getOption("shape_simModel") == "RMF"){
        # I assign in here so that the ancestral value is based solely on the distance and not the distribution
        # this makes our calculations of relative fitness, in RMF landscape models, feasible
        # -- as magnitude was found to at times (when near 0) to be a problem to work with
        # We also suppres the distAsS so that the uniform zero isn't changed to a -1 value.
        options("shape_const_ancestFitness" = fitnessLandscape("",
                                                               getOption("shape_const_ancestFitness"),
                                                               landscapeModel = getOption("shape_simModel"),
                                                               relativeFitness = FALSE,
                                                               func_distribution = "Uniform",
                                                               func_distParameters = c(0,0),
                                                               func_distAsS = FALSE))
      }
    }

    # This saves the parameters object file
    runTime <- 0
    if(!getOption("shape_run_isRecycling")["Parameters"] ||
       c(getOption("shape_thisRep") == getOption("shape_recycle_repStart") &&
         !file.exists(paste(getOption("shape_outDir"),getOption("shape_save_batchIndex")["Parameters"],".RData",sep=""))) ){
      # This builds a list of the parameter objects we want saved so the values can easily be recovered later.
      runParameters <- sapply(getOption("shape_saveParameters"), function(x){
        sapply(x,function(y){
          getOption(y)
        },simplify=FALSE)
      },simplify=FALSE)
      # We call this here since I can't merge I below elsewise that I know of
      shape_fileName_dataBase <- getOption("shape_fileName_dataBase")
      save(runParameters, runTime, connections_dataBase, shape_fileName_dataBase,
           file=paste(getOption("shape_outDir"), save_fileName,sep=""))
      # If we are recycling parameters and we've already saved the repStart parameter file, we'll just load all that information into the global space.
    } else if(getOption("shape_run_isRecycling")["Parameters"] &&
              file.exists(paste(getOption("shape_outDir"),getOption("shape_save_batchIndex")["Parameters"],".RData",sep="")) ){
      # Similar to previous instance of loading, this saves some objects to be restored after the load.
      tmp_realReplicate <- getOption("shape_thisRep")
      tmp_realDist <- getOption("shape_track_distSize")
      tmp_real_batchStrings <- getOption("shape_save_batchIndex")
      tmp_realNames <- getOption("shape_fileName_dataBase")

      # We create an environment, load the values, assign what is needed
      tmp_loadFile <- paste(getOption("shape_outDir"),
                            paste(getOption("shape_save_batchString"),
                                  "Parameters",
                                  getOption("shape_recycle_repStart"),
                                  sep=getOption("shape_sepString")),
                            ".RData",
                            sep="")
      if(file.exists(tmp_loadFile)){
        load(tmp_loadFile, envir= tmpEnvir_recycleParms )
      } else {
        stop("There was a problem trying to load the seed runs' parameters during call to recycle Landscape.  Please review")
      }
      # Even if we're not re-cycling the specific parameters with a logical call, we may need to recover the fitness landscape values.
      loadSet <- if(!getOption("shape_run_isRecycling")["Parameters"]){
        c("FitnessLandscape","DFE","NK_modelElements", "RMF_modelElements")
      } else {
        names(getOption("shape_saveParameters"))
      }
      # We now simply go through the sets to load an place them into the global environment
      for(thisSet in loadSet){
        for(thisNamed in names(tmpEnvir_recycleParms$runParameters[[thisSet]])){
          options(setNames(list(tmpEnvir_recycleParms$runParameters[[thisSet]][[thisNamed]]), thisNamed))
        }
      }
      # We re-update our real replicate value - this would only matter if shape_run_isRecycling["Parameters"] == FALSE
      options(list("shape_thisRep" = tmp_realReplicate,
                   "shape_track_distSize" = tmp_realDist,
                   "shape_save_batchIndex" = tmp_real_batchStrings,
                   "shape_fileName_dataBase" = tmp_realNames))
      # We now remove our loaded parameters to save workspace.
      rm(list=ls(, envir = tmpEnvir_recycleParms), envir = tmpEnvir_recycleParms)
    }
    # At this point we check what was the max previously recorded genotypeID and pass this to the run
    options("shape_nextID" = max(unlist(lapply(RSQLite::dbListTables(connections_dataBase$genotypeSpace), function(thisTable){
      dbGetQuery(connections_dataBase$genotypeSpace,
                 paste("SELECT MAX(genotypeID) FROM ", thisTable,sep=""))
    }))) + 1)

    # We now close the connections to avoid database malformations, recall it was opened to initialise the population.
    for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }

    # We initialise our disturbance counter and a starting step object that is not currently usefull outside of troubleshooting.
    startingStep <- 1
    options("shape_counter_logSteps" = startingStep)

    ################################################################################################################################
    ################################################### ACTUAL ITERATIVE RUN BODY ##################################################
    ################################################################################################################################
    # Now we run through the steps of evolution, this loop is wherein we find the real magic.
    for(thisStep in startingStep:(getOption("shape_numGenerations") / getOption("shape_size_timeStep"))){
      # This is simply for reporting reasons, will be found in the stdout
      if(thisStep %% 100 == 0){
        print(paste("We're now on thisStep: ",thisStep,sep=""))
        print( proc.time() - startTime )
      }
      # Now, if the last step's states are not in the current_workMatrix, then we load it from the DB file
      # NOTE: This is expected to occur only when recovering from mid run crashes re-starts.
      if(!is.element(nameTable_step(thisStep - 1),names(current_workMatrix))){
        # We open the connections, load the population states of the last step, then close the connection.
        connections_dataBase <- sapply(names(getOption("shape_fileName_dataBase")),function(x){
          reset_shapeDB(paste(getOption("shape_outDir"), getOption("shape_fileName_dataBase")[[x]],sep=""),
                        func_type = "connect")
        })
        # We load the SHAPE community state of the last step
        current_workMatrix <- c(current_workMatrix,
                                setNames(list(dbReadTable(connections_dataBase$timeStep_States, nameTable_step(thisStep - 1))),
                                         nameTable_step(thisStep - 1)))
        # We close the connections to the DB to avoid corruption. (it's been found to happen)
        for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
      }
      # This runs one step of a SHAPE simulation
      current_workMatrix <- runReplicate(func_inputFrames = current_workMatrix,
                                         func_currStep = thisStep,
                                         func_stepCounter = getOption("shape_counter_logSteps"))

      # We now advance our step counter to track the disturbance calls.
      options("shape_counter_logSteps" = getOption("shape_counter_logSteps") + 1)
    }

    ###################################################################################################
    ########################################## END OF BODY ############################################
    ###################################################################################################



    ###################################################################################################
    ################################### BEGIN OF SUMMARY TOOLS ########################################
    ###################################################################################################

    # Track the run time value so it can be reported.
    runTime <- proc.time() - startTime

    # This saves the parameters, again, only difference here is that the runTime is now completed.
    # I admit this step is probably not required and if anything may be a confusion since the file
    # is already saved prior to run BUT now the proper run time is saved...
    if(!getOption("shape_run_isRecycling")["Parameters"] ||
       getOption("shape_thisRep") == getOption("shape_recycle_repStart")){
      # This is the list that holds the parameter information
      runParameters <- sapply(getOption("shape_saveParameters"), function(x){
        sapply(x,function(y){
          getOption(y)
        },simplify=FALSE)
      },simplify="list")
      shape_fileName_dataBase <- getOption("shape_fileName_dataBase")
      save(runParameters, runTime, connections_dataBase, shape_fileName_dataBase,
           file=paste(getOption("shape_outDir"), save_fileName,sep=""))
    }

    connections_dataBase <- sapply(names(getOption("shape_fileName_dataBase")),function(x){
      reset_shapeDB(paste(getOption("shape_outDir"), getOption("shape_fileName_dataBase")[[x]],sep=""),
                    func_type = "connect")
    })
    # Now for all unexplored mutational spaces around existing genotypes, we define that neighbourhood.
    tmpTables <- RSQLite::dbListTables(connections_dataBase$timeStep_States)
    all_lastGenotypes <- dbReadTable(connections_dataBase$timeStep_States,
                                     RSQLite::dbListTables(connections_dataBase$timeStep_States)[which.max(as.numeric(nameTable_step(tmpTables, funcSplit = TRUE)))])
    # We check if there are any established lineages that need to have their local neighbourhood defined.
    establishedGenotypes <- querryEstablished(func_inMatrix= all_lastGenotypes,
                                              func_estProp = getOption("shape_const_estProp"))
    if(nrow(establishedGenotypes) > 0){
      # If a mutant is not found within the genotypeID's of the database of stepChanges,
      # then we need to create is local nieghbourhood so that our final reporting can consider that space.
      for(thisGenotype in 1:nrow(establishedGenotypes)){
        tmp_genotypeInfo <- retrieve_binaryString(func_genotypeID = establishedGenotypes[thisGenotype,"genotypeID"],
                                                  func_numMuts = establishedGenotypes[thisGenotype,"numMuts"],
                                                  func_subNaming = getOption("shape_db_splitTables"),
                                                  func_landscapeCon = connections_dataBase$genotypeSpace)
        # We ask if this genotype has been explored
        if(!as.logical(tmp_genotypeInfo[1,"isExplored"])){
          createGenotypes(tmp_focalGenotype = tmp_genotypeInfo[1,"binaryString"],
                          tmp_focalFitness= if(is.element(getOption("shape_simModel"),c("NK","Additive"))){
                            getOption("shape_const_siteBystate_fitnessMat")
                          } else {
                            establishedGenotypes[thisGenotype,"fitness"]
                          },
                          maxHamming = getOption("shape_max_numMutations"),
                          tmp_landModel = getOption("shape_simModel"),
                          tmp_relativeFitness = getOption("shape_const_relativeFitness"),
                          tmpDistribution = getOption("shape_constDist"),
                          tmpParameters = getOption("shape_const_distParameters"),
                          tmp_genCon = connections_dataBase$genotypeSpace,
                          tmp_tableSplit = getOption("shape_db_splitTables"),
                          tmp_distAsS = getOption("shape_const_distAsS"))
        }
      }
    }
    for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }

    # This is a reporting call, it will be in the stdout.
    if(getOption("shape_toggle_forceCompletion")  &&
        thisStep != (getOption("shape_numGenerations") / getOption("shape_size_timeStep"))){
        stopError("There was a problem and the main body run did not complete, please review")
    } else if(thisStep == (getOption("shape_numGenerations") / getOption("shape_size_timeStep"))) {
        print("Simulations completed without crashing ")
    } else {
        print("Simulation code reached end of loop ")
    }

    ######### Now the results of the run are processed. ##############
    connections_dataBase <- sapply(names(getOption("shape_fileName_dataBase")),function(x){
      reset_shapeDB(paste(getOption("shape_outDir"),
                          getOption("shape_fileName_dataBase")[[x]],
                          sep=""),
                    func_type = "connect")
    })

    # Now we run the processing function noting that the <shape_processedData_fileName> and <shape_processedObjects>
    # are actually loaded from the source file's parameters
    tmpReport <- runProcessing(func_saveFile = getOption("shape_processedData_fileName"),
                               func_subNaming = getOption("shape_db_splitTables"),
                               func_stepsCon = connections_dataBase[["timeStep_States"]],
                               func_landscapeCon = connections_dataBase[["genotypeSpace"]],
                               func_hoodCon = connections_dataBase[["nearestNeighbours"]],
                               func_estProp = getOption("shape_const_estProp"),
                               func_size_timeStep = getOption("shape_size_timeStep"),
                               func_processObjects = getOption("shape_processedObjects"),
                               func_hoodPriority = getOption("shape_const_hoodDepth"))
    # This gets printed to stdout and is simply a message reporting on if processing occured properly
    print( tmpReport )
    # Now that we've pre-processed we delete the Steps file if asked to do so, but only if we've been
    # given a flag that the processing was completed
    if(getOption("shape_results_removeSteps") &&
       grepl("Processing completed", tmpReport) &&
       file.exists(getOption("shape_processedData_fileName"))){
      # If all those conditions are true, we have every expectation the the run was proccessed correctly.
      # The purpose of deleting this is that it can require a lot of disk space and once processed
      # is not likely to contain information of interest for your experiment.
      dbDisconnect(connections_dataBase[["timeStep_States"]])
      connections_dataBase[["timeStep_States"]] <- NULL
      file.remove(paste(getOption("shape_outDir"),
                        getOption("shape_fileName_dataBase")[["timeStep_States"]],
                        sep=""))
    }
    for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
    # Now if we're not recycling our neighbourhood database we'll be removing it from existence
    if(!getOption("shape_run_isRecycling")["Neighbourhood"]){
      file.remove(paste(getOption("shape_outDir"),
                        getOption("shape_fileName_dataBase")[["nearestNeighbours"]],
                        sep=""))
    }

    ###################################################################################################
    #################################### END OF SUMMARY TOOLS #########################################
    ###################################################################################################



    ###################################################################################################
    ################################### BEGIN OF SELFING TOOLS ########################################
    ###################################################################################################
    # Normally, you'd run all replicates within a single call and there would be no need to call
    # new R runs.  However, when I was developing SHAPE, the remote server with which I was working
    # had a maximum wall time for any single proccess/job.  This meant it became important that I develop
    # a means for each replicate of SHAPE runs to be independent.  This section will be important
    # if you've identified a similar circumstance where you need to create independent jobs in order to
    # avoid wall-time issues.  It is namely controlled by the logical toggle "shape_externalSelfing"
    if(getOption("shape_thisRep") < getOption("shape_maxReplicates") && getOption("shape_externalSelfing")){
      # Then we read in the selfing template script, but we check that it exsits.
      if(!file.exists(getOption("shape_tmp_selfScript"))){
        stopError("There was a problem while trying to self, as the template did not exist.  Please review")
      } else {
        # Now it depends on if we're selfing on the remote CAC server or doing a local job
        # NOTE: It's assumed that the shape_tmp_selfScript is keyed for the proper type of job
        # This is the name for the job we want to create
        tmp_nextRep <- getOption("shape_thisRep") + 1
        tmp_jobName <- name_batchString(funcBase = getOption("shape_save_batchString"),
                                        func_repID = tmp_nextRep)
        # I simply update the selfing scripts
        tmpScript <- readLines(getOption("shape_tmp_selfScript"))
        # We update the getOption("shape_thisRep") information and write out the file
        tmpScript <- updateLines(func_inLines = tmpScript,
                                 func_searchPattern = list("rep"=c("shape_thisRep=","shape_thisRep=([[:digit:]])+")),
                                 func_values = list("rep"=paste("shape_thisRep=",tmp_nextRep,sep="")))
        # We write out the lines and make the script executable
        writeLines(tmpScript, con = getOption("shape_tmp_selfScript"))
        Sys.chmod(getOption("shape_tmp_selfScript"), mode="0777")
        # Now, this matters only if the job that is currently runing has some code such that
        # once it is complete, it will submit a new job and pass in the outside argument of the replicate.
      }
      # To continue the run we also remove the stop trigger, because as per the selfing script I've written
      # if will continue to submit new jobs until the stopFile exists at the expected location.
      tmp_stopFile <- paste(getOption("shape_finalDir"),getOption("shape_external_stopFile"), sep="")
      if(file.exists(tmp_stopFile)){
        file.remove(tmp_stopFile)
      }
    }
    # If we're selfing internally, but using a remote server to farm, where the compute nodes
    # are not trust worthy to hold the intermediate files indefinitely (which was my case and
    # I would presume yours as well - PS don't have the run call from the compute node to DB
    # files on the host node as the latency is BRUTAL), then this will ensure that the compute
    # node's locally created files are regularly passed back to a relative filepath on the host node.
    if(getOption("shape_serverFarm") && !getOption("shape_externalSelfing")) {
      # Every twenty replicates we'll make a call to transfer out files back to the original directory of the run
      if(getOption("shape_thisRep") %% 20 == 0){
        system(paste("cd ",getOption("shape_outDir")," ; cp *.RData *.sqlite *.Rout *.o ",getOption("shape_finalDir"),sep=""))
      }
    }

    ###################################################################################################
    #################################### END OF SELFING TOOLS #########################################
    ###################################################################################################

  } # This closes the for loop of potential internal replicating as define by the vector workingReplicates
  # and the logical toggle shape_externalSelfing
  invisible( NULL )
}



#' This is a function used to read the SHAPE_experimentalDesign type input file
#' and then build a SHAPE experiment by creating all the folder structure, .R and .sh scripts
#' required to programatically run your experiment -- excluding post-analysis, that's a you problem.
#'
#' @param func_filepath_toDesign This is the absolute filepath which points to the SHAPE_experimentalDesign like
#' template you'd like used to identify parameter combinations for building your experiment.
#' @param func_templateDir This is the absolute filepath to a directory on your machine where the SHAPE template
#' scripts/files have been saved.  They are used by this function to help build your experiment.
#' @param func_maxGrouped_perShell Integer value defining the maximal number of jobs that an output shell script
#' will try to have run in parallel once executed.  This is related to your parallel computing potential.
#' @param func_filePath_R This is the absolute path to the R application on the system where SHAPE would be run
#' via BATCH MODE, its value is applied in shell scripts written for running the experiment.  If left NULL then
#' this function will try to use standard R install paths of which I'm aware.
#' @param func_baseCall This is a string element of arguments when calling BATCH MODE if R via shell script.
#' @param func_rArgs This is a character string which represent additional arguments to be passed via shell
#' script BATCH mode call of R.  I consider it most practicable to set the replicate and output directory of SHAPE.
#' @param func_remoteLocation The filepath of the compute node on a remote server where your job would be run.  The default
#' is based on the environment variable value used in CAC's SLURM submission system.
#' @param func_submitArgs This is information concerning sheel script lines for automatic submission of jobs to
#' a remote server's submission system.  I'm basing this off of the SLURM system of the Center for Advanced Computing
#' Queen's University computing platform.  If your system is different you may need to tweak this.  Sorry?
#' This should be a vector of arguments passed for job submissions on a remote server
#' The example here would call 1 core with 8 Gb RAM and a wall time of 14 days and an outFile be named
#' You can add more arguments if your server requires this, they'll get used.
#' BUT where the job's name MUST be identified as ---  fakeJob   ----  and the output log as  --- fakeOut ---, you can change the argument queues
#' I also assume your remote server will create a local directory on the compute nodes whre your job once submitted,
#' and that there will be the location defined by func_remoteLocation.
#' @param func_processingCores This is the number of parallel cores you would like the summairseExperiment() to call
#' when trying to process your experimental output.
#' @param func_suppressOld_summaryFiles Logical flag controlling if your summariseExperiment() will delete old output
#' summary files.  setting to FALSE (default) is ideal if you could ever expect you might need to restart whereas
#' TRUE becomes practical if you are worried you'd have updated output to process and you want to ensure a fresh processing start.
#'
#' @return If no error is encountered, a message will be returned suggesting the build was successful.  SHAPE makes
#' no effort to perform validation of this effort to build the experiment and presumes no fatal errors is sufficient evidence.
#'
#' @examples
#' # This function relies on script templates which can be found at:
#' # 'https://github.com/JDench/SHAPE_library/tree/master/SHAPE_templates'
#' # Once these have been downloaded you can pass the appropriate filepath values
#' # to the first two arguments. For this example, I'll assume you've installed
#' # them to a folder position that is now just under the root of your
#' # R-environment working directory.
#' # However, before runing the function we need to parameterise your run of SHAPE,
#' # here I call the default parameters:
#' defineSHAPE()
#' # Now using the default templates we design an experiment folder complete with
#' # shell scripts to submit our work programatically.
#' # NOTE: Again, this example assumes you've downloaded the templates and placed
#' #        them at the next filepath and directory-path locations
#' \donttest{shapeExperiment(func_filepath_toDesign = "~/SHAPE_templates/SHAPE_experimentalDesign.v.1.r",
#'                              func_templateDir = "~/SHAPE_templates/")}
#' # You should be greeted with a message suggesting your experiment is built.
#' # You can find the files now at that script's SHAPE workingDirectory.
#' list.files(getOption("shape_workDir"))
#' # Voila!  You can go see the spread of variable evolutionary parameters that were
#' # considered by looking at -- yourJob_parameterCombos.table -- which is a tab
#' # delimated file.
#' # Lastly, you may have R installed elsewhere and so want to have that noted while
#' # your experiment is built because the shell scripts will need to point to the correct place.
#' \donttest{shapeExperiment(func_filepath_toDesign = "~/SHAPE_templates/SHAPE_experimentalDesign.v.1.r",
#'                              func_templateDir = "~/SHAPE_templates/",
#'                              func_filePath_R = "~/your_R_folder/R_app/bin/R")}
#' # Now obviously the above location likely is not where you installed R,
#' # but ideally you get the point. The difference is in how the shell scripts were written.
#'
#' @export
shapeExperiment <- function(func_filepath_toDesign, func_templateDir,
                            func_maxGrouped_perShell = 2,
                            func_filePath_R = NULL, func_baseCall = "CMD BATCH",
                            func_rArgs = '"--args shape_thisRep=1 shape_outDir=\'fake_serverPath/fakeDir/\'"',
                            func_remoteLocation = "$TMPDISK",
                            func_submitArgs = c("number_ofCores"='-c 1',
                                                 "memory"='--mem=8192',
                                                 "jobName"='-J fakeJob',
                                                 "wallTime"='-t 14-00:00:00',
                                                 "fileOut"='-o fakeOut'),
                            func_processingCores = 1,
                            func_suppressOld_summaryFiles = FALSE){
  #### DID YOU GRAB THE SHAPE TEMPLATE FILES FROM:
  #### IF NOT, THEN THIS ISN'T GOING TO WORK, IF SO THEN THE FOLDER
  #### TO WHERE YOU SAVED THEM SHOULD BE THE ARGUMENT:  func_templateDir
  #### AND YOUR func_filepath_toDesign ARGUMENT FILE WILL BE A MODIFIED
  #### (TO REPRESENT YOUR EXPERIMENT )FORM OF THE FILE FROM THOSE TEMPLATES.

  # We check that the design template file exists, otherwise this won't matter
  if(!file.exists(func_filepath_toDesign)){
    stopError(paste("Could not find the experimental design file at: ",func_filepath_toDesign,sep=""))
  }
  # If not specified, then we look for R installed in common locations
  if(is.null(func_filePath_R)){
    func_tmpVersion <- strsplit(version[['version.string']], "[[:space:]]+")[[1]][3]
    # These are the standard install filepaths as far as I know them
    func_filePath_R <- if(grepl("^window",Sys.info()["sysname"],ignore.case=TRUE)){
                          paste("C:/Program Files/R/R-",func_tmpVersion,"/bin/R",sep="")
                       } else {
                          "/usr/bin/R"
                       }
  }

  # This is a reference object that links that farmerParms input names to the object names used in SHAPE
  inputReference <- list("shape_serverFarm",
                         "shape_results_removeSteps",
                         "shape_externalSelfing",
                         "shape_toggle_forceCompletion",
                         "shape_workDir",
                         "shape_postDir",
                         "shape_save_batchBase",
                         "shape_save_batchJob",
                         "shape_maxReplicates",
                         "uniqueReplicates",
                         "shape_const_distType",
                         "shape_init_distPars",
                         "shape_const_growthGenerations",
                         "shape_numGenerations",
                         "shape_size_timeStep",
                         "shape_genomeLength",
                         "shape_const_focal_popValue",
                         "shape_const_birthProb",
                         "shape_const_growthRate",
                         "shape_const_deathProb",
                         "shape_death_byDensity",
                         "shape_death_densityCorrelation",
                         "shape_death_densityCap",
                         "shape_scaleGrowth_byDeaths",
                         "shape_const_mutProb",
                         "shape_muts_onlyBirths",
                         "shape_allow_backMutations",
                         "shape_const_growthForm",
                         "shape_simModel",
                         "shape_const_ancestFitness",
                         "shape_constDist",
                         "shape_const_distParameters",
                         "shape_const_distAsS",
                         "shape_const_RMF_initiDistance",
                         "shape_const_RMF_theta",
                         "shape_const_numInteractions",
                         "shape_const_fixedFrame",
                         "shape_includeDrift",
                         "shape_track_asWhole",
                         "shape_const_estProp",
                         "shape_const_hoodThresh")

  # This is a reference of the input parameters that are not used in making factorial combinations for an experiment
  notExpanded_reference <- c("shape_serverFarm",
                             "shape_results_removeSteps",
                             "shape_externalSelfing",
                             "shape_toggle_forceCompletion",
                             "shape_workDir",
                             "shape_save_batchBase",
                             "shape_maxReplicates",
                             "uniqueReplicates")

  # This is a reference of those run parameters that used 1:1 to make a parameter set.
  comboReference <- c("shape_simModel",
                      "shape_const_ancestFitness",
                      "shape_constDist",
                      "shape_const_distParameters",
                      "shape_const_distAsS")

  # This is a reference of parameters that only need to be made factorial given particular conditions
  conditionReference <- list("RMF"=c("shape_const_RMF_theta",
                                     "shape_const_RMF_initiDistance"),
                             "NK"=c("shape_const_numInteractions"),
                             "Fixed"=c("shape_const_fixedFrame"))
  # These are the string literal expressions that can help this script find the appropriate templates in the func_templateDir folder
  fileName_templates <- c("body"="SHAPE_runBody.v.",
                          "serverSubmit"="SHAPE_serverTemplate.v.",
                          "localSubmit"="SHAPE_localTemplate.v.",
                          "parameter_output" = "SHAPE_parameters.v.")
  # Step 1: find the template files and report any that are missing
  allFiles <- list.files(path = func_templateDir, recursive = FALSE)
  for(thisFile in names(fileName_templates)){
    tmpFile <- allFiles[which(grepl(fileName_templates[thisFile], allFiles,fixed=TRUE))]
    # We check if there is a unique file
    if(length(tmpFile) == 1){
      fileName_templates[thisFile] <- paste(func_templateDir,tmpFile,sep="")
      if(thisFile == "functions"){ source(fileName_templates[thisFile]) }
    } else {
      stopError(paste("Could not find a unique file like ",fileName_templates[thisFile],
                      " in ",func_templateDir," please review and ensure the template exists",sep=""))
    }
  }
  # Step 2: read the input parameters and begin by building the experiment's folders
  inputParms <- readLines(func_filepath_toDesign)
  inputParms <- inputParms[intersect(intersect(which(!grepl("#",inputParms,fixed=TRUE)),
                                               which(!grepl("[[:space:]]",substr(inputParms,
                                                                                 sapply(nchar(inputParms),function(x){ min(1,x) }),
                                                                                 sapply(nchar(inputParms),function(x){ min(1,x) })) ))),
                                     which(nchar(inputParms) > 0))]
  # We now parse the input into values, this involves supressing all leading and trailing spaces.
  inputParms <- lapply(strsplit(inputParms," <-"),function(x){
                    # There is one condition where we split by c(), the rest by commas
                    return( c(x[1],
                              if(!is.element(x[1],c("shape_const_distParameters","shape_init_distPars"))){
                                tmpReturn <- unlist(strsplit(x[2],',',fixed=TRUE))
                                # We now go and try to evaluate
                                tmpReplace <- suppressWarnings(as.numeric(tmpReturn))
                                tmpReturn[which(!is.na(tmpReplace))] <- tmpReplace[which(!is.na(tmpReplace))]
                                tmpSpaces <- gregexpr("[[:space:]]",tmpReturn)
                                for(thisItem in 1:length(tmpSpaces)){
                                  if(tmpSpaces[[thisItem]][1] != -1){
                                    # This removes leading spaces
                                    tmpStart <- if(tmpSpaces[[thisItem]][1] == 1){
                                      tmpMax <- which(diff(tmpSpaces[[thisItem]]) != 1)
                                      if(length(tmpMax) == 0){
                                        tmpMax <- 0
                                      } else {
                                        tmpMax <- tmpSpaces[[thisItem]][max(tmpMax)] - 1
                                      }
                                      2+tmpMax
                                    } else {
                                      1
                                    }
                                    # This removes trailing spaces
                                    tmpEnd <- if(tmpSpaces[[thisItem]][length(tmpSpaces[[thisItem]])] == nchar(tmpReturn[thisItem])){
                                      tmpMax <- which(diff(tmpSpaces[[thisItem]]) != 1)
                                      if(length(tmpMax) == 0){
                                        tmpMax <- nchar(tmpReturn[thisItem])
                                      } else {
                                        tmpMax <- tmpSpaces[[thisItem]][max(tmpMax)+1]
                                      }
                                      tmpMax -1
                                    } else {
                                      nchar(tmpReturn[thisItem])
                                    }
                                    tmpReturn[thisItem] <- substr(tmpReturn[thisItem],tmpStart,tmpEnd)
                                  }
                                }
                                tmpReturn
                              } else {
                                tmpReturn <- unlist(strsplit(gsub("[[:space:]]","",x[2]),'c(',fixed=TRUE))
                                tmpReturn <- tmpReturn[which(nchar(tmpReturn) > 0)]
                                unlist(lapply(tmpReturn,function(x){
                                  paste('c(',
                                        if(substr(x,nchar(x),nchar(x)) == ','){
                                          substr(x,1,nchar(x)-1)
                                        } else {
                                          x
                                        },
                                        sep="")
                                }))
                              } ))
                  }) # This closes out my parsing the inputParameters lines
  # I now name the inputParameters given the first item in their list
  names(inputParms) <- unlist(lapply(inputParms,function(x){x[1]}))
  for(thisList in 1:length(inputParms)){ inputParms[[thisList]] <- inputParms[[thisList]][-1] }
  # Now the inputParms of I've allowed the user to input "targetNumber" if they want the "const_focal_popValue"
  if(is.null(inputParms$shape_death_densityCap) || inputParms$shape_death_densityCap == "NULL"){
    inputParms$shape_death_densityCap <- inputParms$shape_const_focal_popValue
  }
  # Now also, if is.null(inputParms[["shape_workDir"]]), we set to tempdir()
  if(is.null(inputParms[["shape_workDir"]]) || inputParms[["shape_workDir"]] == "NULL"){
    inputParms[["shape_workDir"]] <- paste('"',gsub('\\','/',tempdir(),fixed=TRUE),'"',sep="")
  }

  # This ensures that file pathing objects have terminal slashes
  for(thisObject in c("shape_workDir")){
    tmpString <- inputParms[[thisObject]]
    if(!grepl('/\"$',tmpString)){
      inputParms[[thisObject]] <- paste(substr(tmpString,1,nchar(tmpString)-1),
                                        "/",
                                        substr(tmpString,nchar(tmpString),nchar(tmpString)),
                                        sep="")
    }
  }
  # Now we assess if this is a server run and build the appropriate directories, getting the true path requires
  # removing the leading and trailing escape character and quotation marks
  tmpDir <- paste(trimQuotes(inputParms[["shape_workDir"]]),collapse="")
  if(!dir.exists(tmpDir)){
    dir.create(tmpDir, recursive = TRUE)
  }

  # This is the point where I build a matrix of the parameter combinations
  parameterCombos <- shapeCombinations(func_inLines = inputParms,
                                       func_comboRef = comboReference,
                                       func_indepRef = notExpanded_reference,
                                       func_condRef = conditionReference)
  # We write out the parameter combos as a tab  separated file
  write.table(parameterCombos,
              file = paste(tmpDir,trimQuotes(inputParms$shape_save_batchBase),"_parameterCombos.table",sep=""),
              append = as.numeric(parameterCombos[1,"shape_save_batchJob"]) != 1,
              quote = FALSE,
              row.names = FALSE,
              qmethod = "double",
              sep = "\t",
              dec =".",
              fileEncoding = "cp1252")

  # Now for each possible parameter combination, and unique statistical replicate, we build directories
  tmp_newDirs <- paste(tmpDir,
                       name_batchString(funcBase = trimQuotes(inputParms$shape_save_batchBase),
                                        func_setID = unlist(lapply(1:as.numeric(inputParms$uniqueReplicates),rep,"times"=nrow(parameterCombos))),
                                        func_jobID = parameterCombos$shape_save_batchJob,
                                        func_sepString = getOption("shape_sepString")),
                       sep="")
  for(thisDir in tmp_newDirs){
    if(!dir.exists(thisDir)){ dir.create(thisDir, recursive = TRUE) }
  }

  # I now copy in templates for all template files, updating the plotting and parameter output as required
  for(thisTemplate in names(fileName_templates)){
    if (thisTemplate == "parameter_output") {
      writeParameters(func_infile = fileName_templates[thisTemplate],
                      func_inParms = inputParms,
                      func_inCombos = parameterCombos,
                      func_outDir = tmpDir,
                      func_bodyScript = fileName_templates["body"],
                      func_ExternalStopper = getOption("shape_external_stopFile"))
    } else if(as.logical(inputParms$shape_serverFarm) && thisTemplate == "serverSubmit" ||
              !as.logical(inputParms$shape_serverFarm) && thisTemplate == "localSubmit") {
      write_subScript(func_subScipt = fileName_templates[thisTemplate],
                      func_outDir = tmpDir,
                      func_inCombos = parameterCombos,
                      func_inParms = inputParms,
                      func_maxJobs= func_maxGrouped_perShell,
                      func_appLocation = func_filePath_R,
                      func_remoteLocation = func_remoteLocation,
                      func_commonArgs = func_baseCall,
                      func_submitArgs = func_submitArgs,
                      func_passedArgs = func_rArgs)
    }
  }

  # We write this simple R script to run post analysis, which is basically a setup the the minimum necessary
  # rSHAPE environmental variables for the experiment and the to run the function.
  func_processingScript <- paste(trimQuotes(inputParms[["shape_workDir"]]),"summariseExperiment_",
                                 trimQuotes(inputParms[["shape_save_batchBase"]]),".r",sep="")
  cat(c("library(rSHAPE)",
        paste("defineSHAPE(shape_save_batchSet = ",if(!is.null(inputParms$uniqueReplicates)){1}else{NULL},",\n ",
              "shape_save_batchJob = ",if(!is.null(inputParms$shape_save_batchJob)){1}else{NULL},",\n ",
              'shape_sepString = "',getOption("shape_sepString"),'",\n ',
              'shape_sepLines = "',getOption("shape_sepLines"),'",\n ',
              "shape_workDir = ",inputParms[["shape_workDir"]],",\n ",
              'shape_postDir = "',paste(trimQuotes(inputParms[["shape_workDir"]]),
                                        trimQuotes(inputParms[["shape_postDir"]]),
                                        sep=""),'",\n ',
              "shape_save_batchBase = ",inputParms$shape_save_batchBase,",\n ",
              'shape_string_lineDescent = "',getOption("shape_string_lineDescent"),'")',
              sep=""),
        paste("summariseExperiment(func_numCores = ",func_processingCores,", ",
              "func_suppressOld = ",func_suppressOld_summaryFiles,
              ")",sep=""),
        'q(save="no")'),
      file = func_processingScript,
      sep = "\n",
      append = FALSE)
  # Now we write a shell submission script for the experiment's post analysis
  func_processing_subScript <- paste(trimQuotes(inputParms[["shape_workDir"]]),"summariseExperiment_",
                                     trimQuotes(inputParms[["shape_save_batchBase"]]),".sh",sep="")
  cat(c("#!/bin/bash",
        paste(func_filePath_R,
              func_baseCall,
              func_processingScript,
              paste(func_processingScript,".Rout",sep=""),
              sep=" ")),
      file = func_processing_subScript,
      sep = "\n",
      append = FALSE)

  # I make all .sh files executable
  Sys.chmod(list.files(path=tmpDir,pattern='.sh',recursive=TRUE,full.names=TRUE),
            mode="0777")

  # I now save an image of this workSpace so that future scripts can load into what has been built here
  fileName_farmingImage <- paste(trimQuotes(inputParms$shape_save_batchBase),"_shapeExperiment_image.RData",sep="")
  save.image(paste(tmpDir,fileName_farmingImage,sep=""))

  return( "Function completed to end, experiment is likely built!" )
}


# This is a function for writting out a job's submission script, it allows a submission script to be built for
# either remote or local job submits.
#' This function is used to programatically take vectors of paramters and write suites of R parameter scripts
#' that will form part of a SHAPE experiment that is being built for running.  This is a wrapper for writting out
#' the suite of necessary scripts to form a run.
#'
#' @param func_subScipt This is the template script that needs to be replicated
#' @param func_outDir This is the filepath directory where output should be placed
#' @param func_inCombos This is the combinations of parameters that are to be used in the experiment.
#' @param func_inParms # These are additional parameters to be implemented in writing out combinations.
#' @param func_maxJobs This is the maximum number of individual R jobs that should be called at once
#' by the shell submission scripts, it can affect both local and remote server calls.
#' @param func_appLocation This is the filepath for R so that batch mode runs can be called.
#' @param func_commonArgs These are common arguments important when running the batch mode
#' @param func_submitArgs These are common arguments important when submitting the batch mode
#' @param func_remoteLocation This is a remote server location where an experiment built is to be run
#' it affects the filepathing called by submission scripts and the associated batch mode runs performed.
#' @param func_passedArgs These are arguments passed through this wrapper to inner functions.
#' @param func_externalStopper This is a file which exists as a flag for stopping SHAPE from trying to create
#' additional replicates.
#' @param func_sepString This is the common string used to collapse information.
#'
#' @return A character string that should indicate the experiment has been created.  Otheriwse this has failed.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
write_subScript <- function(func_subScipt, func_outDir, func_inCombos, func_inParms, func_maxJobs,
                            func_appLocation, func_commonArgs, func_submitArgs,
                            func_remoteLocation, func_passedArgs,
                            func_externalStopper = getOption("shape_external_stopFile"),
                            func_sepString = getOption("shape_sepString")){
  func_submitLines <- readLines(func_subScipt,warn=FALSE)
  # If we're doing external selfing we add in some lines to the start and end of the script
  if(as.logical(func_inParms$shape_externalSelfing)){
    # I find the bash line and add a first line after it but before the rest
    func_tmpLine <- which(grepl("#!/bin/bash",func_submitLines,fixed=TRUE))
    func_submitLines <- c(func_submitLines[1:func_tmpLine],
                          paste('echo "" > fake_workDir/',func_externalStopper,sep=""),
                          func_submitLines[(func_tmpLine+1):length(func_submitLines)],
                          paste('if ! [ -e "fake_workDir/',func_externalStopper,'" ]',sep=""),
                          'then',
                          '	fake_workDir/fakeSubmit',
                          'fi')
  }
  # We build the general outline of the server submission script file, but not all elements are included if we're not on the server
  func_submitLines <- gsub("fake_appPath",func_appLocation,func_submitLines)
  func_submitLines <- gsub("fake_commandArgs",func_commonArgs,func_submitLines)
  func_submitLines <- gsub("fake_passedArgs",func_passedArgs,func_submitLines)
  if(as.logical(func_inParms$shape_serverFarm)){
    # We now add in the submission commands to this templated call
    func_tmpLine <- which(grepl("fake_subMit_command",func_submitLines))
    func_submitLines <- c(func_submitLines[1:(func_tmpLine-1)],
                          paste("#SBATCH ",func_submitArgs,sep=""),
                          func_submitLines[(func_tmpLine+1):length(func_submitLines)])
    func_submitLines <- gsub("fake_serverPath",func_remoteLocation,func_submitLines)
  }

  # now since we're building a submission sccript for each of the job, we loop through them, but tracking how many
  # are placed into any one file so that each submission script has no more than maxJobs worth
  func_tmpCounter <- 0
  func_jobBatch <- 1
  func_masterCon <- file(description = paste(func_outDir,name_batchSubmit(func_jobBatch),sep=""), open = "w")
  for(thisCombo in 1:nrow(func_inCombos)){
    for(shape_thisRep in 1:as.numeric(func_inParms$uniqueReplicates)){
      tmp_jobString <- name_batchString(funcBase = trimQuotes(func_inParms$shape_save_batchBase),
                                        func_setID = shape_thisRep,
                                        func_jobID = func_inCombos[thisCombo,"shape_save_batchJob"],
                                        func_sepString = func_sepString)
      # I define this job's working directory
      tmp_jobDir <- paste(func_outDir,
                          tmp_jobString,
                          "/",sep="")
      tmpSubmit_fileName <- name_subScript(tmp_jobString)
      # I now also build the submission script file for this job, some elements are for the server only
      tmp_submit <- gsub("fake_workDir",tmp_jobDir,func_submitLines)
      tmp_submit <- gsub("fake_tmpScript.r",name_bodyScript(tmp_jobString),tmp_submit)

      if(as.logical(func_inParms$shape_serverFarm)){
        tmp_submit <- gsub("fakeOut",paste(tmp_jobString,".o",sep=""),tmp_submit)
        tmp_submit <- gsub("fakeJob",tmp_jobString,tmp_submit)
        tmp_submit <- gsub("fakeDir",tmp_jobString,tmp_submit)
      } else {
        tmp_submit <- gsub("fake_serverPath/fakeDir/",tmp_jobDir,tmp_submit)
      }
      # If we're selfing externally then I update the submission script name
      if(as.logical(func_inParms$shape_externalSelfing)){
        tmp_submit <- gsub("fakeSubmit",tmpSubmit_fileName,tmp_submit)
      }

      # I now write out the job's submission script
      writeLines(tmp_submit, con = paste(tmp_jobDir,tmpSubmit_fileName,sep=""))
      # Now I write to the jobBatch submission script
      func_tmpCounter <- func_tmpCounter + 1
      # This updates our counters if need be and differ based on local vs server farming.
      if(as.logical(func_inParms$shape_serverFarm)){
        if(func_tmpCounter > func_maxJobs){
          func_tmpCounter <- func_tmpCounter - (func_maxJobs + 1)
          func_jobBatch <- func_jobBatch + 1
          # We close the old batch script and open a new one
          close(func_masterCon)
          func_masterCon <- file(description = paste(func_outDir,name_batchSubmit(func_jobBatch),sep=""), open = "w")
        }
      } else {
        if(func_tmpCounter > ceiling((nrow(func_inCombos) * as.numeric(func_inParms$uniqueReplicates))/func_maxJobs)){
          func_tmpCounter <- 0
          func_jobBatch <- func_jobBatch + 1
          close(func_masterCon)
          func_masterCon <- file(description = paste(func_outDir,name_batchSubmit(func_jobBatch),sep=""), open = "w")
        }
      }
      tmp_masterSubmit <- paste(if(as.logical(func_inParms$shape_serverFarm)){"sbatch "}else{NULL},
                                tmp_jobDir,tmpSubmit_fileName,sep="")
      writeLines(tmp_masterSubmit, con = func_masterCon, sep="\n")
    }
  }
  # We close the connection and report completing the tasks.
  close(func_masterCon)
  return( "Done writting out submission scripts" )
}



# This will write out a parameter file for each job, these will be inherited by the script as it runs
#' This is a file for updating the post analysis plotting script and creating an updated copy in the experiment's folder
#'
#' @param func_infile This is the filepath location for the template script to be writte in.
#' @param func_inParms These are the parameters to be updated in the plotting file
#' @param func_inCombos This is the combination of parameters to be written
#' @param func_outDir This is the director filepath to which output should be written.
#' @param func_bodyScript This is a run body of SHAPE script to be read in as template
#' @param func_ExternalStopper This is a file placed externally used as a logical flag that SHAPE should
#' stop trying to seed new replicates to be run.
#' @param func_sepString This is the common string for collapsing information.
#'
#' @return A character string indicating that the plotting file-s- have been written
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
writeParameters <- function(func_infile, func_inParms, func_inCombos, func_outDir, func_bodyScript,
                            func_ExternalStopper, func_sepString = getOption("shape_sepString")){
  # Step 1 we read in the lines of the script we're looking to update
  func_parmLines <- readLines(func_infile, warn=FALSE)
  func_bodyLines <- readLines(func_bodyScript, warn=FALSE)
  # So, for each set of parameters, we go and update lines.  Starting with the parameter combinations that don't
  # change we set out template items to follow those values
  func_constParms <- c("shape_serverFarm","shape_results_removeSteps","shape_externalSelfing","shape_toggle_forceCompletion",
                       "shape_save_batchBase", "shape_maxReplicates","shape_workDir")
  func_parmLines <- updateLines(func_inLines = func_parmLines,
                                func_searchPattern = sapply(names(func_inParms[func_constParms]),paste," <-","sep"=""),
                                func_values = func_inParms[func_constParms])
  # Right, now we need to update and write out a parameter file and body file for each parameter combination.
  for(thisCombo in 1:nrow(func_inCombos)){
    tmp_parmLines <- updateLines(func_inLines = func_parmLines,
                                 func_searchPattern = sapply(colnames(func_inCombos),paste," <-","sep"=""),
                                 func_values = func_inCombos[thisCombo,colnames(func_inCombos)])
    # Now we cycle through each independent statistical replicate and write out files
    for(shape_thisRep in 1:as.numeric(func_inParms$uniqueReplicates)){
      tmp_jobString <- name_batchString(funcBase = trimQuotes(func_inParms$shape_save_batchBase),
                                        func_setID = shape_thisRep,
                                        func_jobID = func_inCombos[thisCombo,"shape_save_batchJob"],
                                        func_sepString = func_sepString)
      # I define this job's working directory
      tmp_jobDir <- paste(func_outDir,
                          tmp_jobString,
                          "/",sep="")
      tmpParm_fileName <- name_parameterScript(tmp_jobString)
      tmpBody_fileName <- name_bodyScript(tmp_jobString)
      # I now update the body and parm lines with respect to the working directory and jobSet
      tmp_bodyLines <- updateLines(func_inLines =func_bodyLines,
                                   func_searchPattern = list("source"="sourceParms <-"),
                                   func_values = list("source"= addQuotes(paste(tmp_jobDir,tmpParm_fileName,sep=""))))

      tmp_jobParms <- updateLines(func_inLines =tmp_parmLines,
                                  func_searchPattern = list("shape_save_batchSet"="shape_save_batchSet <-",
                                                            "shape_save_batchJob"="shape_save_batchJob <-"),
                                  #func_values = list("shape_workDir"= addQuotes(paste(sub(func_outDir,trimQuotes(func_inParms$shape_workDir),tmp_jobDir),sep="")),
                                  func_values = list("shape_save_batchSet"= shape_thisRep,
                                                     "shape_save_batchJob"=func_inCombos[thisCombo,"shape_save_batchJob"]))
      # Now I write out both files into their respective job directory
      writeLines(tmp_bodyLines, con = paste(tmp_jobDir,tmpBody_fileName,sep=""))
      writeLines(tmp_jobParms, con = paste(tmp_jobDir,tmpParm_fileName,sep=""))
    }
  }
  return( "Have written out all job's parameters and body scripts" )
}


#' This is a function which is used to update lines that are searched and replace in a manner conditional to this script's circumstances
#' The input lines can be a vector of any length, and the search patterns can be a list of any length where each list vector is used together.
#' The values should be a list of information used as replacement info.
#'
#' @param func_inLines These are the lines that are to be updated before output
#' @param func_searchPattern These are the string-s- to be searched for replacement
#' @param func_values These are the values that are to replace the searched strings.
#'
#' @return A vector of character strings that has now been updated.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
updateLines <- function(func_inLines,func_searchPattern, func_values){
  # There is an asssumption that the searchPattern and values are equally named
  if(!all(is.element(names(func_searchPattern),names(func_values)))){
    stopError("There were not similar searchPattern and values passed to updateLines, please review")
  }
  # For each of the names search patterns we'll try to update the first instance of the line in the code
  func_updateLines <- NULL
  for(thisUpdate in names(func_searchPattern)){
    func_updateLines <- func_inLines
    # The way lines are update depends on the length of the searchPattern information
    func_tmpLine <- which(grepl(func_searchPattern[[thisUpdate]][1],func_updateLines))[1]
    if(length(func_searchPattern[[thisUpdate]]) == 1){
      # This means we simply search for the line and replace it as the pattern and the value pasted
      func_updateLines[func_tmpLine] <- paste(func_searchPattern[[thisUpdate]],func_values[[thisUpdate]][1],sep=" ")
    } else if(length(func_searchPattern[[thisUpdate]]) == 2){
      # If there are two search patterns, the second is used to inform what is replaced in the line
      func_updateLines[func_tmpLine] <- sub(func_searchPattern[[thisUpdate]][2],func_values[[thisUpdate]][1],func_updateLines[func_tmpLine])
    } else if(length(func_searchPattern[[thisUpdate]]) == 3){
      # This means we'll be using the first element to find the line, and the other two to know how to position our change.
      # We use the first and second to identify positions in the string, and then use the region earliest in the string
      # but after the first search term.  I build these as matrices for convenience
      func_tmpPositions <- list(regexpr(func_searchPattern[[thisUpdate]][2],func_updateLines[func_tmpLine],fixed=TRUE),
                                gregexpr(func_searchPattern[[thisUpdate]][3],func_updateLines[func_tmpLine]))
      func_tmpPositions <- lapply(func_tmpPositions,function(x){
        matrix(c(unlist(x),attr(if(is.list(x)){x[[1]]}else{x},"match.length")),nrow=2,
               byrow = TRUE,dimnames=list(c("location","length"),NULL))
      })
      # We check that the two positions were found
      if(!all(sapply(func_tmpPositions,function(x){ any(x["location",] > 0) }))){
        stopError(paste("Could not find the proper replacement positions for ",thisUpdate," please review",sep=""))
      }
      # The positions to be used will be the first of the second searches which is larger than the first
      func_usePosition <- which(func_tmpPositions[[2]]["location",] > func_tmpPositions[[1]]["location",1])[1]
      func_tmpPositions <- func_tmpPositions[[2]][,func_usePosition]

      # This means the line will become what was previously there, cut around the positions, with the values inserted
      func_updateLines[func_tmpLine] <- paste(substr(func_updateLines[func_tmpLine],1,func_tmpPositions[1]-1),
                                              func_values[[thisUpdate]],
                                              substr(func_updateLines[func_tmpLine],sum(func_tmpPositions),nchar(func_updateLines[func_tmpLine])),
                                              sep="")
    }
    func_inLines <- func_updateLines
  }
  return( func_inLines )
}


#' This is a function to take the input parameters and build the parameter combinations
#'
#' @param func_inLines These are the template lines of text to be updated.
#' @param func_comboRef This is the reference identifiers for grouped as pairwise parameter combinations
#' @param func_indepRef This is the reference identifiers for independent parameter values not to be done pairwise
#' @param func_condRef This is the reference indetifiers for grouped parameter combinations which are conditional on others.
#'
#' @return A table of parameter combinations which represents the combination of experimental parameters for a SHAPE experiment.
#'
#' @section Note:
#' There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
#'
#' @export
shapeCombinations <- function(func_inLines, func_comboRef, func_indepRef, func_condRef){
  # Step one, we'll build all possible combinations of parameters that are going to be expanded
  func_parmGrid <- expand.grid(func_inLines[setdiff(names(func_inLines),c(func_comboRef,
                                                                          func_indepRef,
                                                                          unname(unlist(func_condRef)),
                                                                          "shape_save_batchJob"))],
                               stringsAsFactors=FALSE,
                               KEEP.OUT.ATTRS = FALSE)
  # Now we build all the combinatorial references but start by ensuring the vectors are of the same length, they should be,
  # but if the user has not done this we'll recycle the smaller vectors until they are of length for the longest one
  func_tmpLengths <- sapply(func_inLines[func_comboRef],length)
  if(!all(max(func_tmpLengths) == func_tmpLengths)){
    func_tmpUpdate <- which(func_tmpLengths != max(func_tmpLengths))
    for(thisUpdate in func_tmpUpdate){
      func_inLines[func_comboRef[thisUpdate]] <- list(unname(unlist(rep(func_inLines[func_comboRef[thisUpdate]],max(func_tmpLengths))))[1:max(func_tmpLengths)])
    }
  }
  # We now build on the conditional parameter combinations into a matrix with the conditional combinations
  func_comboAdd <- matrix(unlist(func_inLines[func_comboRef]),nrow=max(func_tmpLengths),dimnames=list(NULL,func_comboRef))
  for(thisCond in names(func_condRef)){
    # We record the columnames prior to updating this is ued in the later steps of this loops
    func_startCols <- colnames(func_comboAdd)
    # We find which rows in the combinations hold the simulation model related to our conditional
    func_tmpLines <- which(grepl(thisCond,func_comboAdd[,"shape_simModel"]))
    # We build the combinations of parameters for this conditional
    func_condAdd <- as.matrix(expand.grid(func_inLines[func_condRef[[thisCond]]],
                                          stringsAsFactors=FALSE, KEEP.OUT.ATTRS = FALSE))
    # we now add the first condAdd elements to the func_comboAdd matrix
    for(thisAddition in 1:ncol(func_condAdd)){
      # We add the first conditional combination, then any additional required
      func_comboAdd <- cbind(func_comboAdd,func_condAdd[1,thisAddition])
    }
    colnames(func_comboAdd)[(ncol(func_comboAdd)-(ncol(func_condAdd)-1)):ncol(func_comboAdd)] <- colnames(func_condAdd)
    if(length(func_tmpLines) > 0){
      if(nrow(func_condAdd) > 1){
        for(thisRow in 2:nrow(func_condAdd)){
          for(thisUnique in func_tmpLines){
            # We need to add a number of rows equal to the number of func_tmpLines
            func_comboAdd <- rbind(func_comboAdd, c(unlist(func_comboAdd[thisUnique,func_startCols]),unlist(func_condAdd[thisRow,])))
          }
        }
      }
    }
  }
  # Now for each row in func_comboAdd, we need to have nrow(func_parmGrid) replicates of that row
  func_tmpNew1 <- NULL
  for(thisRow in 1:nrow(func_comboAdd)){
    func_tmpNew1 <- rbind(func_tmpNew1,matrix(rep(func_comboAdd[thisRow,],nrow(func_parmGrid)),
                                              nrow=nrow(func_parmGrid), byrow = TRUE,
                                              dimnames=list(NULL,colnames(func_comboAdd))))
  }
  # We now replicate the func_parmGrid a number of times equal to the rows of the other combinations
  func_tmpNew2 <- NULL
  for(thisAddition in 1:nrow(func_comboAdd)){
    func_tmpNew2 <- rbind(func_tmpNew2,func_parmGrid)
  }
  # We now stick the pieces together
  func_parmGrid <- cbind(func_tmpNew2,func_tmpNew1,stringsAsFactors=FALSE)
  # We only keep unique combinations
  func_parmGrid <- unique(func_parmGrid)
  # I now build the shape_save_batchJob vector, which is the "job" ID for these parameters
  func_parmGrid <- cbind(func_parmGrid,as.numeric(func_inLines$shape_save_batchJob):
                           (nrow(func_parmGrid)+(as.numeric(func_inLines$shape_save_batchJob)-1)),
                         stringsAsFactors=FALSE)
  colnames(func_parmGrid)[ncol(func_parmGrid)] <- "shape_save_batchJob"

  # I now return this
  return( func_parmGrid )
}







#' This function is a wrapper for getting a summary of the results of an rSHAPE run and/or experiment as a whole.  The
#' former is presumed to be of greater use but either is fine as per your needs.  This wrapper will cause
#' RData files to be created which contain the summarised experimental details that you can then use more easily
#' for analysis.
#'
#' @param func_processingTypes A vector of character strings which define the type of processing to be performed
#' when callign this experimental analysis wrapper function.  At present, the types include:
#' "fileList", "parameters", "popDemographics","repeatability" as per the rSHAPE option - shape_procExp_filenames
#' @param func_numCores Integer number of computer cores to be requested for performing parallel processing
#' of experiment files.  It defaults as 1, which effectively means in tandem - ie: not parallel.
#' @param func_suppressOld This is a logical toggle if files which exist in the expected location should be deleted.
#' Default is FALSE and the function will simply not process alraedy processed output.  TRUE might be useful as a means
#' to forcibly re-run the summary fresh.
#'
#' @return A message detailing if the requested processed files can be found, either affirmative for all
#' or a note when at least one is missing.
#'
#' @section Note:
#' There is no example as this cannot work without a complete rSHAPE experiment to be analysed.
#'
#' @export
summariseExperiment <- function(func_processingTypes = c("fileList", "parameters", "popDemographics","repeatability"),
                                func_numCores = 1, func_suppressOld = FALSE){
  # For convenience we grab the values of the processing filenames
  func_processingOptions <- getOption("shape_procExp_filenames")
  # If the user is having a laugh and passing nothing for the processing we'll drop out
  func_processingTypes <- func_processingTypes[which(is.element(func_processingTypes,
                                                                names(func_processingOptions)))]
  if(length(func_processingTypes) == 0 ){
    return("No valid func_processingTypes strings passed, not trying to run experiment processing.")
  }
  # Now we check if we're supposed to delete any files
  if(func_suppressOld){
    # We check which files exist to be removed
    func_tmpRemove <- func_processingOptions[which(file.exists(func_processingOptions))]
    if(length(func_tmpRemove) > 0){ file.remove(func_tmpRemove) }
  }

  # This function is the one which can take advantage of parallelisation.  So here is where we'll
  # register multiple cores or make a cluster for processing.
  func_useCluster <- parallel::makeCluster(func_numCores, type="PSOCK")
  doParallel::registerDoParallel(func_useCluster)


  # This is a regular expression like string that can be used to find files.
  funcSave_jobExpression <- name_batchString(funcBase = getOption("shape_save_batchBase"),
                                             func_setID = if(!is.null(getOption("shape_save_batchSet"))){
                                                             paste('[^',getOption("shape_sepString"),']+',sep="")
                                                           }else{
                                                             NULL
                                                           },
                                             func_jobID = if(!is.null(getOption("shape_save_batchJob"))){
                                                            paste('[^',getOption("shape_sepString"),']+',sep="")
                                                          }else{
                                                            NULL
                                                          },
                                             func_sepString = getOption("shape_sepString"))

  # Now we call the different functions for summary post-analysis
  if(is.element("fileList",names(func_processingOptions)) &&
     !file.exists(func_processingOptions["fileList"])){
        # We call the function which has standard shape options as all arguments.
        summarise_experimentFiles()
  }
  if(is.element("parameters",names(func_processingOptions)) &&
            !file.exists(func_processingOptions["parameters"])){
        # We call the function which has standard shape options as all arguments.
        summarise_experimentParameters()
  }
  if(is.element("popDemographics",names(func_processingOptions)) &&
            !file.exists(func_processingOptions["popDemographics"])){
        # We call the function which takes only a single non-default argument based on reference strings
        summarise_popDemographics(funcSave_jobExpression = funcSave_jobExpression)
  }
  if(is.element("repeatability",names(func_processingOptions)) &&
            !file.exists(func_processingOptions["repeatability"])){
        # We call the function which takes only a single non-default argument based on reference strings
        summarise_evolRepeatability(funcSave_jobExpression = funcSave_jobExpression)
  }

  # We close the cluster that was built
  stopCluster(func_useCluster)

  # We check which of the output summaries exist
  func_existingOutput <- file.exists(func_processingOptions[func_processingTypes])
  # We return a message based on our search for output files.
  return( ifelse(all(func_existingOutput),
                     "Reached end of summariseExperiment, output summaries found.",
                     "Missing at least one of the requested summary outputs from summariseExperiment.") )
}



#' This function will find all initially processed output files from individual replicates and return summary information.
#' That information is saved to an RData file which will contain 3 objects: all_proccessedFiles, all_jobInfo, all_dividedFiles
#'
#' @param func_experimentDir This is the filepath to the root directoy under which all your experimental files can
#' be found.
#' @param func_saveFile This is the filepath and filename (ending in .RData please) to which the results of this
#' step will be saved.
#' @param func_search_filePattern This is a string which can be used to search and find the files which relate to
#' the processed output of individual replicates rSHAPE runs.
#' @param func_sepString This is the character string which was used for commonly collapsing elements in the rSHAPE run.
#'
#' @section Note:
#' There is no example as this cannot work without a complete rSHAPE experiment to be analysed.
#'
#' @export
summarise_experimentFiles <- function(func_experimentDir = getOption("shape_workDir"),
                                      func_saveFile = getOption("shape_procExp_filenames")["fileList"],
                                      func_search_filePattern = getOption("shape_processedData_filePattern"),
                                      func_sepString = getOption("shape_sepString")){
  # This will go and find every single processed output file that can be used for the experimental summary.
  all_proccessedFiles <- list.files(path= func_experimentDir,
                                    pattern= func_search_filePattern,
                                    recursive = TRUE, full.names=TRUE)
  # Now I know what the epxression is for jobs, but better yet I know how the jobs were built and so can extract all my info
  all_jobInfo <- name_batchString(funcBase = paste(getOption("shape_save_batchBase"),
                                                   unlist(lapply(strsplit(all_proccessedFiles, func_search_filePattern),function(thisFile){
                                                     gsub('.RData',"",thisFile[2],fixed=TRUE)
                                                   })),
                                                   sep=""),
                                  func_setID = !is.null(getOption("shape_save_batchSet")),
                                  func_jobID = !is.null(getOption("shape_save_batchJob")),
                                  func_repID = TRUE,
                                  funcSplit = TRUE,
                                  func_sepString = func_sepString)
  # We now find unique jobs by looking at the information just after the func_search_filePattern
  # The regular expressions are so that I can have subsetted replicates of jobs yet they be considered part of a single set.
  all_uniqueJobs <- name_batchString(funcBase = getOption("shape_save_batchBase"),
                                     func_setID = paste('[^',func_sepString,']+',sep=""),
                                     func_jobID = unique(all_jobInfo["jobID",]),
                                     func_sepString = func_sepString)
  # We now subdivide our all_proccessedFiles into a list of the different unique jobs
  # If you're not aware sets are used by rSHAPE to create replicates of identical parameter
  # combinations (which are jobs - and each job can have many replicates).
  all_dividedFiles <- sapply(all_uniqueJobs,function(x){
                          return( all_proccessedFiles[which(grepl(paste(x,"(.+)",sep=func_sepString), all_proccessedFiles))] )
                        },simplify=FALSE)
  # As a sanity check we ensure that all the all_proccessedFiles have been placed in our dividedfiles object
  if(length(all_proccessedFiles) != length(unlist(all_dividedFiles))){
    stopError(" There was a problem when dividing the all_proccessedFiles into divdedFiles, please review")
  }
  save(all_proccessedFiles, all_jobInfo, all_dividedFiles,
       file= func_saveFile)
  # Nothing is returned, we'll have saved an output
  invisible( NULL )
}



#' This function will use output from summarise_experimentFiles to locate all parameter files and then
#' report on all those parameters for the jobs that were run.  This will save an RData file which will
#' contain one object: all_parmInfo
#'
#' @param func_workEnvir This is an environment used to load files with the load function.  It's used to encapsulate
#' the loaded information to a controlled space.
#' @param func_saveFile This is the filepath and filename (ending in .RData please) to which the results of this
#' step will be saved.
#' @param func_experimentDir This is the filepath to the root directoy under which all your experimental files can
#' be found.
#' @param func_refFile This is the filepath to the reference file that contains information regarding all the
#' processed files for the rSHAPE experiment.
#'
#' @section Note:
#' There is no example as this cannot work without a complete rSHAPE experiment to be analysed.
#'
#' @export
summarise_experimentParameters <- function(func_workEnvir = new.env(),
                                           func_saveFile = getOption("shape_procExp_filenames")["parameters"],
                                            func_experimentDir = getOption("shape_workDir"),
                                            func_refFile = getOption("shape_procExp_filenames")["fileList"]){
  # We start by loading the reference file into this workspace, but if it does not exist we complain
  if(!file.exists(func_refFile)){
    stopError("Did not find the reference file needed for summarise_experimentParameters")
  }
  load(func_refFile, envir = func_workEnvir)

  # We find all the directories which use the names of our divided files
  func_tmpDirs <- list.dirs(path=func_experimentDir,full.names=FALSE,recursive=FALSE)
  func_tmpDirs <- unique(unlist(lapply(names(func_workEnvir$all_dividedFiles), function(thisName){
                                    return( func_tmpDirs[which(grepl(thisName,func_tmpDirs))] )
                                  })))

  # I need to create an outer list object which stores the information for each job
  all_parmInfo <- NULL
  # Now we break the processing into chunks to be added up and collected
  # This done so as to reduce the size of the write packets.
  tmp_chunkSize <- 100
  tmp_numChunks <- ceiling(length(func_tmpDirs)/tmp_chunkSize)
  for(thisChunk in 1:tmp_numChunks){
    tmpAdd <- NULL
    for(thisJob in func_tmpDirs[(1+((thisChunk - 1)*tmp_chunkSize)):min((thisChunk*tmp_chunkSize),length(func_tmpDirs))]){
      # We load the parameters for thisJob
      tmp_parmFile <- list.files(path=paste(func_experimentDir,thisJob,"/",sep=""), pattern = "Parameters")
      if(length(tmp_parmFile) == 1){
        load(paste(func_experimentDir,thisJob,"/",tmp_parmFile,sep=""),envir= func_workEnvir)
      } else {
        stopError(paste("Could not find the parameters files to load for ",thisJob,sep=""))
      }
      # I now build an object which captures the runParms of interest, this should be a list object
      tmp_parmsRef <- c(func_workEnvir$runParameters$Population,
                        func_workEnvir$runParameters$Growth_Disturbance[which(!is.element(names(func_workEnvir$runParameters$Growth_Disturbance),
                                                                                          c("shape_const_growthGenerations",
                                                                                            "shape_init_distPars",
                                                                                            "shape_track_distSize")))],
                        "distFactor"= unname(func_workEnvir$runParameters$Growth_Disturbance$shape_init_distPars["factor"]),
                        "distSpread"= unname(func_workEnvir$runParameters$Growth_Disturbance$shape_init_distPars["random"]),
                        "mean_realisedDilution"=mean(func_workEnvir$runParameters$Growth_Disturbance$shape_track_distSize[,"factor"]),
                        func_workEnvir$runParameters$FitnessLandscape,
                        func_workEnvir$runParameters$DFE,
                        "db_splitTables"= func_workEnvir$runParameters$DataManagement$shape_db_splitTables)
      tmp_parmsRef$shape_const_distParameters <- paste(tmp_parmsRef$shape_const_distParameters,collapse=",")

      # Now I remove the runParameters from the workspace as a sanity check
      rm(list=c("runParameters"), envir = func_workEnvir)
      # We now assign this information into our list object, the jobs set looks for which of the names(func_workEnvir$all_dividedFiles)
      # associates with thisJob for latter grouping of data.  To assign a unique jobSet we look if there is the jobString
      # otherwise use the directory name.
      tmpReturn <- data.frame("jobName"=thisJob,
                              "jobSet"= if(!is.null(getOption("shape_save_batchSet"))){
                                            # We then check which of the func_workEnvir$all_dividedFiles names could refer to this
                                            # We then choose the longest one that matches, as that will be the most complete
                                            tmpNames <- names(func_workEnvir$all_dividedFiles)[sapply(names(func_workEnvir$all_dividedFiles),grepl,"x"=thisJob)]
                                            tmpNames[which.max(nchar(tmpNames))]
                                          } else {
                                            thisJob
                                          },
                              t(unlist(tmp_parmsRef)),
                              stringsAsFactors = FALSE)
      # We now just reset the data types
      for(thisCol in names(tmp_parmsRef)){ mode(tmpReturn[,thisCol]) <- mode(tmp_parmsRef[[thisCol]]) }
      # Build up the chunk of info.
      tmpAdd <- rbind(tmpAdd, tmpReturn)
    }
    # Build on the chunk of info
    all_parmInfo <- rbind(all_parmInfo, tmpAdd)
  }
  # We save the gathered information
  save(all_parmInfo, file= func_saveFile)

  # We now silently return nothing
  invisible( NULL )
}





#' This function will use output from summarise_experimentFiles and summarise_experimentParameters
#' to help with expectations concerning run output and handling.  This will save an RData file which
#' will contain one object: all_popSets, which is a list of relevant control information about I/O
#' and then a series of other RData files which contain the demographics information as a matrix with the
#' mean and standard deviation of demographics for all replicates.
#'
#' @param funcSave_jobExpression This is a string expression that can be used to find elements of the experiment
#' being analysed.  It should be some robust unique string or regular expression.
#' @param func_saveFile This is the filepath and filename (ending in .RData please) to which the results of this
#' step will be saved.
#' @param func_experimentDir This is the filepath to the root directoy under which all your experimental files can
#' be found.
#' @param func_saveDir This is the directory to which output will be saved.
#' @param func_refFile This is the filepath to the reference file that contains information regarding all the
#' processed files for the rSHAPE experiment.
#' @param func_workEnvir This is an environment used to load files with the load function.  It's used to encapsulate
#' the loaded information to a controlled space.
#' @param func_objPrefix This is a character string for programatic naming of objects of this type.
#'
#' @section Note:
#' There is no example as this cannot work without a complete rSHAPE experiment to be analysed.
#'
#' @export
summarise_popDemographics <- function(funcSave_jobExpression,
                                      func_saveFile = getOption("shape_procExp_filenames")["popDemographics"],
                                      func_experimentDir = getOption("shape_workDir"),
                                      func_saveDir = getOption("shape_postDir"),
                                      func_refFile = getOption("shape_procExp_filenames")[c("fileList","parameters")],
                                      func_workEnvir = new.env(),
                                      func_objPrefix = "popDemo_"){
  # This step can take time with larger experiments so we print a flag
  print("Starting popDemographics")
  # We start by loading the reference file into this workspace, but if it does not exist we complain
  if(any(!file.exists(func_refFile))){
    stopError("Did not find the reference file needed for summarise_experimentParameters")
  }
  for(thisFile in func_refFile){
    load(thisFile, envir = func_workEnvir)
  }
  # We create an object that tracks the names of save files, objects and jobs for each set
  all_popSets <- sapply(unique(func_workEnvir$all_parmInfo$jobSet), function(thisSet){
                      tmp_thisName <- nameObject(func_inString = thisSet,
                                                 func_inPrefix = func_objPrefix)
                      return( list("objName"=tmp_thisName,
                                   "saveFile"= paste(sub(funcSave_jobExpression,"",
                                                          sub(getOption("shape_save_batchBase"),
                                                              tmp_thisName,
                                                              func_saveFile),
                                                         fixed=TRUE),sep=""),
                                   "setJobs" = unique(func_workEnvir$all_parmInfo$jobName[which(func_workEnvir$all_parmInfo$jobSet == thisSet)])) )

                    }, simplify = FALSE)
  # For all_uniqueSets, we find the ones where the files don't exist so those are the only we build
  tmp_missingSets <- names(all_popSets)[unlist(lapply(all_popSets,function(thisSet){
                                                  !file.exists(paste(func_saveDir,thisSet[["saveFile"]],sep=""))
                                                }))]
  # This assignment is to make check() shut up despite no issue with global binding within the downstream foreach context.
  thisCheck <- NULL

  # Now for each main job we'll gather some information and save it into it's own list object
  tmpCheck <- foreach(thisCheck = tmp_missingSets, .combine="c",
                      .packages = c("rSHAPE", "DBI","RSQLite")) %dopar% {

    tmpList <- vector(mode="list",
                      length=length(all_popSets[[thisCheck]][["setJobs"]]))
    names(tmpList) <- all_popSets[[thisCheck]][["setJobs"]]
    # Ok here is the level at which we'll place our foreach looping to create a list of results
    for(thisJob in all_popSets[[thisCheck]][["setJobs"]]){
       # We find the job files associated to thisJob
       tmpJobs <- func_workEnvir$all_dividedFiles[[thisCheck]][which(grepl(thisJob, func_workEnvir$all_dividedFiles[[thisCheck]]))]
       if(length(tmpJobs) > 0){
         # I need to create an outer list object which stores the information for each job
         tmp_replicateInfo <- NULL
         # Now we break the processing into chunks to be added up and collected
         tmp_chunkSize <- 100
         tmp_numChunks <- ceiling(length(tmpJobs)/tmp_chunkSize)
         for(thisChunk in 1:tmp_numChunks){
           tmpAdd <- NULL
           # Now for each and every replicate within this job we want to grab the:
           # (1) the population wide - through time - min, mean,max fitness. - I'll return the whole demoMat
           for(thisFile in tmpJobs[(1+((thisChunk - 1)*tmp_chunkSize)):min((thisChunk*tmp_chunkSize),length(tmpJobs))]){
             load(thisFile, envir = func_workEnvir)
             # Ok we'll now return this information along with the demoMat
             tmpAdd <- c(tmpAdd, list("popDemographics"= func_workEnvir$runDemographics$demoMat))
             rm(list=c("info_estLines","runDemographics" ), envir = func_workEnvir)
             # Clear up memory
             gc()
           } # This closes out the for loop gathering information about the replicates for this job
           tmp_replicateInfo <- c(tmp_replicateInfo, tmpAdd)
           rm(tmpAdd)
           gc()
         }

         #tmp_all_popDemos <- array(sapply(tmp_replicateInfo,function(x){ x[["popDemographics"]] }),
         tmp_all_popDemos <- array(sapply(tmp_replicateInfo,function(x){ x }),
                                   dim=c(dim(tmp_replicateInfo[[1]]),length(tmp_replicateInfo)),
                                   dimnames = list(rownames(tmp_replicateInfo[[1]]),
                                                   colnames(tmp_replicateInfo[[1]]),
                                                   NULL))
         # For the population demograpihcs we collapse this into the mean values through time across all replicates
         tmp_all_popDemos <- matrix(apply(tmp_all_popDemos,MARGIN=2,function(thisCol){
                                         # Now we return the mean and sd values for each column
                                         c(apply(thisCol,MARGIN=1,mean),apply(thisCol,MARGIN=1,sd))

                                       }),nrow=nrow(tmp_all_popDemos),ncol=ncol(tmp_all_popDemos)*2,
                                       dimnames=list(rownames(tmp_all_popDemos),
                                                     unlist(lapply(colnames(tmp_all_popDemos),function(x){ return(c(x,paste("sd",x,sep="_")))}))) )
         # We now assign this information into our list object
         tmpList[[thisJob]] <- tmp_all_popDemos
         rm(tmp_all_popDemos, tmp_replicateInfo)
         gc()
       } # This closes out the logical conditional that there is at least some replicates in thisDir
       # This closes out the loop for the different jobs
    }
    # We can now save this object to a file
    assign(all_popSets[[thisCheck]][["objName"]], tmpList, pos= func_workEnvir)
    save(list = all_popSets[[thisCheck]][["objName"]],
         file= all_popSets[[thisCheck]][["saveFile"]],
         envir = func_workEnvir)
    rm(list= all_popSets[[thisCheck]][["objName"]], envir = func_workEnvir)
    rm(tmpList)
    # Clear up memory
    gc()
    # We return a confirmation that the job completed
    return(  nameObject(func_inString = all_popSets[[thisCheck]][["objName"]],
                        func_inPrefix = func_objPrefix,
                        func_splitStr = TRUE)  )
  }
  # We now check if all the jobs completed
  tmp_checkCompleted <- sapply(tmp_missingSets,is.element,"set"=tmpCheck)
  if(all(tmp_checkCompleted)){
    # We now save the object which stores all the population demographics data object names and files locations
    save(all_popSets, file= func_saveFile)
  } else {
    stop(paste("There was a problem with the sets of: ",paste(names(tmp_checkCompleted)[which(!tmp_checkCompleted)],collapse=" ",sep=" "),sep=""))
  }

  # We silently return nothinig
  invisible( NULL )
}



#' This function will use output from summarise_experimentFiles and summarise_experimentParameters
#' to help with expectations concerning run output and handling.  This will save an RData file which
#' will contain one object: all_popSets, which is a list of relevant control information about I/O
#' and then a series of other RData files which contain the demographics information as a matrix with the
#' mean and standard deviation of demographics for all replicates.
#'
#' @param funcSave_jobExpression This is a string expression that can be used to find elements of the experiment
#' being analysed.  It should be some robust unique string or regular expression.
#' @param func_saveFile This is the filepath and filename (ending in .RData please) to which the results of this
#' step will be saved.
#' @param func_experimentDir This is the filepath to the root directoy under which all your experimental files can
#' be found.
#' @param func_saveDir This is the directory to which output will be saved.
#' @param func_refFile This is the filepath to the reference file that contains information regarding all the
#' processed files for the rSHAPE experiment.
#' @param func_workEnvir This is an environment used to load files with the load function.  It's used to encapsulate
#' the loaded information to a controlled space.
#' @param func_objPrefix This is a character string for programatic naming of objects of this type.
#' @param func_sepString This is rSHAPE's sepString option but here to be passed into foreach
#' @param func_string_line_ofDescent This is rSHAPE's option of similar name to be passed into foreach
#' @param func_processedPattern This is rSHAPE's option of the similar name to be passed into foreach
#' @param func_sepLines This is rSHAPE's option of the similar name passed into foreach
#'
#' @section Note:
#' There is no example as this cannot work without a complete rSHAPE experiment to be analysed.
#'
#' @export
summarise_evolRepeatability <- function(funcSave_jobExpression,
                                        func_saveFile = getOption("shape_procExp_filenames")["repeatability"],
                                        func_experimentDir = getOption("shape_workDir"),
                                        func_saveDir = getOption("shape_postDir"),
                                        func_refFile = getOption("shape_procExp_filenames")[c("fileList","parameters")],
                                        func_workEnvir = new.env(),
                                        func_objPrefix = "Repeat_",
                                        func_sepString = getOption("shape_sepString"),
                                        func_string_line_ofDescent = getOption("shape_string_lineDescent"),
                                        func_processedPattern = getOption("shape_processedData_filePattern"),
                                        func_sepLines = getOption("shape_sepLines")){
  # This step can take time with larger experiments so we print a flag
  print("Starting repeatability")
  # We start by loading the reference file into this workspace, but if it does not exist we complain
  if(any(!file.exists(func_refFile))){
    stopError("Did not find the reference file needed for summarise_experimentParameters")
  }
  for(thisFile in func_refFile){
    load(thisFile, envir = func_workEnvir)
  }

  # We create an object that tracks the names of save files, objects and jobs for each set
  all_repSets <- sapply(unique(func_workEnvir$all_parmInfo$jobSet), function(thisSet){
                      tmp_thisName <- nameObject(func_inString = thisSet,
                                                 func_inPrefix = func_objPrefix)
                      return( list("objName"=tmp_thisName,
                                   "saveFile"= paste(sub(funcSave_jobExpression,"",
                                                         sub(getOption("shape_save_batchBase"),
                                                             tmp_thisName,
                                                             func_saveFile),
                                                         fixed=TRUE),sep=""),
                                   "setJobs" = unique(func_workEnvir$all_parmInfo$jobName[which(func_workEnvir$all_parmInfo$jobSet == thisSet)])) )

                    }, simplify = FALSE)

  # For all_uniqueSets, we find the ones where the files don't exist so those are the only we build
  tmp_missingSets <- names(all_repSets)[unlist(lapply(all_repSets,function(thisSet){
                                                    !file.exists(paste(getOption("shape_postDir"),thisSet[["saveFile"]],sep=""))
                                                  }))]
  # Now for each main job we'll gather some information and save it into it's own list object
  tmpCheck <- foreach(thisSet = tmp_missingSets, .combine="c",
                      .packages = c("rSHAPE", "DBI","RSQLite")) %dopar% {
    # We define an environment character string for the purposes of this foreach
    tmp_foreachString <- "_foreach_"

    tmpList <- vector(mode="list",length=length(all_repSets[[thisSet]][["setJobs"]]))
    names(tmpList) <- all_repSets[[thisSet]][["setJobs"]]
    # Ok here is the level at which we'll place our foreach looping to create a list of results
    for(thisJob in all_repSets[[thisSet]][["setJobs"]]){
      # We also create a parameter space for this job
      load(paste(func_experimentDir, thisJob,"/",thisJob,"_Parameters_1.RData",sep=""),
           envir = func_workEnvir)

      # This builds a connection required for fixing the absRank measurements
      connections_dataBase <- list("genotypeSpace" = reset_shapeDB(paste(func_experimentDir,thisJob,"/",
                                                                              func_workEnvir$runParameters$DataManagement$shape_fileName_dataBase[["genotypeSpace"]],
                                                                              sep=""),
                                                                    func_type = "connect"))

      func_landscapeCon <- connections_dataBase$genotypeSpace
      func_landTables <- RSQLite::dbListTables(func_landscapeCon)
      allow_backMutations <- func_workEnvir$runParameters$FitnessLandscape$shape_allow_backMutations
      genomeLength <- func_workEnvir$runParameters$Population$shape_genomeLength
      db_splitTables <- func_workEnvir$runParameters$DataManagement$shape_db_splitTables

      # We find the job files associated to thisJob
      tmpJobs <- func_workEnvir$all_dividedFiles[[thisSet]][which(grepl(thisJob, func_workEnvir$all_dividedFiles[[thisSet]]))]
      # I now build a list for storing the information of our replicate
      tmp_replicateInfo <- NULL
      if(length(tmpJobs) > 0){
        # Now for each and every replicate within this job we want to grab the information concerning transitions and
        # the fitness landscape of the parent to child.  We don't define .combine to gets list return
        # I need to create an outer list object which stores the information for each job
        tmp_replicateInfo <- NULL
        # Now we break the processing into chunks to be added up and collected
        tmp_chunkSize <- 100
        tmp_numChunks <- ceiling(length(tmpJobs)/tmp_chunkSize)
        for(thisChunk in 1:tmp_numChunks){
          tmpAdd <- NULL
          # Now for each and every replicate within this job we want to grab the:
          # (1) the population wide - through time - min, mean,max fitness. - I'll return the whole demoMat
          for(thisFile in tmpJobs[(1+((thisChunk - 1)*tmp_chunkSize)):min((thisChunk*tmp_chunkSize),length(tmpJobs))]){
            tmp_fileString <- nameEnviron(strsplit(strsplit(thisFile,"/")[[1]][length(strsplit(thisFile,"/")[[1]])],".",fixed=TRUE)[[1]][1],
                                          funcBase = tmp_foreachString)
            assign(tmp_fileString, new.env(), envir = func_workEnvir)
            load(thisFile, envir = func_workEnvir[[tmp_fileString]])
            # We extract what were the established lineages
            tmpEstablished <- func_workEnvir[[tmp_fileString]]$runDemographics$vec_estLineages
            # the final lineages which existed can be found by checking the last row of info_estLines[["lineDemo"]][,,as.character(tmpEstablished)]
            lineDemo_maxStep <- which.max(as.numeric(nameTable_step(rownames(func_workEnvir[[tmp_fileString]]$info_estLines$lineDemo),
                                                                    funcSplit = TRUE,
                                                                    func_sepString = func_sepString)))
            # We extract what were the established lineages
            tmp_final_estLineage <- data.frame("genotypeID"=tmpEstablished[which(func_workEnvir[[tmp_fileString]]$info_estLines$lineDemo[lineDemo_maxStep,"isEstablished",as.character(tmpEstablished)] == 1)])
            tmp_final_estLineage$binaryString <- rep("",nrow(tmp_final_estLineage))
            tmp_final_estLineage[,"binaryString"] <- unlist(lapply(tmp_final_estLineage[,"genotypeID"],function(x){
                                                                retrieve_binaryString(func_genotypeID = x,
                                                                                      func_landscapeCon = func_landscapeCon,
                                                                                      func_subNaming = db_splitTables)[,"binaryString"]
                                                              }))

            # The final dominant lineage is the transitioned lineage which had the largest population on the last step.  In the event of a tie
            # we look for the lineage which transitioned latest and break subsequent ties by highest fitness.  All to ensure a single return.
            tmp_transitionMat <- func_workEnvir[[tmp_fileString]]$runDemographics$transitionMat
            tmp_final_domLineage <- cbind("popSize"= func_workEnvir[[tmp_fileString]]$info_estLines$lineDemo[lineDemo_maxStep,"popSize",as.character(tmp_transitionMat[,"genotypeID"])],
                                          matrix(tmp_transitionMat[,c("Step","fitness","genotypeID")],
                                                 nrow=nrow(tmp_transitionMat), dimnames = list(NULL,c("Step","fitness","genotypeID"))))
            tmp_maxSized = which(tmp_final_domLineage[,"popSize"] == max(tmp_final_domLineage[,"popSize"]))

            # As a final safety, but unlikely requirement, if multiple genotypeID's transitioned at the same time, have the same popSize,
            # and have the same fitness, we simply take the first value.  This should never matter.
            tmp_final_domLineage <- if(length(tmp_maxSized) != 1){
                                      tmpReturn <- which(tmp_final_domLineage[tmp_maxSized,"Step"]==max(tmp_final_domLineage[tmp_maxSized,"Step"]))
                                      if(length(tmpReturn) > 1){
                                        tmpReturn <- tmpReturn[which(tmp_final_domLineage[tmp_maxSized[tmpReturn],"fitness"]==max(tmp_final_domLineage[tmp_maxSized[tmpReturn],"fitness"]))][1]
                                      }
                                      tmp_final_domLineage[tmp_maxSized[tmpReturn],"genotypeID"]
                                    } else {
                                      tmp_final_domLineage[tmp_maxSized,"genotypeID"]
                                    }

            # We now build this information into an object shaped as intended for use.
            tmp_final_domLineage <- c("genotypeID"=unname(tmp_final_domLineage),
                                      "binaryString"="")
            tmp_final_domLineage["binaryString"] <- tmp_final_estLineage[which(tmp_final_estLineage[,"genotypeID"] == tmp_final_domLineage["genotypeID"]),"binaryString"]

            # We now also want to track the line of descent for the dominant lineage and it's first transition
            tmp_line_of_descent <- as.character(func_workEnvir[[tmp_fileString]]$info_estLines$end_Lines_of_Descent[[tmp_final_domLineage["genotypeID"]]])
            # We now try and extract the first transition in this line of descent
            tmp_dom_transitions <- strsplit(tmp_line_of_descent,func_string_line_ofDescent)[[1]]
            # Now if there was no transition we should be left with the WT genotype ID of "0" as the onyl element... otherwise we have a first step.
            tmp_firstTransition <- if(length(tmp_dom_transitions) == 1){
                                      paste(rep(tmp_dom_transitions[1],2),collapse=func_string_line_ofDescent)
                                    } else {
                                      paste(tmp_dom_transitions[1:2],collapse=func_string_line_ofDescent)
                                    }
            # We now add more information to the domLineage vector
            tmp_final_domLineage <- c(tmp_final_domLineage,
                                      "line_ofDescent"= paste(tmp_line_of_descent,collapse= func_sepLines),
                                      "firstStep"=tmp_firstTransition)

            # We build now the transition mat into something which shows the steps taken on the fitness landscape
            # In order to know the transitions, we look at the transition matrix and for a lineage, but starting with the first
            # that is found in the tranisition matrix which was not established on Step_0
            tmp_initialEst <- dimnames(func_workEnvir[[tmp_fileString]]$info_estLines$lineDemo)[[3]]
            tmp_initialEst <- unique(c(tmp_initialEst[which.max(func_workEnvir[[tmp_fileString]]$info_estLines$lineDemo[nameTable_step(0, func_sepString = func_sepString),
                                                                                                                        "popSize", tmp_initialEst])],
                                       tmp_initialEst[which(func_workEnvir[[tmp_fileString]]$info_estLines$lineDemo[nameTable_step(0, func_sepString = func_sepString),
                                                                                                                    "isEstablished", tmp_initialEst] == 1)]))
            # If there are any genotypes which were not in the initial set we go further, otherwise we're returning some null stuff
            tmp_returnMat <- NULL
            if(any(!is.element(tmp_transitionMat[,"genotypeID"],as.numeric(tmp_initialEst))) && length(tmp_initialEst) > 0){
              # Ok this means that some transition occured.  I now just want to find what transitions occured
              tmp_landscapeTopology <- func_workEnvir[[tmp_fileString]]$info_estLines$landscapeTopology
              tmpSteps <- strsplit(rownames(tmp_landscapeTopology),func_string_line_ofDescent)
              # This is hte first row of our transitionMat at which there was a change
              tmp_minRow <- min(which(!is.element(tmp_transitionMat[,"genotypeID"],as.numeric(tmp_initialEst))))
              # Ok this is what will be assigned as the return matrix, it should inform on the
              # genotypeID_1, fitness_1, absRank, numMuts_1, transition, transitionStep
              for(thisRow in tmp_minRow:nrow(tmp_transitionMat)){
                tmpID <- tmp_transitionMat[thisRow, ]
                # We look to the landscapeTopology object to get transition information, what we are looking for
                # is an instance where the second element (offspring) is the same as our tmpID, and the progenitor
                # is some element of the pedigree for our lineage.
                tmp_possProgenitors <- which(sapply(tmpSteps,function(x){
                                                  as.numeric(x[2]) == tmpID["genotypeID"]
                                                }))
                # If this transition is simply due to a shift in dominance, rather than a new type whic has emerged,
                # due to mutation, then there may be no progenitors.  In this case we skip it as it is not a mutational step.
                if (length(tmp_possProgenitors) == 0){ next }
                tmp_possProgenitors <- matrix(tmp_landscapeTopology[tmp_possProgenitors,],
                                              nrow = length(tmp_possProgenitors),
                                              dimnames=list(rownames(tmp_landscapeTopology)[tmp_possProgenitors],
                                                            colnames(tmp_landscapeTopology)))
                tmp_possProgenitors <- cbind(tmp_possProgenitors,
                                             t(sapply(strsplit(rownames(tmp_possProgenitors), func_string_line_ofDescent), function(x){
                                                    x <- as.numeric(x)
                                                    return( c("progenitorID"=x[1],
                                                              "offspringID"=x[2]) )
                                                  })))
                # We return now the information
                tmp_returnMat <- rbind(tmp_returnMat,
                                       data.frame("progenitor_genotypeID"= tmp_possProgenitors[,"progenitorID"],
                                                  "offspring_genotypeID"= tmp_possProgenitors[,"offspringID"],
                                                  "progenitor_fitness"= tmp_possProgenitors[,"progenitor_fitness"],
                                                  "offspring_fitness"= tmp_possProgenitors[,"offspring_fitness"],
                                                  "absRank" = tmp_possProgenitors[,"absRank"],
                                                  "hoodSize" = tmp_possProgenitors[,"hoodSize"],
                                                  "hoodMin"= tmp_possProgenitors[,"hoodMin"],
                                                  "hoodMax"= tmp_possProgenitors[,"hoodMax"],
                                                  "num_altPaths"=  tmp_possProgenitors[,"num_altPaths"],
                                                  "relFit_altPaths" = tmp_possProgenitors[,"relFit_altPaths"],
                                                  "prop_maxFit" = tmp_possProgenitors[,"prop_maxFit"],
                                                  "progenitor_numMuts"= tmp_possProgenitors[,"progenitor_numMuts"],
                                                  "offspring_numMuts"= tmp_possProgenitors[,"offspring_numMuts"],
                                                  "transition"=rownames(tmp_possProgenitors),
                                                  "transitionStep"= tmpID["Step"],
                                                  row.names = rownames(tmp_possProgenitors),
                                                  stringsAsFactors = FALSE) )
              }
              # Now it's possible the same genotype transition occured multiple times which woould result in
              # replicated rows OTHER than the transitionStep value, we search for this
              tmp_checkCols <- which(colnames(tmp_returnMat) != "transitionStep")
              tmpDuplicated <- duplicated(tmp_returnMat[,tmp_checkCols])
              if(any(tmpDuplicated)){
                # We find the paired rows for each duplicated, I setup a working reference matrix due
                # to observed issues comparing rows with my apply method.
                tmp_workMat <- matrix(gsub("[[:space:]]","",as.character(unname(unlist(tmp_returnMat[,tmp_checkCols])))),
                                      ncol = ncol(tmp_returnMat[,tmp_checkCols]))
                tmpDuplicated <- lapply(which(tmpDuplicated),function(tmpRow){
                  which(apply(tmp_workMat,MARGIN=1,function(x){
                    all(x == tmp_workMat[tmpRow,tmp_checkCols])
                  }))
                })
                # So for each unique set in the list we trim our return matrix
                tmp_removeRows <- NULL
                for(thisSet in unique(tmpDuplicated)){
                  tmp_returnMat[thisSet[1],"transitionStep"] <- paste(tmp_returnMat[thisSet,"transitionStep"],collapse= func_sepString)
                  tmp_removeRows <- c(tmp_removeRows,thisSet[-1])
                }
                if(!is.null(tmp_removeRows)){ tmp_returnMat <- tmp_returnMat[-tmp_removeRows,] }
              }
            } else {
              # this is the rather null case of no steps being taken
              tmp_returnMat <- data.frame("progenitor_genotypeID"= tmp_transitionMat[1,"genotypeID"],
                                          "offspring_genotypeID"= NA,
                                          "progenitor_fitness"=tmp_transitionMat[1,"fitness"],
                                          "offspring_fitness"= NA,
                                          "absRank" = NA,
                                          "hoodSize" = NA,
                                          "hoodMin"= NA,
                                          "hoodMax"= NA,
                                          "num_altPaths"=  NA,
                                          "relFit_altPaths" = NA,
                                          "prop_maxFit" = NA,
                                          "progenitor_numMuts"= tmp_transitionMat[1,"numMuts"],
                                          "offspring_numMuts"= NA,
                                          "transition"=paste(tmp_transitionMat[1,"genotypeID"]),
                                          "transitionStep"=0,
                                          stringsAsFactors = FALSE)
            }  # This closes out building the tmp_returnMat object

            # Ok we'll now return this information
            tmpReturn <- list(list("final_estLineages"= tmp_final_estLineage,
                                   "final_domLineage"= tmp_final_domLineage,
                                   "transitions"= tmp_returnMat))
            names(tmpReturn)[length(tmpReturn)] <- nameEnviron(func_Index = tmp_fileString,
                                                               funcSplit = TRUE,
                                                               funcBase = tmp_foreachString)

            tmpAdd <- c(tmpAdd, tmpReturn)
            # Now we remove the environment of the file and it's data
            rm(list=c(tmp_fileString),pos= func_workEnvir)
          }
          tmp_replicateInfo <- c(tmp_replicateInfo, tmpAdd)
          rm(tmpAdd)
          gc()
        } # This closes out the for loop gathering information about the replicates for this job
      } # This closes out the conditional if that there is at least one file for this job
      tmpList[[thisJob]] <- tmp_replicateInfo
      # We remove the objects we created in this section
      for(thisConnection in connections_dataBase){ RSQLite::dbDisconnect(thisConnection) }
      # Clean up space
      rm(tmp_replicateInfo, connections_dataBase,func_landscapeCon,
         allow_backMutations, genomeLength,db_splitTables,func_landTables)
    } # This closes out the for loop for different jobs

    # We can now save this object to a file
    assign(all_repSets[[thisSet]][["objName"]], tmpList, pos= func_workEnvir)
    save(list = all_repSets[[thisSet]][["objName"]],
         file= all_repSets[[thisSet]][["saveFile"]],
         envir = func_workEnvir)
    rm(list= all_repSets[[thisSet]][["objName"]], envir = func_workEnvir)
    # Clean up space
    gc()
    # We return a confirmation that the job completed
    return(  nameObject(func_inString = all_repSets[[thisSet]][["objName"]],
                        func_inPrefix = func_objPrefix,
                        func_splitStr = TRUE)  )
  } # This closes out the foreach loop

  # We now check if all the jobs completed
  tmp_checkCompleted <- sapply(tmp_missingSets,is.element,"set"=tmpCheck)
  if(all(tmp_checkCompleted)){
    # We now save the object which stores all the repeatability data object names and files locations
    save(all_repSets,file= func_saveFile)
  } else {
    stop(paste("There was a problem with the sets of: ",
               paste(names(tmp_checkCompleted)[which(!tmp_checkCompleted)],collapse=" ",sep=" "),
               sep=""))
  }

  # We silently return nothinig
  invisible( NULL )
}
