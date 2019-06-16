# This is the SHAPE body file but should not itself be called directly.  Instead use the farming file to actually create and submit jobs for running.
# The program will simulate the population demographics and genotype evolution
# A broad range of parameters are currently supported to allow comparison of dynamics across various models and assumptions

############ DEPENDENCIES ############
# It uses the parameter input file named by the sourceParms object
# Running this script should pass along elements which can be captured by the commandArgs() function
# which should return at least a 'getOption("shape_thisRep")' argument, and possibly an 'shape_outDir'
# then just the libraries within the source files called below

#### NOTE: This is only established to work for when mutational steps are single mutants, this could be fixed if we address how neighbours are found (see function)

# We tidy up our space and then load libraries
rm(list=ls())

# We load the SHAPE library and then initialise the basic parameters.
library(rSHAPE)
defineSHAPE()

# This is the file path to the parameters that will be used for this run, this will over-write many of the default parameters.
sourceParms <- "C:/Users/Jonathan/Documents/Programming/MyScripts/SHAPE/compileSHAPE_aux/SHAPE_parameters.v.1.r"
if(!file.exists(sourceParms)){
     stop(paste("Could not find the parameters file located at:",sourceParms,sep=" "))
     q(save="no")
} else {
     source(file=sourceParms)
}

# We create a vector that controls which replicates this job should be cycling across,
# if we replicate internally this is between the getOption("shape_thisRep") and shape_maxReplicates value
# otherwise it's simply the getOption("shape_thisRep")
workingReplicates <- if(getOption("shape_externalSelfing")){
                             getOption("shape_thisRep")
                        } else {
                             seq(getOption("shape_thisRep"), getOption("shape_maxReplicates"),by=1)
                        }

# We set the strings as factors option
options(stringsAsFactors = getOption("shape_stringsAsFactors"))

# This calls SHAPE to run given the parameters loaded.
runSHAPE()

q(save="no")



