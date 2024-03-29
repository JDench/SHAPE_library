% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sourceSHAPE.R
\name{runProcessing}
\alias{runProcessing}
\title{This is a wrapper function to process a SHAPE run and extract meaningful summary information.}
\usage{
runProcessing(func_saveFile, func_subNaming, func_stepsCon,
  func_landscapeCon, func_hoodCon, func_estProp, func_size_timeStep,
  func_processObjects = getOption("shape_processedObjects"),
  func_hoodPriority = getOption("shape_const_hoodDepth"))
}
\arguments{
\item{func_saveFile}{This is the filepath where the SHAPE run processed objects are to be saved.}

\item{func_subNaming}{This is a logical which controls if the tables which report on all genotypes with X mutations should be
forced into a single table or it SHAPE is allowed to split them into multiple tables.}

\item{func_stepsCon}{This is the filepath to an SQLite database storing information for the stepwise changes of a SHAPE run.}

\item{func_landscapeCon}{This is the filepath to an SQLite database storing information for the complete explored and neighbouring fitness landscape of a SHAPE run.}

\item{func_hoodCon}{This is the filepath to an SQLite database storing information for high priority mutational neighbourhood information}

\item{func_estProp}{This value is used to define the threshold size required for a population before it is considered established.}

\item{func_size_timeStep}{This is the proportion of a standard biological generation being considered to be within a single time step.}

\item{func_processObjects}{This is a vector of character strings which define the names of what objects will be produced and creates a global objects.  DO NOT CHANGE THESE VALUES.}

\item{func_hoodPriority}{This is an object to control which strains we get deep neighbourhood information for
It should be one of \strong{"none","limited","priority","full"}
setting this higher will cost more and more in post analysis runtime.}
}
\value{
This returns a string vector stating the result of trying to process for the specified filepath.
}
\description{
This is a wrapper function to process a SHAPE run and extract meaningful summary information.
}
\section{Note}{

There is no example as this cannot work outisde of a runSHAPE call, it requires data produced by the simulation experiment.
}

