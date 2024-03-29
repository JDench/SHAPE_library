% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sourceSHAPE.R
\name{summarise_evolRepeatability}
\alias{summarise_evolRepeatability}
\title{This function will use output from summarise_experimentFiles and summarise_experimentParameters
to help with expectations concerning run output and handling.  This will save an RData file which
will contain one object: all_popSets, which is a list of relevant control information about I/O
and then a series of other RData files which contain the demographics information as a matrix with the
mean and standard deviation of demographics for all replicates.}
\usage{
summarise_evolRepeatability(funcSave_jobExpression,
  func_saveFile = getOption("shape_procExp_filenames")["repeatability"],
  func_experimentDir = getOption("shape_workDir"),
  func_saveDir = getOption("shape_postDir"),
  func_refFile = getOption("shape_procExp_filenames")[c("fileList",
  "parameters")], func_workEnvir = new.env(),
  func_objPrefix = "Repeat_",
  func_sepString = getOption("shape_sepString"),
  func_string_line_ofDescent = getOption("shape_string_lineDescent"),
  func_processedPattern = getOption("shape_processedData_filePattern"),
  func_sepLines = getOption("shape_sepLines"))
}
\arguments{
\item{funcSave_jobExpression}{This is a string expression that can be used to find elements of the experiment
being analysed.  It should be some robust unique string or regular expression.}

\item{func_saveFile}{This is the filepath and filename (ending in .RData please) to which the results of this
step will be saved.}

\item{func_experimentDir}{This is the filepath to the root directoy under which all your experimental files can
be found.}

\item{func_saveDir}{This is the directory to which output will be saved.}

\item{func_refFile}{This is the filepath to the reference file that contains information regarding all the
processed files for the rSHAPE experiment.}

\item{func_workEnvir}{This is an environment used to load files with the load function.  It's used to encapsulate
the loaded information to a controlled space.}

\item{func_objPrefix}{This is a character string for programatic naming of objects of this type.}

\item{func_sepString}{This is rSHAPE's sepString option but here to be passed into foreach}

\item{func_string_line_ofDescent}{This is rSHAPE's option of similar name to be passed into foreach}

\item{func_processedPattern}{This is rSHAPE's option of the similar name to be passed into foreach}

\item{func_sepLines}{This is rSHAPE's option of the similar name passed into foreach}
}
\description{
This function will use output from summarise_experimentFiles and summarise_experimentParameters
to help with expectations concerning run output and handling.  This will save an RData file which
will contain one object: all_popSets, which is a list of relevant control information about I/O
and then a series of other RData files which contain the demographics information as a matrix with the
mean and standard deviation of demographics for all replicates.
}
\section{Note}{

There is no example as this cannot work without a complete rSHAPE experiment to be analysed.
}

