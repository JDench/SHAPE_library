% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sourceSHAPE.R
\name{fitnessLandscape}
\alias{fitnessLandscape}
\title{This function will calculate the fitness values for genotypes being newly recorded to the fitness landscape.}
\usage{
fitnessLandscape(tmpGenotypes, tmp_focalFitness, landscapeModel = "HoC",
  tmp_ancestralFitness = getOption("shape_const_ancestFitness"),
  tmp_weightsRMF = getOption("shape_const_RMF_theta"),
  tmp_optimaRMF = getOption("shape_const_RMF_globalOptima"),
  tmp_correlationsNK = getOption("shape_const_NK_interactionMat"),
  tmp_const_numInteractionsNK = getOption("shape_const_numInteractions"),
  tmp_NK_ancestDep = getOption("shape_const_DepbySite_ancestFitness"),
  relativeFitness = TRUE,
  func_genomeLength = getOption("shape_genomeLength"),
  func_distribution = getOption("shape_constDist"),
  func_distParameters = getOption("shape_const_distParameters"),
  func_distAsS = getOption("shape_const_distAsS"),
  func_sepString = getOption("shape_sepString"))
}
\arguments{
\item{tmpGenotypes}{This is a vector of the binaryString values that represent the genotype(s) for which you want
to calculate new fitness values.}

\item{tmp_focalFitness}{This argument has different meaning depending upon the fitness landscape model being simulated.
It can be a vector of fitness values, a matrix, a single value, etc...}

\item{landscapeModel}{This is the character string that defines the fitness landscape model being simulated in this SHAPE run.
At present it can be one of: \strong{Additive, Fixed, HoC, NK, RMF}}

\item{tmp_ancestralFitness}{This is the fitness value of the pure WT genotype, it does not always have meaning.}

\item{tmp_weightsRMF}{This is the weighting of the constant/deterministic term calculated in the RMF fitness landscape equation.}

\item{tmp_optimaRMF}{This is the binary string genotype of the optimal genotype in the current RMF fitness landscape.  It needn't yet have
been yet explored, it is simply the genotype that will be the deterministic global optimum.}

\item{tmp_correlationsNK}{This is the matrix of fitness values and interactions between mutational states for the NK fitness lanscape model}

\item{tmp_const_numInteractionsNK}{This is the "K" value of the NK fitness landscape value and represents the number of other sites
correlated to the fitness of a focal site.}

\item{tmp_NK_ancestDep}{This is the fitness value of the WT mutant for an NK fitness landscape, it is passed as a computational
ease so that it needn't be calculated each time this function is called.}

\item{relativeFitness}{This is a logical toggle controlling if the fitness values returned should be relative fitness values}

\item{func_genomeLength}{This is the genome length of individuals.}

\item{func_distribution}{This is a character string representing which of the allowed distribution functions can be called
for draws of stochastic values when calculating fitness values.  See fitnessDist for those implemented.}

\item{func_distParameters}{This is a vector of the ordered distribution parameters expected by the distribution
referenced by func_distribution}

\item{func_distAsS}{This is a logical toggle to control in the final returned values should be considered as selection coefficients,
which is achieved by subtracting the calculated value by 1.}

\item{func_sepString}{This is a character string used for collapsing vectors of information, and expanding the collpased information back into
a vector of values.}
}
\value{
A vector of fitness values to be assgined for each of the newly explored genotypes defined in the vector tmpGenotypes
}
\description{
This function will calculate the fitness values for genotypes being newly recorded to the fitness landscape.
}
\section{Note}{

There is no example as this does not have meaning outisde of a runSHAPE call.
}

