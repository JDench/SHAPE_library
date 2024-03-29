% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sourceSHAPE.R
\name{nameTable}
\alias{nameTable}
\title{This is a standardising function which allows SHAPE to programatiically name tables for the fitness landscape OR split
a named table and extract the embedded information from its naming.}
\usage{
nameTable(func_tmpMutations, func_tmpIndex = NULL,
  func_baseString = getOption("shape_string_tableNames"),
  func_sepString = getOption("shape_sepString"),
  func_splitName = FALSE,
  func_subNaming = getOption("shape_db_splitTables"))
}
\arguments{
\item{func_tmpMutations}{Integer value(s) for the number of mutations to be expected in mutants stored within the named tables.}

\item{func_tmpIndex}{An optinal element that will be used to insert a unique vector ID}

\item{func_baseString}{This is the standard prefix character string used in table naming.}

\item{func_sepString}{This is a character string used to collapse vectors of characters.}

\item{func_splitName}{A logical toggle to control if this function is splitting a named table or not.  So, FALSE (default)
means we're creating a table name whereas TRUE is splitting a named table into it's parts.}

\item{func_subNaming}{This is a logical which controls if the tables which report on all genotypes with X mutations should be
forced into a single table or it SHAPE is allowed to split them into multiple tables.}
}
\value{
If func_splitName is TRUE, then a vector of table names is returned, it would be best practice to not assume recycling of passed
elements and so pass equally lengthed vectors as input.  If FALSE, we split the table and return the data detailing the number of mutations
which ought to be present for genotypes stored in the named table.
}
\description{
This is a standardising function which allows SHAPE to programatiically name tables for the fitness landscape OR split
a named table and extract the embedded information from its naming.
}
\examples{
# This creates a table name in a standard way, it can also split table names to extract info.
defineSHAPE()
nameTable(2,1,"myTest","_",FALSE,FALSE)
nameTable("myTest_2",func_splitName = TRUE)

}
