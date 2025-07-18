% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Help_functions.R
\name{matrix_vector}
\alias{matrix_vector}
\title{Extract Upper Triangular Values from a Matrix as a Vector}
\usage{
matrix_vector(mat, ind = NULL)
}
\arguments{
\item{mat}{A numeric square matrix.}

\item{ind}{Optional integer vector of indices to exclude from the extracted upper triangular values.}
}
\value{
A numeric vector containing the upper triangular elements of \code{mat}, with specified indices removed if provided.
}
\description{
Extracts the elements above the main diagonal (upper triangular part) of a square matrix as a vector.
Optionally, specified indices can be excluded from the extracted values.
}
