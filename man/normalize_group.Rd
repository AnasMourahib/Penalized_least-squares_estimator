% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Help_functions.R
\name{normalize_group}
\alias{normalize_group}
\title{Normalize Elements of a Vector by Groups}
\usage{
normalize_group(v, group_size)
}
\arguments{
\item{v}{A numeric vector to be normalized.}

\item{group_size}{An integer specifying the size of each group.}
}
\value{
A numeric vector of the same length as \code{v}, where each group of \code{group_size} elements is normalized to sum to 1 (unless the group sum is zero).
}
\description{
Splits a vector into groups of equal size and normalizes each group so that its elements sum to 1. If a group's sum is zero, it is left unchanged.
}
