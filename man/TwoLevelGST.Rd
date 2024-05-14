% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pop_Gen_Functions.R
\name{TwoLevelGST}
\alias{TwoLevelGST}
\title{Calculate Two-Level GST}
\usage{
TwoLevelGST(G = matrix)
}
\arguments{
\item{G}{A square matrix representing a genetic distance matrix.}
}
\value{
A list containing the components of two-level GST including within-group gene identity, between-group gene identity, and GST values.
}
\description{
This function calculates two-level GST (Nei's standard genetic distance) measure from a genetic distance matrix.
}