% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lf_to_tdf.R
\name{data_manipulation}
\alias{data_manipulation}
\title{Data Manipulation for Marker Data}
\usage{
data_manipulation(marker, threshold = 500)
}
\arguments{
\item{marker}{A data frame containing marker data, where each row represents
a marker and each column represents a sample.}

\item{threshold}{The threshold value for peak height. Peaks below this
threshold will be replaced with 0.}
}
\value{
A formatted data frame where at least one peak for each sample is
greater than the specified threshold.
}
\description{
This function ensures that at least one peak for each sample is greater than
a specified threshold (default: 500) and then formats the data frame for
the next steps in the analysis.
}
\examples{

marker_data <- data.frame(
Sample1 = c(400, 600, 700,0),
Sample2 = c(450, 550, 480,0),
Sample3 = c(300, 200, 400,200),
Sample4 = c(0,0,0,0),
row.names=c(185,188,191,194)
)
data_manipulation(marker_data,threshold=500)
}
