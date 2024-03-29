% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poly.R
\name{twoStage.poly}
\alias{twoStage.poly}
\title{Polynomial fit for seasonality estimation}
\usage{
twoStage.poly(d, freq, nrDays)
}
\arguments{
\item{d}{data whose seasonality is to be estimated.}

\item{freq}{period of the seasonality to be estimated.}

\item{nrDays}{number of days in a week in the input data.}
}
\description{
Estimates seasonality using polynomial fitting. Currently
    this works only for daily and weekly level of seasonality.
}
\examples{

}
