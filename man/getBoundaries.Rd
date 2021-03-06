% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Segmentation.R
\name{getBoundaries}
\alias{getBoundaries}
\title{Candidate Differential Boundaries}
\usage{
getBoundaries(B, h, w, t)
}
\arguments{
\item{B}{Matrix containing the following columns:
1 & 2. Chromosome, 3. Focus, 4. Group, 5. Border score}

\item{h}{Window half-width.}

\item{w}{Minimum boundary separation.}

\item{t}{Threshold local Z-score for candidate boundaries.}
}
\description{
Candidate Differential Boundaries
}
