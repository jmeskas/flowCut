\name{removeLowDensSections}
\alias{removeLowDensSections}

\title{
    Remove low density regions of a flowFrame.
}

\description{
    Remove low density regions of a flowFrame.
}

\usage{
    removeLowDensSections (f, LowDensityRemoval=0.1, Verbose=FALSE)
}
\arguments{
    \item{f}{
        A single flowFrame to be processed.
    }
    \item{LowDensityRemoval}{
        An numeric value between 0 and 1 that specifies the proportion cut-off of the low density region. All events with less density than LowDensityRemoval will be removed. The default is 0.1.
    }
    \item{Verbose}{
        A logical value that dictates if computational information is printed while the code is running. The default is FALSE.
    }
}

\details{
    removeLowDensSections removes low density regions of a flowFlame. The flowFrame must have a time channel.
}

\value{
    A list containing two elements. The first element, $frame, is the flowFrame file returned after removeLowDensSections has removed the low density sections. The second element, $rem.ind, is a vector containing indices of events being removed.
}

\author{
    Maintainer: Justin Meskas <jmeskas@bccrc.ca>

    Authors: Justin Meskas, Sherrie Wang
}


\examples{
data(flowCutData)
res_lowDens <- removeLowDensSections(flowCutData[[1]], LowDensityRemoval=0.1)
library(flowDensity)
plotDens(flowCutData[[1]], c("Time","FL1-H"))
points(exprs(flowCutData[[1]])[res_lowDens$rem.ind, c("Time","FL1-H")], pch=".", col="grey")
}

\keyword{FlowCytData}