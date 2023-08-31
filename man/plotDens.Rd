\name{plotDens}
\docType{methods}
\alias{plotDens}
\alias{plot, flowFrame, CellPopulation-method}
\title{
  Plot flow cytometry data with density-based colors
}

\description{
    Generate a scatter dot plot with colors based on the distribution of
    the density of the provided channels.
}
\usage{
plotDens(obj, channels,node=NA ,col, main, xlab, ylab, xlim,ylim, pch=".", density.overlay=c(FALSE,FALSE),count.lim=20, dens.col=c("grey48","grey48"),cex=1,verbose=TRUE,
dens.type=c("l","l"),transparency=1, adjust.dens=1,show.contour=F, contour.col="darkgrey", ...)
}
\arguments{
  \item{obj}{
  	a 'FlowFrame','GatingHierarchy' or 'cellPopulation' object.
  }
  \item{channels}{
    a vector of two channel names or their corresponding indices in the 'flow.frame'.
  }
  \item{node}{
   Only used if gatingHierarchy is provided as an input.It subtractes the flowFrame from the node of gatingHierarchy. see getNodes() in flowWorkspace.
  }
  
  \item{col}{
 	A specification for the default plotting color: see '?par'. 	
  }
  \item{main}{
  	an overall title for the plot: see '?plot'
  }
  \item{xlab}{
  	a title for the x axis: see '?plot'
  }
  \item{ylab}{
  	a title for the y axis: see '?plot'  	
  }
\item{xlim}{
  	a range for the x axis: see '?plot'
  }
\item{ylim}{
  	a range for the y axis: see '?plot'
  }
  \item{cex}{
  	A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.: see '?par'
  }
  \item{density.overlay}{
  	Logical vector of length 2, to plot density overlays on the x and y axes. Default is c(FALSE,FALSE).  	
  }
 \item{dens.type}{
  	2-character string giving the type of plot desired, see type in '?plot'
  }
 \item{dens.col}{
  	2-character string giving the color of plot desired for density curves.
  }
\item{transparency}{
      Transparency of the bi-variate plot, to see the densitu curves in the background. The lower it is, the more transparent the plot is.
  	
  }
\item{show.contour}{
      Default is FALSE. It add the contourLines to plot.
}
\item{contour.col}{
    Color for contourLines. Default is darkgrey.
}
 \item{adjust.dens}{
    The smoothness of density in [0,Inf] to be used in density(.). The default value is 1 and should not be changed unless necessary
  }
\item{count.lim}{
   Cutoff for number of events to set color. Default is 20. Samples with less than 20 cells will be plotted in black.
}
\item{verbose}{
   Default is True. It will add that the sample has 0 cells in the plot title.
}
  \item{pch}{
  	Either an integer specifying a symbol or a single character to be used as the default in plotting points: see '?par'.
  }
  \item{...}{
    can be used to provide desired arguments for the plot() function used to
    plot the output results.
  }
}


\value{
  a scatter dot plot  with density-based colors, along with density overlays if desired. Set xlim and ylim when plotting if you would like to have all your plots to have same range on the axes (specially when density.overlay=TRUE)
}

\author{
Mehrnoush Malek \email{<mmalekes@bccrc.ca>}
  Jafar Taghiyar \email{<jtaghiyar@bccrc.ca>}
  
}

\examples{

# data_dir <- system.file("extdata", package = "flowDensity")
# load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
#Plot CD3 vs. CD19 to see the distribution of cell populations and their density
#plotDens(f,c("V450-A","PerCP-Cy5-5-A"))
}
\keyword{ Automated gating }

\keyword{ FlowCytData }