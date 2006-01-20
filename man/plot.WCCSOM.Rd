\name{plot.WCCSOM}
\alias{plot.WCCSOM}
\title{Plot WCCSOM object}
\description{
  Plot self-organising map, obtained from function WCCSOM. Several types
  of plots are supported.
}
\usage{
plot.WCCSOM(x, type = c("changes", "codes", "counts", "mapping",
                        "property", "quality"),
            classif = x$unit.classif, obj, labels=NULL, pchs=NULL,
            bgcol=heat.colors(25)[6:25],
            zlim=NULL, property, heatkey=TRUE, \dots)
}
\arguments{
  \item{x}{SOM object.}
  \item{type}{Type of plot. (Wow!)}
  \item{classif}{classification object, as returned by \code{wccassign}.
    Only needed if \code{type} equals '"mapping"', '"quality"',
    and '"counts"'. Also a vector of class numbers may be given. If the
    network was trained with \code{keep.data} equal to TRUE, then the
    object already contains this information.}
  \item{obj}{A new object to be mapped. Only needed if \code{type}
    equals '"property"'. }
  \item{labels}{labels to plot when \code{type} equals
    '"mapping"'.}
  \item{pchs}{plotting symbols to use when \code{type} equals
    '"mapping"'.}
  \item{bgcol}{colors to use as unit background for '"codes"',
    '"property"' and '"counts"' plotting types.}
  \item{zlim}{Optional range for color coding of unit backgrounds.}
  \item{property}{Values to use if the '"property"' plotting type. Can be used
    for colouring units in general.}
  \item{heatkey}{Whether or not to generate a heatkey at the left side
    of the plot in the '"quality"', '"counts"', and '"property"'
    plotting types.}
  \item{\dots}{Other graphical parameters, e.g. colours of labels in the
    '"mapping"' plotting type.}
}
\details{Several different types of plots are supported:
  \describe{
    \item{'"changes"'}{Shows the mean change in similarity to the best
      matching codebook vector for each epoch. Since codebook vectors
      become more similar to the data that are mapped to them, the
      changes should always be positive. Upon convergence, the changes
      should be very small.}
    \item{'"codes"'}{Shows the codebook vectors.}
    \item{'"counts"'}{Shows the number of objects mapped to the
      individual units. Empty units are depicted in gray.}
    \item{'"mapping"'}{Show where a set of objects is mapped. It
      needs a '"labels"' argument: a string name for each object.}
    \item{'"property"'}{The standard application is to show the similarity
      of one object to all units in the network; this is done by supplying
      only the '"obj"' argument. The WCC with all units in the network
      is calculated and used to colour the units. The unit with the
      highest WCC value is indicated with blue crosshairs.
      In other cases, these similarities may already be precomputed and
      then they can be provided as the '"property"' argument; the '"obj"'
      argument is unnecessary in that case. Also other quantities may be
      used to colour the units: see the example of distances below. If
      the '"property"' argument is given, no blue crosshairs are shown.}
    \item{'"quality"'}{Shows the units coloured according to the mean
      agreement (WCC) of mapped objects to the unit vector. A colour key is
      plotted on the left. The variation in the WCCs of the mapped
      objects is indicated by the blue line: if it is pointing downwards it
      indicates low variation, if pointing upwards large variation.
      Empty units are depicted in gray.}
  }
}
\value{If \code{type} equals '"property"', the wcc values for all units
  are returned.}

\author{Ron Wehrens}
\seealso{\code{\link{WCCSOM}}}
\examples{
data(cepha)
gr <- somgrid(3, 3, "hexagonal")
set.seed(7)
x <- WCCSOM(cepha$patterns, grid=gr, trwidth=20, rlen=100)

par(mfrow=c(3,2))
plot(x, type="codes", main="Codebook vectors")
plot(x, type="changes", main="Convergence")
plot(x, type="counts", main="Mapping counts")
plot(x, type="quality", main="Mapping quality")
plot(x, type="mapping", main="Mapping", 
     labels=cepha$class.names, col=as.integer(factor(cepha$class.names)))
plot(x, type="mapping", main="Mapping",
     pchs=as.integer(factor(cepha$class.names)),
     col=as.integer(factor(cepha$class.names)))

obj1.wccs <- plot(x, type="property", obj=cepha$patterns[1,])
obj1.wccs

par(mfrow=c(1,2))
unit.dists <- unit.distances(x$grid, toroidal=FALSE)
plot(x, type="property", property=unit.dists[1,],
     main="Distances to unit 1", zlim=c(0,4))
unit.dists <- unit.distances(x$grid, toroidal=TRUE)
plot(x, type="property", property=unit.dists[1,],
     main="Toroidal distances to unit 1", zlim=c(0,4))
}
\keyword{classif}
