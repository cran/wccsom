.First.lib <- function(lib, pkg) library.dynam("wccsom",pkg,lib)
.Last.lib <- function(libpath) library.dynam.unload("wccsom", libpath)
