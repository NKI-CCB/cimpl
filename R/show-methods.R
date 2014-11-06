# Show methods of the cimpl package
# 
# Author: Jelle ten Hoeve
###############################################################################

setMethod("show", signature("cimplAnalysis"), function(object) {
	cat("Cimpl analysis object\n")
	cat(paste("\tData:", object@data.name, "\n"))
	cat(paste("\tChromosomes:", paste(object@chromosomes, collapse=', '),"\n"))
	cat(paste("\tScales:", paste(object@scales, collapse=', '),"\n"))
	cat(paste("\tSystem:", object@system, "\n"))
	cat(paste("\tSpecificity pattern:", object@specificity.pattern, "\n"))
	cat(paste("\tLocal hopping correction:", object@lhc.method, "\n"))
	cat(paste("\tNumber of iterations:", object@n_iterations, "\n"))
	cat(paste("\tSample period:", object@sample_period, "\n"))
	cat(paste("\tMaximum sample points:", object@max_sample_points, "\n\n"))
	cat("Use 'plot' to display the results visually\n")
})
