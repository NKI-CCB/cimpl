setMethod("plot", signature(x="cimplAnalysis", y="missing"), function(x, y, type=c('kse', 'null.cdf', 'scale.space'), chr=x@chromosomes[1], scale=x@scales[1], alpha=0.05, mul.test=TRUE, bpLim, plot.tumor.densities=FALSE, interactive=TRUE, ...) {
	chr.idx <- which(x@chromosomes == chr)
	kw.idx <- which(x@scales == scale)
	co <- x@cimplObjects[[chr.idx]][[kw.idx]]
	
	type <- match.arg(type)
	
	if (mul.test) {
		n_tests <- sum(sapply(x@cimplObjects, function(y) { y[[kw.idx]]@n_peaks }))  # total number of peaks for the selected scale
	} else {
		n_tests <- 1
	}

	if (missing(bpLim)) {
		bpLim <- range(co@kse$x)
	}
	
	if (type == 'scale.space') {
		
		plotfunc <- function(lim) {
			if (length(x@scales) > 1) {
				rh <- min(x@scales[-1] - x@scales[-length(x@scales)])  # rectangle height
			} else {
				rh <- 1
			}
			
			plot(1, type='n', xlim=lim$x, ylim=range(x@scales) + c(-rh/2, rh/2), xlab='Genomic position', ylab='Scale', main='Scale space')
			
			regions <- lapply(x@cimplObjects[[chr.idx]], function(co) {
				if (mul.test) {
					alpha <- alpha / n_tests
				}
				regions <- .getRegionsAboveThreshold(co, alpha)
	
				vis.idx <- regions$start_pos <= lim$x[2] & regions$end_pos >= lim$x[1]
				if (any(vis.idx)) {
					rect(xleft=regions$start_pos[vis.idx], ybottom=co@scale-rh/2, xright=regions$end_pos[vis.idx], ytop=co@scale+rh/2, col='green', border=NA)
				}
			})
		}
		
		extent <- list(x=bpLim, y=range(x@scales))
		
		if (interactive) {
			invisible(.zoom(plotfunc, extent=extent))
		} else {
			plotfunc(lim=extent)
		}
		
	} else {
		plot(co, type=type, plot.tumor.densities=plot.tumor.densities, alpha=alpha, n_tests=n_tests, bpLim=bpLim, interactive=interactive, ...)
	}
})

setMethod("plot", signature(x="cimplObject", y="missing"), function(x, y, type=c('kse', 'null.cdf'), alpha=0.05, n_tests=1, bpLim, plot.tumor.densities=FALSE, interactive=TRUE, ...) {
	type <- match.arg(type)
	info <- paste('chromosome:', x@chromosome, ', scale:', x@scale, ', # peaks:', x@n_peaks)
	
	if (type == 'kse') {
		
		plotfunc <- function(lim) {
			ylim <- c(0, max(x@kse$y[x@kse$x >= lim$x[1]-x@scale & x@kse$x <= lim$x[2]+x@scale]))
			
			# a zoom factor based on the kernel width
			zf <- (lim$x[2] - lim$x[1]) / x@scale

			plot(x@kse, xlim=lim$x, ylim=ylim, type='l', xlab='genome', ylab='density', col='blue', main=paste(info, '\n', x@chromosome, ':', floor(lim$x[1]), '-', ceiling(lim$x[2]), sep=''), ...)
			if (plot.tumor.densities & zf < 100) {
				tumor.densities <- .getTumorDensities(x, cumulative=TRUE)
				for (j in dim(tumor.densities)[2]:1) {
					polygon(x@kse$x, tumor.densities[,j], col=j, border=1)
				}
				if (dim(tumor.densities)[2] < 15)
					legend('topright', legend=colnames(tumor.densities), col=1:dim(tumor.densities)[2], lty=1)
			}
			lines(x@bg_density, col='green')
			points(x@peaks, col='red')
			
			lines(.CISThresholdLine(x, alpha=alpha), col='red')
			if (n_tests > 1) {
				lines(.CISThresholdLine(x, alpha=alpha / n_tests), col='red', lty=2)
			}
			
			points(x=x@data$location, y=rep(0,sum(dim(x@data)[1])), col=match(x@data[[.getSampleIDColName(x@data)]], unique(x@data[[.getSampleIDColName(x@data)]])), pch='|')
			
			if (!is.null(x@data$hop)) {
				# 'X' indicates a hop
				points(x=x@data$location[x@data$hop], y=rep(0,sum(x@data$hop)), pch='x')
			}
		}
		
		if (missing(bpLim)) {
			bpLim <- range(x@kse$x)
		}
		extent <- list(x=bpLim, y=c(0, max(x@kse$y)))
		
		if (interactive) {
			invisible(.zoom(plotfunc, extent=extent))
		} else {
			plotfunc(lim=extent)
		}
	} else if (type == 'null.cdf') {
		persp( x@null_cdf, ticktype="detailed", theta=60, phi=30, expand=0.5, shade=0.5, col="cyan", ltheta=-30, xlab='bg density', ylab='peak height', zlab='cumulative density', main=info, ...)
	}
})

.zoom <- function(plotfunc, extent){
	# (C) 2007 GPL by Huidae Cho <http://geni.ath.cx>
	# Simple R function for interactive zooming/panning
	#
	# Example:
	# data <- runif(100)*10
	# extent <- list(x=c(1, 100), y=c(0, 10))
	# plotfunc <- function(lim){
	# 	plot(data, xlim=lim$x, ylim=lim$y)
	# 	abline(mean(data), 0, col="red")
	# }
	# zoom(plotfunc, extent)
 
	print("Zoom in:     Click two corners", quote=FALSE)
	print("Zoom out:    Click above plot", quote=FALSE)
	print("Prev extent: Click left of plot", quote=FALSE)
	print("Next extent: Click right of plot", quote=FALSE)
	print("Full extent: Click below plot", quote=FALSE)
	print("Pan:         Double click", quote=FALSE)
	print("Quit:        Right button", quote=FALSE)
 
	lim <- extent
	lim.stack <- c(lim$x, lim$y)
	lim.depth <- 1
	lim.cur <- 1
 
	repeat{
		plotfunc(lim)
 
		l <- locator(1)
		if(is.null(l))
			break
		ext <- par()$usr
		if(l$x < ext[1] || l$x > ext[2]){
			cur <- lim.cur
			lim.cur <- if(l$x < ext[1]) max(lim.cur-1, 1)
				else min(lim.cur+1, lim.depth)
			if(lim.cur != cur)
				lim <- list(x=lim.stack[lim.cur, 1:2],
					y=lim.stack[lim.cur, 3:4])
			next
		}
		if(l$y < ext[3])
			lim <- extent
		else
		if(l$y > ext[4]){
			cx <- (lim$x[1] + lim$x[2]) / 2
			cy <- (lim$y[1] + lim$y[2]) / 2
			w <- lim$x[2] - lim$x[1]
			h <- lim$y[2] - lim$y[1]
			lim <- list(x=c(cx-w, cx+w), y=c(cy-h, cy+h))
		}else{
			l2 <- locator(1)
			if(is.null(l2))
				break
			if(sum(l$x == l2$x) || sum(l$y == l2$y)){
				w <- lim$x[2] - lim$x[1]
				h <- lim$y[2] - lim$y[1]
				lim <- list(x=c(l2$x-w/2, l2$x+w/2),
					y=c(l2$y-h/2, l2$y+h/2))
			}else
				lim <- list(x=sort(c(l$x, l2$x)),
					y=sort(c(l$y, l2$y)))
		}
		if(lim.cur < lim.depth){
			lim.stack <- lim.stack[-((lim.cur+1):lim.depth),]
			lim.depth <- lim.cur
		}
		lim.stack <- rbind(lim.stack, c(lim$x, lim$y))
		lim.depth <- lim.depth + 1
		lim.cur <- lim.cur + 1
	}
 
	lim
}

