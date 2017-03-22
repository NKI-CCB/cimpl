doCimplAnalysis <- function(
	data,
	scales = (1:15) * 1e4,
	n_iterations = 100,
	sample_period = .1,
	max_sample_points = 2^19,
	chromosomes = unique(data$chr),
	BSgenome,
	system = c('MMTV', 'MuLV', 'SB', 'PB'),
	specificity.pattern,
	lhc.method = c('none', 'exclude'),  # local hopping correction method
	verbose=TRUE,
	cores=1
	) {

	if (any(!chromosomes %in% data$chr)) {
		warning(sprintf('Dropped chromosomes not present in data: %s',
						 paste(chromosomes[!chromosomes %in% data$chr], collapse=', ')))
		chromosomes = chromosomes[chromosomes %in% data$chr]
	}

	# processing arguments
	if (!missing(system)) {
		system <- match.arg(system)
		specificity.pattern <- switch(system,
			MMTV     = '',
			MuLV     = '',
			SB       = 'TA',
			PB       = 'TTAA')
	} else {
		system <- ''
		if ( missing(specificity.pattern) ) {
			specificity.pattern <- ''
			warning('No system and no specificity.pattern defined.')
		}
	}

	lhc.method <- match.arg(lhc.method)

	#print(match.call(expand.dots=TRUE))
	# starting the analysis
	if (verbose) cat('Starting CIMPL analysis...\n')

	cimplObjects <- lapply(chromosomes, function(chr) {
		if (verbose) cat(paste('>>>> chromosome', chr, '<<<<\n'))

		if ( specificity.pattern != '' ) {

			strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

			background <- c(matchPattern(specificity.pattern, BSgenome[[chr]])@ranges@start, matchPattern(strReverse(specificity.pattern), BSgenome[[chr]])@ranges@start)
                        # for R 2.11.1
                        #background <- c(matchPattern(specificity.pattern, BSgenome[[chr]])@start, matchPattern(strReverse(specificity.pattern), BSgenome[[chr]])@start)

		} else {
			background <- NULL
		}

		chr_length <- length(BSgenome[[chr]])

		lapply(scales, function(h) {
			if (verbose) cat(paste('h = ',h , ', ', sep=''))
			n_sample_points <- chr_length / sample_period / h
			if (n_sample_points > max_sample_points) {
#				if (verbose) cat(paste('Warning: undersampling for h =', h , '!! Sample rate set to', max_sample_points, '. '))
				n_sample_points <- max_sample_points
			}

			chr_data <- data[data$chr == chr, ]
			if (lhc.method != 'none') {
				chr_data$hop <- .isHop(chr_data, h * 3)
#				chr_data$hop_weights <- .hopWeights(chr_data, h * 3)
				n_insertions <- sum(!chr_data$hop)
				if (verbose) cat( paste(sum(chr_data$hop), 'hops found...') )
			} else {
				n_insertions <- dim(chr_data)[1]
			}

			if (is.null(background)) {

				#if (verbose) cat('null-densities...')
				#null_densities <- .generateNullDensities(chr_length, h, n_insertions, n_iterations, n_sample_points, verbose=verbose)

				#if (verbose) cat('null-peaks...')
				#null_peaks <- .getPeakDistribution(null_densities)

				if (verbose) cat('null-peaks...')
				null_peaks <- .generatePeakDistribution(chr_length, h, n_insertions, n_iterations, n_sample_points, verbose=verbose, cores=cores)

				if (verbose) cat('cimpl_object...')
				cimplObject <- .getCimplObject(chr_data, chr, null_peaks, h=h, n_sample_points=n_sample_points, verbose=verbose)

				if (verbose) cat('\n')
				return(cimplObject)

			} else {

				#if (verbose) cat('null-densities...')
				#null_densities <- .generateNullDensitiesFromBackground(background, h, n_insertions, n_iterations, n_sample_points, verbose=verbose)

				#if (verbose) cat('null-peaks...')
				#null_peaks <- .getPeakDistribution(null_densities)

				if (verbose) cat('null-peaks...')
				null_peaks <- .generatePeakDistributionFromBackground(background, h, n_insertions, n_iterations, n_sample_points, verbose=verbose, cores=cores)

				if (verbose) cat('at-density...')
				bg_density <- density(background, bw=h, n=n_sample_points)

				if (verbose) cat('cimpl_object...')
				cimplObject <- .getCimplObject(chr_data, chr, null_peaks, bg_density, h, n_sample_points, verbose=verbose)

				if (verbose) cat('\n')
				return(cimplObject)
			}
		})
	})

	new('cimplAnalysis',
		cimplObjects = cimplObjects,
		chromosomes = chromosomes,
		scales = scales,
		system = system,
		specificity.pattern = specificity.pattern,
		lhc.method = lhc.method,
		n_iterations = n_iterations,
		sample_period = sample_period,
		max_sample_points = max_sample_points,
		data.name = deparse(substitute(data)),
		call = match.call()
	)
}

.generatePeakDistribution <- function(chr_length, h, n_insertions, n_iterations, n_sample_points, verbose=TRUE, cores=1) {

	apply_func <- function(i) {
		r_insertions <- sort(floor(runif(n_insertions, min=1, max=chr_length)))
		d <- density(r_insertions, bw=h, n=n_sample_points)  # calculate density
		lms <- .getLocalMaxima(d$y)                          # get peaks
		list(max_pos=d$x[lms$max_pos], max_val=lms$max_val)
	}

	if (cores > 1) {
		suppressMessages(library(parallel))
		local_maxima_list <- mclapply(1:n_iterations, apply_func, mc.cores=cores)
	} else {
		local_maxima_list <- lapply(1:n_iterations, apply_func)
	}

	x <- do.call('c', lapply(local_maxima_list, function(lm) lm$max_pos))
	y <- do.call('c', lapply(local_maxima_list, function(lm) lm$max_val))

	# sort ascending
	srt <- sort.int(x, index.return=TRUE)
	x <- x[srt$ix]
	y <- y[srt$ix]

	return( list( x=x, y=y ) )
}


.generatePeakDistributionFromBackground <- function(background, h, n_insertions, n_iterations, n_sample_points, verbose=TRUE, cores=1) {

	apply_func <- function(i) {
		r_insertions <- sort(sample(background, n_insertions, replace=TRUE))
		d <- density(r_insertions, bw=h, n=n_sample_points)  # calculate density
		lms <- .getLocalMaxima(d$y)                          # get peaks
		list(max_pos=d$x[lms$max_pos], max_val=lms$max_val)
	}

	if (cores > 1) {
		suppressMessages(library(parallel))
		local_maxima_list <- mclapply(1:n_iterations, apply_func, mc.cores=cores)
	} else {
		local_maxima_list <- lapply(1:n_iterations, apply_func)
	}

	x <- do.call('c', lapply(local_maxima_list, function(lm) lm$max_pos))
	y <- do.call('c', lapply(local_maxima_list, function(lm) lm$max_val))

	# sort ascending
	srt <- sort.int(x, index.return=TRUE)
	x <- x[srt$ix]
	y <- y[srt$ix]

	return( list( x=x, y=y ) )
}


.getSampleIDColName <- function(data, possible_sampleID.colNames = c('tumorID', 'tumourID', 'sampleID') ) {
	idx <- which( possible_sampleID.colNames %in% colnames(data) )
	if ( length(idx) > 0 )
		possible_sampleID.colNames[idx[1]]
	else
		stop('data contains no sample id column, use: \'tumorID\', \'tumourID\' or \'sampleID\'')
}


.isHop <- function(data, max.hopping.distance ) {
	hops <- rep(FALSE, dim(data)[1])
	for ( sampleID in unique(data[[.getSampleIDColName(data)]]) ) {
		idx <- which( data[[.getSampleIDColName(data)]] == sampleID )
		# find the hops per sample
		hops.idx <- .getHopIndices(data$location[idx], data$contig_depth[idx], max.hopping.distance)
		if (length(hops.idx) > 0)
			hops[idx[hops.idx]] <- TRUE
	}
	hops
}

.getHopIndices <- function(location, contig_depth, max.hopping.distance) {
    if (is.null(contig_depth)) {
        stop('Missing contig_depth column for hop exclusion')
    }

	# compute distance matrix between locations
	dist_m <- as.matrix(dist(location))

	# which insertions (pairs) are within hopping range
	pairs <- which(dist_m <= max.hopping.distance, arr.ind=TRUE)

	# remove diagonal and duplicates
	pairs <- pairs[pairs[, 1] < pairs[, 2], , drop=FALSE]
	#distances <- dist_m[pairs]

	# for each pair, the insertion with the smallest contig depth is the hopped insertion
	idx1 <- contig_depth[pairs[, 1]] < contig_depth[pairs[, 2]]
	unique ( c ( pairs[idx1, 1], pairs[!idx1, 2] ) )
}


.hopWeights <- function(data, max.hopping.distance ) {
	weights <- rep(1, dim(data)[1])
	for ( sampleID in unique(data[[.getSampleIDColName(data)]]) ) {
		idx <- which( data[[.getSampleIDColName(data)]] == sampleID )
		# get the weightd
		weights[idx] <- .getHopWeights(data$location[idx], data$contig_depth[idx], max.hopping.distance)
	}
	weights
}

.getHopWeights <- function(location, contig_depth, max.hopping.distance) {
	# compute distance matrix between locations
	dist_m <- as.matrix(dist(location))

	# which insertions (pairs) are within hopping range
	pairs <- which(dist_m <= max.hopping.distance, arr.ind=TRUE)

	# remove diagonal and duplicates (upper part)
	pairs <- pairs[pairs[, 1] < pairs[, 2], , drop=FALSE]
	#distances <- dist_m[pairs]

	weights <- rep(1, length(location))
	for (i in 1:dim(pairs)[1]) {
		i1 <- pairs[i, 1]
		i2 <- pairs[i, 2]

		w1 <- weights[i1]
		w2 <- weights[i2]

		cd1 <- contig_depth[i1]
		cd2 <- contig_depth[i2]

		weights[i1] <- cd1/(cd1+cd2) * w1
		weights[i2] <- cd2/(cd1+cd2) * w2
	}
	weights
}

.generateNullDensities <- function(chr_length, h, n_insertions, n_iterations, n_sample_points, verbose=TRUE) {
	# compute the null-distribution
	null_densities <- lapply(1:n_iterations, function(i) {
		if (verbose) cat('.')
		# draw from uniform background
		r_insertions <- sort(floor(runif(n_insertions, min=1, max=chr_length)))

		# calculate density
		density(r_insertions, bw=h, n=n_sample_points)
	})

	return( null_densities )
}


.generateNullDensitiesFromBackground <- function(background, h, n_insertions, n_iterations, n_sample_points, verbose=TRUE) {
	# compute the null-distribution
	null_densities <- lapply(1:n_iterations, function(i) {
		if (verbose) cat('.')
		# draw from TA/AT sites background
		r_insertions <- sort(sample(background, n_insertions, replace=TRUE))

		# calculate density
		density(r_insertions, bw=h, n=n_sample_points)
	})

	return( null_densities )
}

.getPeakDistribution <- function(density_list) {
	local_maxima_list <- lapply(density_list, function(d) {
		lms <- .getLocalMaxima(d$y)
		list(max_pos=d$x[lms$max_pos], max_val=lms$max_val)
	})

	x <- do.call('c', lapply(local_maxima_list, function(lm) lm$max_pos))
	y <- do.call('c', lapply(local_maxima_list, function(lm) lm$max_val))

	# sort ascending
	srt <- sort.int(x, index.return=TRUE)
	x <- x[srt$ix]
	y <- y[srt$ix]

	return( list( x=x, y=y ) )
}

.filterFakePeaks <- function(peaks, h, n) {
	fake.threshold <- 0.1 * dnorm(0, sd=h) / n
	fake.idx <- which(peaks$y < fake.threshold)
	if (length(fake.idx != 0)) {
		peaks$x <- peaks$x[-fake.idx]
		peaks$y <- peaks$y[-fake.idx]
	}
	peaks
}

.getCimplObject <- function(chr_data, chr, null_peaks, bg_density, h, n_sample_points, verbose=TRUE) {
	# calculate density
	if (is.null(chr_data$hop)) {
		insertions <- chr_data$location
	} else {
		insertions <- chr_data$location[!chr_data$hop]
	}
	d <- density(insertions, bw=h, n=n_sample_points)

	# locate peaks
	lms <- .getLocalMaxima(d$y)
	peaks <- list( x=d$x[lms$max_pos], y=lms$max_val )

	# Due to numerical aberrations there are fake peaks in regions where the kse is almost 0.
	null_peaks <- .filterFakePeaks(null_peaks, h, length(insertions))
	peaks <- .filterFakePeaks(peaks, h, length(insertions))

	if (!missing(bg_density)) {
		# calculate peak priors
		peaks$priors <- approx(bg_density, xout=peaks$x)$y

		# calculate null-peak priors
		null_peaks$priors <- approx(bg_density, xout=null_peaks$x)$y

		# estimate prior-ph joint pdf and make it a joint cdf, conditioned on the prior
		if (length(null_peaks$y) > 50e3) { # use binned-kde to reduce calculation time and memory use
			if (verbose) cat('(bkde)...')
			bw2d <- c(bandwidth.nrd(null_peaks$priors), bandwidth.nrd(null_peaks$y))
			est <- bkde2D(cbind(null_peaks$priors, null_peaks$y ), gridsize=c(100,100), bandwidth=bw2d)
			names(est) <- c('x', 'y', 'z')
		} else {
			est <- kde2d(null_peaks$priors, null_peaks$y, n=100)        # calculate prior-ph density
		}
		est$z <- t(apply( est$z, 1, function(x) cumsum(x)/sum(x) )) # make it a CDF-surface

		# calculate p-vals
		peaks$p_vals <- 1 - .biapprox( est, cbind(peaks$priors, peaks$y) )    # do a bi-linear approximation for the observed ph's - priors.
		peaks$p_vals[is.na(peaks$p_vals)] <- as.numeric( peaks$y[is.na(peaks$p_vals)] < max(est$y) )  # fix peaks that fall outside the density boundaries

		# somehow, numerical aberrations occur resulting in negative p-vals
		peaks$p_vals[peaks$p_vals < 0] <- 0
	} else {
		# calculate p-vals (based on the empirical CDF of peak heights)
		peaks$p_vals <- sapply(peaks$y, function(height) {
			sum(null_peaks$y >= height) / length(null_peaks$y)
		})
	}

	n_peaks <- length(peaks$x)

	if (!missing(bg_density)) {
		cimplObject <- new('cimplObject',
			data = chr_data,
			kse = list(x=d$x, y=d$y),
			peaks = peaks,
			null_peaks = null_peaks,
			null_cdf = est,
			bg_density = list(x=bg_density$x, y=bg_density$y),
			n_peaks = n_peaks,
			scale = h,
			chromosome = chr
		)
	} else {
		cimplObject <- new('cimplObject',
			data = chr_data,
			kse = list(x=d$x, y=d$y),
			peaks = peaks,
			null_peaks = null_peaks,
			n_peaks = n_peaks,
			scale = h,
			chromosome = chr
		)
	}

	return(cimplObject)
}

.getLocalMaxima <- function(x) {
	####### a port of Jeroens matlab code #######
	diff_vec            <- x[-1] - x[-length(x)]
	# Make a tri-state vector indicating the transitions
	tristate <- sign(diff_vec)
	# Remove the flat parts
	non_flat_index      <- which(tristate != 0)
	bistate             <- tristate[tristate != 0]
	# Find the 1-->-1 transitions
	trans_vec           <- bistate[-length(bistate)] - bistate[-1]
	trans_index         <- which(trans_vec == 2);
	if (length(trans_index) == 0) {
		stop('No local maxima!')
	}
	max_pos             <- non_flat_index[trans_index+1];
	max_val             <- x[max_pos];
	#############################################

	return (list(max_pos=max_pos, max_val=max_val))
}

.biapprox <- function (obj, loc) {
	# obj is a surface object like the list for contour or image.
	# loc is a matrix of (x, y) locations
	x <- obj$x
	y <- obj$y
	x.new <- loc[,1]
	y.new <- loc[,2]
	z <- obj$z

	ind.x <- findInterval(x.new, x, all.inside=T)
	ind.y <- findInterval(y.new, y, all.inside=T)

	ex <- (x.new - x[ind.x]) / (x[ind.x + 1] - x[ind.x])
	ey <- (y.new - y[ind.y]) / (y[ind.y + 1] - y[ind.y])

	# set weights for out-of-bounds locations to NA
	ex[ex < 0 | ex > 1] <- NA
	ey[ey < 0 | ey > 1] <- NA

#	return(
#		z[cbind(ind.y,     ind.x)]     * (1 - ex) * (1 - ey) +  # upper left
#		z[cbind(ind.y + 1, ind.x)]     * (1 - ex) * ey       +  # lower left
#		z[cbind(ind.y + 1, ind.x + 1)] * ex       * ey       +  # lower right
#		z[cbind(ind.y,     ind.x + 1)] * ex       * (1 - ey)    # upper right
#	)
	return(
		z[cbind(ind.x,     ind.y)]     * (1 - ex) * (1 - ey) +  # upper left
		z[cbind(ind.x + 1, ind.y)]     * ex       * (1 - ey) +  # lower left
		z[cbind(ind.x + 1, ind.y + 1)] * ex       * ey       +  # lower right
		z[cbind(ind.x,     ind.y + 1)] * (1 - ex) * ey          # upper right
	)
}

.sortDataFrame <- function(x, key, ...) {
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}

.getTumorDensities <- function(cimplObject, cumulative=FALSE) {

	ids <- unique(cimplObject@data[[.getSampleIDColName(cimplObject@data)]])
	tds <- matrix(0, length(cimplObject@kse$x), length(ids))

	for(i in 1:length(ids)) {
		if (is.null(cimplObject@data$hop)) {
			insertions <- cimplObject@data$location[cimplObject@data[[.getSampleIDColName(cimplObject@data)]] == ids[i]]
			d <- density(insertions, cimplObject@scale, n=length(cimplObject@kse$x))
			tds[,i] <- approx(d, xout=cimplObject@kse$x)$y * length(insertions) / dim(cimplObject@data)[1]
		} else {
			insertions <- cimplObject@data$location[cimplObject@data[[.getSampleIDColName(cimplObject@data)]] == ids[i] & !cimplObject@data$hop ]
			d <- density(insertions, cimplObject@scale, n=length(cimplObject@kse$x))
			tds[,i] <- approx(d, xout=cimplObject@kse$x)$y * length(insertions) / sum(!cimplObject@data$hop)
		}
	}
	colnames(tds) <- ids

	if (cumulative) {
		tds[is.na(tds)] <- 0
		for(i in 1:dim(tds)[1]) {
			tds[i,] <- cumsum(tds[i,])
		}
	}

	tds
}

getEnsemblGenes <- function(cimplAnalysis, geneIdentifiers=c('ensembl_gene_id', 'external_gene_id'), mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl")) {
	genes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'), filters='chromosome_name', values=substring(cimplAnalysis@chromosomes, 4), mart=mart)

	# , 'unigene'
	# filter out unknown
	if (!all(i <- geneIdentifiers %in% listAttributes(mart)[,1])) {
		stop(paste(geneIdentifiers[!i], ' is/are no valid biomaRt attributes for this mart.'))
	}

	if (length(geneIdentifiers) > 0) {
		if ( !(length(geneIdentifiers) == 1 & geneIdentifiers[1] == 'ensembl_gene_id' )) {
			annot <- getBM(attributes=unique(c('ensembl_gene_id', geneIdentifiers)), filters='ensembl_gene_id', values=genes$ensembl_gene_id, mart=mart)

			extra_annot <- sapply(genes$ensembl_gene_id, function(ens) {
				idx <- annot$ensembl_gene_id == ens
				sapply(2:ncol(annot), function(col) {
					vals <- annot[idx, col]
					vals[vals==''] <- NA
					paste(unique(na.omit(vals)), collapse='|')
				})
			})

			if (!is.null(dim(extra_annot)))
				extra_annot <- t(extra_annot)

			extra_annot <- data.frame(extra_annot, stringsAsFactors=FALSE)
			colnames(extra_annot) <- colnames(annot)[-1]
			genes <- cbind(genes, extra_annot)
		}
	}
	if (! 'ensembl_gene_id' %in% geneIdentifiers) {
		genes$ensembl_gene_id <- NULL
	}

	attr(genes, 'geneIdentifiers') <- geneIdentifiers

	genes
}

# @deprecated: use getCISs
.getPeaks <- function(cimplAnalysis, alpha=0.05, mul.test=TRUE, mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl"), genes = getEnsemblGenes(cimplAnalysis, mart), order.by='p_value') {

	df <- do.call('rbind', lapply(cimplAnalysis@chromosomes, function(chr) {
		chr.idx <- which(cimplAnalysis@chromosomes == chr)

		chr_genes <- genes[genes$chromosome_name == substring(chr, 4), ,drop=FALSE]

		do.call('rbind', lapply(cimplAnalysis@scales, function(kw) {

			kw.idx <- which(cimplAnalysis@scales == kw)
			if (mul.test) {
				n_tests <- sum(sapply(cimplAnalysis@cimplObjects, function(x) { x[[kw.idx]]@n_peaks }))
				print(n_tests)
			} else {
				n_tests <- 1
			}

			cimpl <- cimplAnalysis@cimplObjects[[chr.idx]][[kw.idx]]
			associated_genes <- sapply(round(cimpl@peaks$x), function(loc) {
				d <- 0
				idx <- chr_genes$start_position - d < loc & chr_genes$end_position + d > loc
				if ( sum(idx) == 0 ) {
					c('', '')
				} else if ( sum(idx) == 1 ) {
					as.character(chr_genes[which(idx) ,c(1, 2)])
				} else {
					# this peak has hit multiple genes!
					c(paste(chr_genes[which(idx) , 1], collapse='|'),
					  paste(chr_genes[which(idx) , 2], collapse='|')  )
				}
			})

			data.frame(
				location         = round(cimpl@peaks$x),
				chromosome       = rep(chr, cimpl@n_peaks),
				ensembl_gene_id  = associated_genes[1,],
				external_gene_id = associated_genes[2,],
				p_value          = cimpl@peaks$p_vals,
				significant      = cimpl@peaks$p_vals < alpha / n_tests,
				scale            = as.numeric(rep(kw, cimpl@n_peaks)),
				stringsAsFactors = FALSE
			)
		}))
	}))

	.sortDataFrame(df, order.by, decreasing=FALSE)
}

.CISThresholdLine <- function(cimplObject, alpha) {

	if ( length(cimplObject@null_cdf) == 0 ) {
		list(x=cimplObject@kse$x, y=rep( quantile(cimplObject@null_peaks$y, prob=1-alpha) , length(cimplObject@kse$x)))
	} else {
		null_cdf <- cimplObject@null_cdf
		cuts <- apply(null_cdf$z, 1, function(cdf) null_cdf$y[which(cdf >= 1-alpha)[1]])
		th_line <- approx(x=null_cdf$x, y=cuts, xout=cimplObject@bg_density$y)
		th_line$x <- cimplObject@bg_density$x
		th_line
	}
}

.getRegionsAboveThreshold <- function(cimplObject, alpha) {
	if ( length(cimplObject@null_cdf) == 0 ) {
		co <- quantile(cimplObject@null_peaks$y, prob=1-alpha)
		segments <- .getSegments(cimplObject@kse$y > co)
	} else {
		priors <- approx(cimplObject@bg_density, xout=cimplObject@kse$x)$y
		kse.p_vals <- 1 - .biapprox(cimplObject@null_cdf, cbind(priors, cimplObject@kse$y) )
		kse.p_vals[is.na(kse.p_vals)] <- as.numeric( cimplObject@kse$y[is.na(kse.p_vals)] < max(cimplObject@null_cdf$y) )  # fix value that fall outside the density boundaries
		segments <- .getSegments(kse.p_vals < alpha)
	}
	list(start_pos=round(cimplObject@kse$x[segments$start_pos]), end_pos=round(cimplObject@kse$x[segments$end_pos]))
}

.getSegments <- function(x) {
	# 'x' is a logical vector
#	above_vec           <- x > co
	trans_vec           <- x[-length(x)] - x[-1]
	# at 0 --> 1 a region start
	# at 1 --> 0 a the region ends
	start_pos <- which(trans_vec == -1)
	end_pos <- which(trans_vec == 1) - 1
	list(start_pos=start_pos, end_pos=end_pos)
}

## Genes can be supplied from Ensembl as:
## 	 mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl")
##   genes = getEnsemblGenes(cimplAnalysis, mart)
getCISs <- function(cimplAnalysis, alpha=0.05, chromosomes=cimplAnalysis@chromosomes, scales=cimplAnalysis@scales, mul.test=TRUE, genes=NULL, order.by=c('p_value', 'n_insertions'), decreasing=c(FALSE, TRUE)) {

	df <- do.call('rbind', lapply(chromosomes, function(chr) {
		chr.idx <- which(cimplAnalysis@chromosomes == chr)
		do.call('rbind', lapply(scales, function(kw) {

			kw.idx <- which(cimplAnalysis@scales == kw)
			if (mul.test) {
				n_tests <- sum(sapply(cimplAnalysis@cimplObjects, function(x) { x[[kw.idx]]@n_peaks }))
			} else {
				n_tests <- 1
			}

			cimplObject <- cimplAnalysis@cimplObjects[[chr.idx]][[kw.idx]]

			regions <- .getRegionsAboveThreshold(cimplObject, alpha / n_tests)
			n_regions <- length(regions$start_pos)

			if (n_regions == 0) {
				return(NULL)
			} else if (n_regions == 1) {
				annDf <- .annotateRegion(chr, kw, regions$start_pos, regions$end_pos, cimplObject, genes)
			} else {
				annDf <- do.call('rbind', lapply(1:n_regions, function(i) {
					.annotateRegion(chr, kw, regions$start_pos[i], regions$end_pos[i], cimplObject, genes)
				}))
			}

			return(annDf)
		}))
	}))

	if (is.null(df)) {
		# no CISs found!
		NULL
	} else if(dim(df)[1] == 0) {
		# no CISs found!
		NULL
	} else {

		# make some nicely formatted CIS ids and set them as rownames
		rownames(df) <- .makeCISID(df$chromosome, df$peak_location, df$scale)

		# sort
		for (i in 1:length(order.by)) {
			df <- .sortDataFrame(df, order.by[i], decreasing=decreasing[i])
		}
		df
	}
}

.makeCISID <- function(chr, loc, scale) {
	paste('CIS', substring(chr, 4), ':', loc, '_', .formatScales(scale), sep='')
}

.annotateRegion <- function(chr, scale, bpStart, bpEnd, cimplObject, genes=NULL) {
	if (!is.null(genes))
		chr_genes <- genes[genes$chromosome_name == substring(chr, 4), ,drop=FALSE]
	else
		chr_genes <- NULL

	peak.idx <- which(round(cimplObject@peaks$x) >= bpStart & round(cimplObject@peaks$x) <= bpEnd)
	n_peaks <- length(peak.idx)

	if (n_peaks == 0) {
		# CIS contains no peaks
		return(NULL)
	} else if (n_peaks == 1) {
		# CIS contains 1 peak
		.annotatePeak(chr, scale, bpStart, bpEnd, peak.idx, cimplObject, chr_genes)
	} else {
		# CIS contains 2 or more peaks
		do.call('rbind', lapply(1:n_peaks, function(i) {
			.annotatePeak(chr, scale, bpStart, bpEnd, peak.idx[i], cimplObject, chr_genes)
		}))
	}
}

.annotatePeak <- function(chr, scale, cisStart, cisEnd, peak.idx, cimplObject, chr_genes=NULL) {
	peakLoc <- round(cimplObject@peaks$x[peak.idx])


	ann <- data.frame(row.names='')

	ann$chromosome <- chr
	ann$peak_location <- peakLoc
	ann$peak_height <- cimplObject@peaks$y[peak.idx] * dim(cimplObject@data)[1] / dnorm(0, 0, sd=scale)
	ann$start <- cisStart
	ann$end <- cisEnd
	ann$width <- cisEnd - cisStart + 1

	snappedLocs <- .snapToPeaks(cimplObject@data$location, round(cimplObject@peaks$x))
	ann$n_insertions <- sum(snappedLocs >= cisStart & snappedLocs <= cisEnd)

	ann$p_value <- cimplObject@peaks$p_vals[peak.idx]
	ann$scale <- scale

	if (!is.null(chr_genes)) {
		ags <- .associateGenes(chr, cisStart, cisEnd, peakLoc, chr_genes)
		for(id in attr(chr_genes, 'geneIdentifiers')) {

			a <- chr_genes[ags$associated, id]
			a <- a[a!='']

			o <- chr_genes[ags$other, id]
			o <- o[o!='']

			ann[paste('associated_', id, sep='')] <- paste(a, collapse='|')
			ann[paste('other_', id, sep='')] <- paste(o, collapse='|')
		}
	}

	return(ann)
}

.associateGenes <- function(chr, cisStart, cisEnd, peakLoc, chr_genes) {
	# 1. find all genes which overlap the CIS region +/- 100KB
	margin <- 100e3
	# the genes within the region + margin
	gene.idx <- which(chr_genes$start_position < (cisEnd + margin) & chr_genes$end_position > (cisStart - margin))

	if (length(gene.idx) == 0) {
		list(
			associated = numeric(0),
			other      = numeric(0)
		)
	} else if (length(gene.idx) == 1) {
		list(
			associated = gene.idx,
			other      = numeric(0)
		)
	} else {
		# 2. calculate the distances
		dists <- pmin(abs(chr_genes$start_position[gene.idx] - peakLoc), abs(chr_genes$end_position[gene.idx] - peakLoc))

		# genes which contain the peak get distance 0
		dists[chr_genes$start_position[gene.idx] <= peakLoc & chr_genes$end_position[gene.idx] >= peakLoc] <- 0

		# increase the dists of non-MGI (automatic) so that they end up at the end of the associated gene list
#		max.dist <- max(dists)
#		automatic.idx <- which(chr_genes$external_gene_db[gene.idx] == 'MGI (automatic)')
#		curated.idx <- which(chr_genes$external_gene_db[gene.idx] == 'MGI (curated)')
#		dists[curated.idx] <- dists[curated.idx] - 20 * max.dist
#		dists[automatic.idx] <- dists[automatic.idx] - 10 * max.dist

		gene.order <- order(dists, decreasing=FALSE)

		min.ties <- length(which(dists == min(dists)))

		list(
			associated = gene.idx[gene.order][1:min.ties],
			other      = gene.idx[gene.order][-(1:min.ties)]
		)
#		list(
#				associated = chr_genes$external_gene_id[gene.idx[gene.order]][1:min.ties],
#				other      = chr_genes$external_gene_id[gene.idx[gene.order]][-(1:min.ties)]
#		)
	}
}

.formatScales <- function(scales) {
	str <- as.character(scales)
	k.idx <- scales %% 1e3 == 0
	m.idx <- scales %% 1e6 == 0
	str[k.idx] <- paste(scales[k.idx] / 1e3, 'k', sep='')
	str[m.idx] <- paste(scales[k.idx] / 1e6, 'm', sep='')
	str
}

.snapToPeaks <- function(locations, peaks) {
	sapply(locations, function(loc) round(peaks[which.min(abs(peaks - loc ))]))
}

getCISMatrix <- function(cimplAnalysis, ciss) {
	df <- do.call('rbind', lapply(cimplAnalysis@chromosomes, function(chr) {
		chr.idx <- which(cimplAnalysis@chromosomes == chr)

		chr_data <- cimplAnalysis@cimplObjects[[chr.idx]][[1]]@data

		cisids <- do.call('cbind', lapply(cimplAnalysis@scales, function(kw) {

			kw.idx <- which(cimplAnalysis@scales == kw)
			cimplObject <- cimplAnalysis@cimplObjects[[chr.idx]][[kw.idx]]

			# snap insertions to peaks (see http://bioinformatics.nki.nl/forum/viewtopic.php?f=4&t=19)
			snappedLocs <- .snapToPeaks(chr_data$location, cimplObject@peaks$x)

			insertion2cis <- rep('', dim(chr_data)[1])

			ciss.idx <- which(ciss$chromosome == chr & ciss$scale == kw)
			for (i in ciss.idx) {
#				insertion2cis[snappedLocs >= ciss$start[i] & snappedLocs <= ciss$end[i]] <- rownames(ciss)[i]
#				insertion2cis[snappedLocs == ciss$peak_location[i]] <- rownames(ciss)[i]
				locs.idx <- snappedLocs >= ciss$start[i] & snappedLocs <= ciss$end[i]

				insertion2cis[locs.idx] <- paste(insertion2cis[locs.idx], rownames(ciss)[i], sep='|')
			}

			substring(insertion2cis, 2)
		}))
		colnames(cisids) <- cimplAnalysis@scales
		data.frame(chr_data, cisids, stringsAsFactors=FALSE)
	}))
}

getInsertions <- function(cimplAnalysis, chr, scale, bpLim) {
	chr.idx <- which(cimplAnalysis@chromosomes == chr)
	scale.idx <- which(cimplAnalysis@scales == scale)
	chr_data <- cimplAnalysis@cimplObjects[[chr.idx]][[scale.idx]]@data
	#chr_data[chr_data$location >= bpLim[1] & chr_data$location <= bpLim[2], ]

	# snap insertions to peaks (see http://bioinformatics.nki.nl/forum/viewtopic.php?f=4&t=19)
	snappedLocs <- .snapToPeaks(chr_data$location, cimplAnalysis@cimplObjects[[chr.idx]][[scale.idx]]@peaks$x)
	chr_data[snappedLocs >= bpLim[1] & snappedLocs <= bpLim[2], ]
}

.nToString <- function(n, nDigits, prefix='', postfix='') {
	out <- NULL
	sapply(n, function(n) {
		while(n > 0) {
			out <- c(n%%10 , out)
			n <- n%/%10
		}

		if (length(out) < nDigits) {
			out <- c( rep(0, nDigits - length(out)) , out)
		}

		paste(prefix, paste(out, collapse=''), postfix, sep='')
	})
}

