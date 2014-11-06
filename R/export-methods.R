
export.wig <- function(cimplAnalysis, chromosomes=cimplAnalysis@chromosomes, scales=cimplAnalysis@scales, file) {
	wig.file <- file(file, 'w')
	
	for (chr in chromosomes) {
		chr.idx <- which(cimplAnalysis@chromosomes == chr)
		for (scale in cimplAnalysis@scales) {
			scale.idx <- which(cimplAnalysis@scales == scale)
			cimplObject <- cimplAnalysis@cimplObjects[[chr.idx]][[scale.idx]]
			
			# export kse
			kse <- cimplObject@kse
			cat(paste('track type=wiggle_0 name="KSE_', scale, '" description="KSE_', scale, '" visibility=full autoScale=on color=0,0,255', '\n', sep=''), file=wig.file)
			cat(paste('variableStep chrom=', chr, '\n', sep=''), file=wig.file)
			cat(paste(format(round(kse$x[kse$x>0]), scientific=FALSE) , kse$y[kse$x>0], sep='\t'), sep='\n', file=wig.file)
			
			# export background density
			bg <- cimplObject@bg_density
			cat(paste('track type=wiggle_0 name="BG_', scale, '" description="BG_', scale, '" visibility=full autoScale=on color=0,255,0', '\n', sep=''), file=wig.file)
			cat(paste('variableStep chrom=', chr, '\n', sep=''), file=wig.file)
			cat(paste(format(round(bg$x[bg$x>0]), scientific=FALSE), bg$y[bg$x>0], sep='\t'), sep='\n', file=wig.file)
		}
	}
	close(wig.file)
}

export.bed <- function(ciss, file) {
	#re-jig the ciss data frame
	#bed <- data.frame(chr=ciss$chr, start=ciss$start, end=ciss$end, CIS=paste(rownames(ciss), ciss$chr, ciss$peak_location,sep='_'), score=ciss$n_insertions, strand='+')
	bed <- data.frame(chr=ciss$chr, start=ciss$start, end=ciss$end, CIS=rownames(ciss), score=ciss$n_insertions, strand='+')
	write.table(bed, file=file,sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

export.html <- function(cimplAnalysis, chromosomes=cimplAnalysis@chromosomes, scales=cimplAnalysis@scales, genes, alpha=0.05, mul.test=TRUE, plot.tumor.densities=FALSE, order.by=c('scale', 'n_insertions', 'p_value'), decreasing=c(FALSE, TRUE, FALSE), dir='.', verbose=TRUE) {
	# find the CISs
	ciss <- getCISs(cimplAnalysis, chromosomes=chromosomes, scales=scales, alpha=alpha, mul.test=mul.test, genes=genes, order.by=order.by, decreasing=decreasing)
	if (is.null(ciss) | dim(ciss)[1] == 0) {
		warning('no CISs found, nothing exported!')
		return()
	}
	
	ciss_col_formats <- c('s', 's', 'd', 'g', 'd', 'd', 'd', 'd', 'g', 'd')
	ciss_col_formats <- c( ciss_col_formats, rep('s', length(attr(genes, 'geneIdentifiers')) * 2) )
	#ciss_col_formats <- c('d', 's', 's', 'd', 'g', 'd', 'd', 'd', 'd', 's', 'g', 'd')
	
	image_dir_name <- 'images'
	details_dir_name <- 'details_files'
	image_dir <- paste(dir, image_dir_name, sep='/')
	details_dir <- paste(dir, details_dir_name, sep='/')

	# make dirs if nessecary
	dir.create(dir)
	dir.create(image_dir)
	dir.create(details_dir)

	# TODO:make chromosome plots
	
	
	# make a details page for each CIS
	figures <- sapply( 1:dim(ciss)[1], function(i) {
		if (verbose) cat(paste('Exporting CIS ', i, '/', dim(ciss)[1], '\r', sep=''))
		# make plots
		plotLims <- c( ciss$start[i] - 3 * ciss$scale[i], ciss$end[i] + 3 * ciss$scale[i])
		
		# kse
		kse_file_name <- paste('cis', i, '.png', sep='')
		png(file=paste(image_dir, kse_file_name, sep='/') )
		plot(cimplAnalysis, type='kse', alpha=alpha, mul.test=mul.test, chr=ciss$chromosome[i], scale=ciss$scale[i], bpLim=plotLims, plot.tumor.densities=plot.tumor.densities, interactive=FALSE)
		dev.off()
		
		# scale space
		scale_space_file_name <- paste('scale_space', i, '.png', sep='')
		png(file=paste(image_dir, scale_space_file_name, sep='/') )
		plot(cimplAnalysis, type='scale.space', alpha=alpha, mul.test=mul.test, chr=ciss$chromosome[i], scale=ciss$scale[i], bpLim=plotLims, interactive=FALSE)
		dev.off()
		
		# write the details file
		details_file_name <- paste('details', i, '.html', sep='')
		
		print(xtable(ciss[i,], display=ciss_col_formats), append=FALSE, type='html', include.rownames=TRUE, sanitize.text.function=function(x) x, file=paste(details_dir, details_file_name, sep='/'))
	
		
		details_file <- file( paste(details_dir, details_file_name, sep='/'), 'a')
		cat(paste('</br><a href=\'', .getUrl(chr=ciss$chromosome[i], bpLim=c(ciss$start[i], ciss$end[i])), '\'>View</a> in Ensembl Genome Browser</br></br>\n', sep=''), file=details_file)
		cat(paste('<img src=\'../', image_dir_name, '/', kse_file_name, '\' /></br>\n', sep=''), file=details_file)
		cat(paste('\tBlue line: kernel smoothed estimate</br>\n'), file=details_file)
		cat(paste('\tGreen line: \'specificity pattern\'-density</br>\n'), file=details_file)
		cat(paste('\tRed line: significance threshold</br>\n'), file=details_file)
		cat(paste('\tDotted red line: significance threshold after multiple testing correction</br>\n'), file=details_file)
		cat(paste('\tRed circle: peaks location</br>\n'), file=details_file)
		cat(paste('\tColored vertical lines: insertion locations color-coded by tumor or sample</br></br>\n\n'), file=details_file)
		
		cat(paste('<img src=\'../', image_dir_name, '/', scale_space_file_name, '\' /></br>\n', sep=''), file=details_file)
		
		
		cat(paste('Insertions located in the CIS:</br>\n'), file=details_file)
		close(details_file)
		
		# append the insertion details
		linked_insertions <- getInsertions(cimplAnalysis, chr=ciss$chromosome[i], scale=ciss$scale[i], bpLim=c(ciss$start[i], ciss$end[i]))
		xt_details <- xtable(linked_insertions)
		
		if (dim(xt_details)[1] != 0)
			print(xt_details, append=TRUE, type='html', include.rownames=FALSE, sanitize.text.function=function(x) x, file=paste(details_dir, details_file_name, sep='/'))
		
		paste(details_dir_name, details_file_name, sep='/')
	})
	if (verbose) cat('\n')

	# add a link to figures (a details column)
	ciss$details <- sapply(figures, function(f) paste('<a href=\'', f, '\'>details</a>', sep=''))
	
	index_file_name <- paste(dir, 'index.html', sep='/')
	# write index.html with the table.
	index_file <- file(index_file_name, 'w')
	#cat(paste('cimpl_', sessionInfo(package='cimpl')[[4]]$cimpl$Version, ' - ', date(), sep=''), file=index_file)
	cat(paste(date()), file=index_file)
	close(index_file)
	
	#xt_index <- xtable(ciss, display=c(ciss_col_formats, 's'), caption='All CISs (for CISs per scale: scroll down)')
	#print(xt_index, append=FALSE, type='html', include.rownames=FALSE, caption.placement='top', sanitize.text.function=function(x) x, file=index_file_name)

	
	# a table by scale:
	for (s in scales) {
		idx <- ciss$scale == s
		if (any(idx)) {
			xt_scale <- xtable(ciss[ciss$scale == s, ], display=c(ciss_col_formats, 's'), caption=paste('CISs for scale', s), label=paste('s', s, sep=''))
			print(xt_scale, append=TRUE, type='html', include.rownames=TRUE, caption.placement='top', sanitize.text.function=function(x) x, file=index_file_name)
		}
	}
	
	# append anlysis info
	index_file <- file(index_file_name, 'a')	
	cat(paste("\tData:", cimplAnalysis@data.name, "</br>\n"), file=index_file)
	cat(paste("\tChromosomes:", paste(cimplAnalysis@chromosomes, collapse=', '),"</br>\n"), file=index_file)
	cat(paste("\tScales:", paste(cimplAnalysis@scales, collapse=', '),"</br>\n"), file=index_file)
	cat(paste("\tSystem:", cimplAnalysis@system, "</br>\n"), file=index_file)
	cat(paste("\tSpecificity pattern:", cimplAnalysis@specificity.pattern, "</br>\n"), file=index_file)
	cat(paste("\tLocal hopping correction:", cimplAnalysis@lhc.method, "</br>\n"), file=index_file)
	cat(paste("\tNumber of iterations:", cimplAnalysis@n_iterations, "</br>\n"), file=index_file)
	cat(paste("\tSample period:", cimplAnalysis@sample_period, "</br>\n"), file=index_file)
	cat(paste("\tMaximum sample points:", cimplAnalysis@max_sample_points, "</br>\n"), file=index_file)
	cat(paste("\tAlpha:", alpha, "</br>\n"), file=index_file)
	cat(paste("\tMultiple testing:", mul.test, "</br>\n"), file=index_file)	
	cat(paste("\tNumber of CISs:", dim(ciss)[1], "</br></br>\n\n"), file=index_file)
	close(index_file)
	
}

.getUrl <- function(chr, bpLim, prefix='http://www.ensembl.org/Mus_musculus/Location/View?r=') {
	# example: "http://www.ensembl.org/Mus_musculus/Location/View?r=18:34210098-34628339"
	paste( prefix, substring(chr, 4), ':', bpLim[1], '-', bpLim[2], sep='' )
}