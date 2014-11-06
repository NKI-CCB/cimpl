#setClass('cXYcoords', representation(x='numeric', y='numeric'))
#setClass('cXYZcoords', representation(x='numeric', y='numeric', z='numeric'))

setClass('cimplObject', representation(
	data = 'data.frame',          # containing $chromosome and $location columns + arbitrary annotation columns 
	kse = 'list',            # x, y
	peaks = 'list',         # x, y, p_vals
	null_peaks = 'list',     # x, y
	null_cdf = 'list',      # x, y, z
	bg_density = 'list',     # x, y
	n_peaks = 'integer',
	scale = 'numeric',
	chromosome = 'character'
))

setClass('cimplAnalysis', representation(
	cimplObjects = 'list',           # a matrix of cimplObjects ( chromosome X scale )
	chromosomes = 'character',
	scales = 'numeric',
	system = 'character',
	specificity.pattern = 'character',
	lhc.method = 'character',
	n_iterations = 'numeric',
	sample_period = 'numeric',
	max_sample_points = 'numeric',
	data.name = 'character',
	call = 'call'
))
