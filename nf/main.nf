params.input_dir = '/Users/james/Documents/avian-influenza/fasta'	

process input_files {

	publishDir "setup"

	input:
	path input_dir
	val n_files

	output:
	path 'concat.fa'

	script:
	"""
	find ${input_dir} -name '*HA_cns.fa' 
	"""
}

// process cat {
// 	input:
// 	path f

// 	script:
// 	"""
// 	cat ${f}
// 	"""
// }

workflow {
	println 'begin'
	def input_dir = file(params.input_dir)
	

	input_files(input_dir, 10)
	// cat('setup/concat.fa')
	println 'done'
}
