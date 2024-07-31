#Configuration file with defined variables
configfile: "config.yaml"

#Parse from config
gSizeL = config["gSizeL"]
gSizeS = config["gSizeS"]

samples = ['YMX005686', 'YMX005664', 'YMX005668']

rule all:
	input:
		expand("assemblies/final/{sample}_assembly.fasta", sample = samples)

rule filter:
	input:
		raw = "data/longReads/{sample}.fastq.gz"
	output:
		filtered = temp("data/filtered/{sample}.fastq.gz")
	params:
		targetSize = gSizeL * 40
	shell:
		"""
		module load filtlong/0.2.0
		filtlong --target_bases {params.targetSize} --min_length 1000 {input.raw} | gzip > {output.filtered}
		"""

rule assembly:
	input:
		filtered = rules.filter.output.filtered
	output:
		assembly = temp("assemblies/{sample}_project/{sample}_project.contigs.fasta")
	params:
		gsize = gSizeS
	shell:
		"""
		module load canu/2.2 gnuplot/5.2.6
		canu -p {wildcards.sample}_project \
		-d assemblies/{wildcards.sample}_project  genomeSize={params.gsize} \
		-nanopore {input.filtered} useGrid=false maxThreads=12
		"""

rule long_polish:
	input:
		assembly = rules.assembly.output.assembly,
		filtered = rules.filter.output.filtered
	output:
		polished = temp("assemblies/{sample}_medaka/consensus.fasta")
	params:
		threads = 4
	shell:
		"""
		module unload python37/3.7.0
		module load anaconda3/2021.05 bcftools/1.17 minimap2/2.24 samtools/1.16.1 htslib/1.12
		source activate medaka1.11.2
		medaka_consensus -t {params.threads} -i {input.filtered} \
		-d {input.assembly} -o assemblies/{wildcards.sample}_medaka
		"""

rule move:
	input:
		polished = rules.long_polish.output.polished
	output:
		final = "assemblies/final/{sample}_polished.fasta"
	shell:
		"""
		mv {input.polished} {output.final}
		"""

rule short_polish:
	input:
		fasta = rules.move.output.final,
		r1 = "data/shortReads/{sample}_R1_clean.fastq.gz",
		r2 = "data/shortReads/{sample}_R2_clean.fastq.gz"
	output:
		polished = "assemblies/final/{sample}/{sample}_final.fasta"
	shell:
		"""
		module load hapog/1.3.8 bwa/0.7.15 samtools/1.9
		/cm/shared/apps/hapog/1.3.8/hapog.py --genome {input.fasta} \
		--pe1 {input.r1} \
		--pe2 {input.r2} \
		--output assemblies/final/{wildcards.sample}_temp \
		--threads 8 -u
		mv assemblies/final/{wildcards.sample}_temp/* assemblies/final/{wildcards.sample}/
		cp assemblies/final/{wildcards.sample}/hapog_results/hapog.fasta {output.polished}
		rm -r assemblies/final/{wildcards.sample}_temp/
		"""

rule report:
	input:
		final = rules.short_polish.output.polished
	output:
		rep = "reports/{sample}/report.pdf"
	params:
		threads = 2
	shell:
		"""
		module load quast/4.6.3
		/cm/shared/apps/quast/quast-4.6.3/quast.py -o reports/{wildcards.sample} -t {params.threads} {input.final}
		"""
