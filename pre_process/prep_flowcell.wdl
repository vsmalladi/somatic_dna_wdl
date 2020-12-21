version 1.0

task AlignBwaMem {
   input {
	# command
	Fastqs fastqs
	BwaReference bwaReference
	String laneBamPath = "~{fastqs.bamBase}.readgroup.bam"
	# resources
	Int mem
	Int threads
	String dockerImage
   }

	command {
		set -e -o pipefail
		bwa mem \
		-Y \
		-K 100000000 \
		-t ~{threads} \
		-R ~{fastqs.rgInfo} \
		~{bwaReference.fasta} \
		~{fastqs.fastqR1} \
		~{fastqs.fastqR2} \
		| samtools view \
		-Shb \
		-o ~{laneBamPath} \
		-
		
	}
	
	output {
		Bam laneBam = object {
				bam : laneBamPath, 
				bamIndex : sub(laneBamPath, ".bam$", ".bai")
			}
	}
	
	runtime {
		cpu : threads
		memory : mem
		docker : dockerImage
	}
}

task ShortAlignMark {
	input {
		# command
		Bam laneBam
		String bamBase
		String laneBamMarkPath = "~{bamBase}.readgroup_mark.bam"
		# resources
		Int mem
		String dockerImage
	}
	
	command {
		set -e -o pipefail
		filter_bam \
		-I ~{laneBam.bam} \
		-A1 30 \
		-A2 30 \
		-o ~{laneBamMarkPath} \
		| samtools view \
		-b \
		-o ~{laneBamMarkPath} \
		-
	}
	
	output {
		Bam laneBamMark = object {
				bam : laneBamMarkPath, 
				bamIndex : sub(laneBamMarkPath, ".bam$", ".bai")
			}
	}
	
	runtime {
		memory : mem
		docker : dockerImage
	}
}

task Fixmate {
	input {
		#command
		Bam laneBamMark
		String bamBase
		String laneFixmateBamPath = "~{bamBase}.readgroup_fixmate.bam"
		# resources
		Int mem
		String dockerImage
	}
	
	command {
		gatk \
		FixMateInformation \
		--java-options \
		-Xmx24g -XX:ParallelGCThreads=1 \
		--MAX_RECORDS_IN_RAM 2000000 \
		--VALIDATION_STRINGENCY SILENT \
		--ADD_MATE_CIGAR true \
		--ASSUME_SORTED true \
		-I laneBamMark.bam \
		-O laneFixmateBamPath
	}
	
	output {
		Bam laneFixmateBam = object {
				bam : laneFixmateBamPath, 
				bamIndex : sub(laneFixmateBamPath, ".bam$", ".bai")
			}
	}
	
	runtime {
		memory : mem
		docker : dockerImage
	}
}

struct Bam {
	File bam
	File bamIndex
	String? md5sum
}

struct BwaReference {
	File fasta
	File sa
	File pac
	File bwt
	File ann
	File alt
	File amb	
}

struct IndexedReference {
	File fasta
	File dict
	File index
}

struct Fastqs {
	File fastqR1
	String bamBase
	String rgInfo
	String? md5sumR1
	File fastqR2
	String? md5sumR1
}