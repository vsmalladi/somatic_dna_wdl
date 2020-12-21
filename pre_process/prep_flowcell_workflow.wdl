version 1.0

import "prep_flowcell.wdl" as prepFlowcell

workflow PrepFlowcell {
	# command
	Array[Fastqs]+ listOfFastqPairs
	BwaReference bwaReference
	# resources
	Int mem
	Int threads
	String dockerImage
	
	scatter(fastqs in listOfFastqPairs {
		call prepFlowcell.AlignBwaMem {
			input:
				fastqs = fastqs,
				bwaReference = bwaReference,
				mem = mem,
				threads = threads,
				dockerImage = dockerImage
		}
		
		call prepFlowcell.ShortAlignMark {
			input:
				laneBam = prepFlowcell.AlignBwaMem.laneBam,
				bamBase = fastqs.bamBase,
				mem = mem,
				dockerImage = dockerImage
		}
		
		call prepFlowcell.Fixmate {
			input:
				laneBamMark = prepFlowcell.ShortAlignMark.laneBamMark,
				bamBase = fastqs.bamBase,
				mem = mem,
				dockerImage = dockerImage
		
		}
	
	output {
		Bam prepFlowcell.Fixmate.laneFixmateBam
	}
}