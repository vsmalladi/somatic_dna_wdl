version 1.0

import "prep_flowcell.wdl" as prepFlowcell
import "merge_flowcells.wdl" as mergeFlowcells

workflow RunFlowcell {
	# command 
	# 	prep flowcell
	Array[Fastqs]+ listOfFastqPairs
	BwaReference bwaReference
	#	command merge flowcell
	String sampleId
	IndexedReference indexedReference
	Directory tempDir
	IndexedVcf MillsAnd1000G
	IndexedVcf Indels
	IndexedVcf DbSnp
	File chromFile
	# resources
	#	prep flowcell
	Int mem
	Int threads
	String bwaDockerImage
	String shortAlignDockerImage
	String gatkDockerImage
	#	merge flowcell
	String novosortDockerImage
	String samtoolsDockerImage
	
	scatter(fastqs in listOfFastqPairs {
		call prepFlowcell.AlignBwaMem {
			input:
				fastqs = fastqs,
				bwaReference = bwaReference,
				mem = mem,
				threads = threads,
				dockerImage = bwaDockerImage
		}
		
		call prepFlowcell.ShortAlignMark {
			input:
				laneBam = prepFlowcell.AlignBwaMem.laneBam,
				bamBase = fastqs.bamBase,
				mem = mem,
				dockerImage = shortAlignDockerImage
		}
		
		call prepFlowcell.Fixmate {
			input:
				laneBamMark = prepFlowcell.ShortAlignMark.laneBamMark,
				bamBase = fastqs.bamBase,
				mem = mem,
				dockerImage = gatkDockerImage
		}
	
	output {
		Bam prepFlowcell.Fixmate.laneFixmateBam
	}
	
	call mergeFlowcells.NovosortMarkDup {
		input:
			laneBams = prepFlowcell.Fixmate.laneFixmateBam,
			laneBamsLists = prepFlowcell.Fixmate.laneFixmateBamPath,
			temp_dir = tempDir,
			sampleId = sampleId
			dockerImage = novosortDockerImage
	}
	
	call mergeFlowcells.IndexBam as novosortMarkDupIndexed {
		input:
			bam = mergeFlowcells.NovosortMarkDup.mergedDedupBamOnly,
			indexedReference = indexedReference,
			dockerImage = samtoolsDockerImage
	}
	
	call mergeFlowcells.Bqsr38 {
		input:
			tempDir = tempDir
			mergedDedupBam = novosortMarkDupIndexed.indexedBam
			MillsAnd1000G = MillsAnd1000G
			Indels = Indels
			DbSnp = DbSnp
			chromFile = chromFile
			sampleId = sampleId
			
	}
	call mergeFlowcells.PrintReads {
		input:
			tempDir = tempDir
			mergedDedupBam = novosortMarkDupIndexed.indexedBam
			recalGrp = mergeFlowcells.Bqsr38.recalGrp
			sampleId = sampleId
			
	}
	output {
		Bam mergeFlowcells.PrintReads.finalBam
	}
}
