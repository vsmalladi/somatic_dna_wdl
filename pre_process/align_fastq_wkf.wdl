version 1.0

import "prep_flowcell.wdl" as prepFlowcell
import "../wdl_structs.wdl"

workflow AlignFastq {
    # command 
    #   align FASTQ files
    input {
        Array[Fastqs]+ listOfFastqPairs
        BwaReference bwaReference
        # resources
        Int mem
        Int threads
        String bwaDockerImage
        String shortAlignDockerImage
        String gatkDockerImage

    }
    
    scatter(fastqs in listOfFastqPairs) {
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
                laneBam = AlignBwaMem.laneBam,
                bamBase = fastqs.bamBase,
                mem = mem,
                dockerImage = shortAlignDockerImage
        }
        
        call prepFlowcell.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.bamBase,
                mem = mem,
                dockerImage = gatkDockerImage
        }
    }
    
    output {
        Array[Bam] laneFixmateBam = Fixmate.laneFixmateBam
        Array[String] laneFixmateBamPath = Fixmate.laneFixmateBamPath
    }
}
