version 1.0

import "align_fastq.wdl" as alignFastq
import "../wdl_structs.wdl"

workflow AlignFastq {
    # command
    #   align FASTQ files
    input {
        Array[Fastqs]+ listOfFastqPairs
        BwaReference bwaReference
        # resources
        Int bwaMem
        Int threads
        String bwaDockerImage
        String shortAlignDockerImage
        String gatkDockerImage

    }

    scatter(fastqs in listOfFastqPairs) {
        call alignFastq.AlignBwaMem {
            input:
                fastqs = fastqs,
                bwaReference = bwaReference,
                mem = bwaMem,
                threads = threads,
                dockerImage = bwaDockerImage
        }

        call alignFastq.ShortAlignMark {
            input:
                laneBam = AlignBwaMem.laneBam,
                bamBase = fastqs.readgroupId,
                mem = 16,
                dockerImage = shortAlignDockerImage
        }

        call alignFastq.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.readgroupId,
                mem = 8,
                dockerImage = gatkDockerImage
        }
    }

    output {
        Array[File] laneFixmateBam = Fixmate.laneFixmateBam
    }
}
