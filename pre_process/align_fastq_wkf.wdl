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
        Int mem
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
                mem = mem,
                threads = threads,
                dockerImage = bwaDockerImage
        }

        call alignFastq.ShortAlignMark {
            input:
                laneBam = AlignBwaMem.laneBam,
                bamBase = fastqs.bamBase,
                mem = mem,
                dockerImage = shortAlignDockerImage
        }

        call alignFastq.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.bamBase,
                mem = mem,
                dockerImage = gatkDockerImage
        }
    }

    output {
        Array[File] laneFixmateBam = Fixmate.laneFixmateBam
        Array[String] laneFixmateBamPath = Fixmate.laneFixmateBamPath
    }
}
