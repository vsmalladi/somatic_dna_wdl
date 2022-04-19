version 1.0

import "align_fastq.wdl" as alignFastq
import "../wdl_structs.wdl"

workflow AlignFastq {
    # command
    #   align FASTQ files
    input {
        Array[Fastqs] listOfFastqPairs
        BwaReference bwaReference
        File adaptersFa
        Boolean trim = true
        # resources
        Int bwaMem
        Int threads
        Int bwaThreads
    }

    Int additionalDiskSize = 100
    scatter(fastqs in listOfFastqPairs) {
        Int fastqsSize = ceil(size(fastqs.fastqR1, "GB") + size (fastqs.fastqR2, "GB"))

        # Assumption: The bam will be about the same size as the input fastqs. So double
        # the size to account for input and output.
        Int diskSize = fastqsSize * 2 + additionalDiskSize
        
        if (trim) {
            call alignFastq.Skewer {
                input:
                    fastqs = fastqs,
                    adaptersFa = adaptersFa,
                    diskSize = diskSize
            }
        }

        Fastqs fastqsAlign = select_first([Skewer.skewerFastqs, fastqs])

        call alignFastq.AlignBwaMem {
            input:
                fastqsAlign = fastqsAlign,
                bwaReference = bwaReference,
                memoryGb = bwaMem,
                threads = threads,
                bwaThreads = bwaThreads,
                diskSize = diskSize
        }

        call alignFastq.ShortAlignMark {
            input:
                laneBam = AlignBwaMem.laneBam,
                bamBase = fastqs.readgroupId,
                memoryGb = 16,
                diskSize = diskSize
        }

        call alignFastq.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.readgroupId,
                memoryGb = 8,
                diskSize = diskSize
        }
        Int fixmate_bam_size = ceil(size(Fixmate.laneFixmateBam, "GB"))
    }

    output {
        Array[File] laneFixmateBam = Fixmate.laneFixmateBam
        Array[Int] laneFixmateBamSizes = fixmate_bam_size
    }
}
