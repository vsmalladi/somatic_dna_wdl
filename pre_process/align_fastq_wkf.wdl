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
        Int minimapMem
        Int threads
        Int minimapThreads
    }

    Int additionalDiskSize = 10
    scatter(fastqs in listOfFastqPairs) {
        Int fastqsSize = ceil(size(fastqs.fastqR1, "GB") + size (fastqs.fastqR2, "GB"))

        # Assumption: The bam will be about the same size as the input fastqs. So double
        # the size to account for input and output.
        Int skewerDiskSize = (2 * fastqsSize) + additionalDiskSize
        Int alignDiskSize = (3 * fastqsSize) + additionalDiskSize
        Int fixmateDiskSize = (4 * fastqsSize) + additionalDiskSize
        
        if (trim) {
            call alignFastq.Skewer {
                input:
                    fastqs = fastqs,
                    adaptersFa = adaptersFa,
                    diskSize = skewerDiskSize
            }
        }

        Fastqs fastqsAlign = select_first([Skewer.skewerFastqs, fastqs])

        # call alignFastq.AlignBwaMem {
        #     input:
        #         fastqsAlign = fastqsAlign,
        #         bwaReference = bwaReference,
        #         memoryGb = bwaMem,
        #         threads = threads,
        #         bwaThreads = bwaThreads,
        #         diskSize = alignDiskSize
        # }

        call alignFastq.AlignMinimap2 {
            input:
                fastqsAlign = fastqsAlign,
                bwaReference = bwaReference,
                memoryGb = minimapMem,
                threads = threads,
                minimapThreads = minimapThreads,
                diskSize = alignDiskSize
        }

        call alignFastq.ShortAlignMark {
            input:
                laneBam = AlignMinimap2.laneBam,
                bamBase = fastqs.readgroupId,
                memoryGb = 16,
                diskSize = alignDiskSize
        }

        call alignFastq.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.readgroupId,
                memoryGb = 8,
                diskSize = fixmateDiskSize
        }
    }

    output {
        Array[File] laneFixmateBam = Fixmate.laneFixmateBam
    }
}
