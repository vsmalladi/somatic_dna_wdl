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
        
        File fastqR1 = select_first([Skewer.fastqOutR1, fastqs.fastqR1])
        File fastqR2 = select_first([Skewer.fastqOutR2, fastqs.fastqR2])
        
        call alignFastq.AlignBwaMem {
            input:
                fastqs = fastqs,
                fastqR1 = fastqR1,
                fastqR2 = fastqR2,
                bwaReference = bwaReference,
                mem = bwaMem,
                threads = threads,
                bwaThreads = bwaThreads,
                diskSize = diskSize
        }

        call alignFastq.ShortAlignMark {
            input:
                laneBam = AlignBwaMem.laneBam,
                bamBase = fastqs.readgroupId,
                mem = 16,
                diskSize = diskSize
        }

        call alignFastq.Fixmate {
            input:
                laneBamMark = ShortAlignMark.laneBamMark,
                bamBase = fastqs.readgroupId,
                mem = 8,
                diskSize = diskSize
        }
        Int fixmate_bam_size = ceil(size(Fixmate.laneFixmateBam, "GB"))
    }

    output {
        Array[File] laneFixmateBam = Fixmate.laneFixmateBam
        Array[Int] laneFixmateBamSizes = fixmate_bam_size
    }
}
