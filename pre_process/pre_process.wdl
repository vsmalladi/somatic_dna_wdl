version 1.0

import "align_fastq_wkf.wdl" as alignFastq
import "merge_bams_wkf.wdl" as mergeBams
import "../wdl_structs.wdl"

workflow Preprocess {
    # command
    #   align FASTQ files
    #   and
    #   merge lane level BAMs
    input {
        Array[Fastqs] listOfFastqPairs
        BwaReference bwaReference
        #    command merge flowcell
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        IndexedTable callRegions
        IndexedReference referenceFa
        # resources
        #    prep flowcell
        Int bwaMem = 24
        Int novosortMem = 80
        Int threads = 8
    }

    call alignFastq.AlignFastq {
        input:
            listOfFastqPairs = listOfFastqPairs,
            bwaReference = bwaReference,
            bwaMem = bwaMem,
            threads = threads
    }

    call mergeBams.MergeBams {
        input:
            laneFixmateBams = AlignFastq.laneFixmateBam,
            sample_bam_sizes = AlignFastq.laneFixmateBamSizes,
            sampleId = sampleId,
            MillsAnd1000G = MillsAnd1000G,
            Indels = Indels,
            dbsnp = dbsnp,
            callRegions = callRegions,
            referenceFa = referenceFa,
            mem = novosortMem,
            threads = threads,
            sample_bam_sizes = size(AlignFastq.laneFixmateBam)
    }

    output {
        Bam mergedDedupBam = MergeBams.mergedDedupBam
        Bam finalBam = MergeBams.finalBam
    }
}
