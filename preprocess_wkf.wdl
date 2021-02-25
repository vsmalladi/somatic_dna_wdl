version 1.0

import "pre_process/align_fastq_wkf.wdl" as alignFastqWf
import "pre_process/merge_bams_wkf.wdl" as mergeBams
import "qc/qc_wkf.wdl" as qc
import "wdl_structs.wdl"

workflow Preprocess {
    # resources
    Int bwaMem=16
    Int threads=8

    input {
        Array[Fastqs] listOfFastqPairs
        BwaReference bwaReference
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        IndexedReference referenceFa
        File hsMetricsIntervals
        File randomIntervals
        File chromLengths
        Int bwaMem = 24
        Int novosortMem = 80
        Int threads = 8
    }

    call alignFastqWf.AlignFastq {
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
            DbSnp = dbsnp,
            indexedReference = referenceFa,
            randomIntervals = randomIntervals,
            mem = novosortMem,
            threads = threads
    }

    call qc.QcMetrics {
        input:
            finalBam = MergeBams.finalBam,
            referenceFa = referenceFa,
            sampleId = sampleId,
            hsMetricsIntervals = hsMetricsIntervals,
            randomIntervals = randomIntervals,
            chromLengths = chromLengths
    }

    output {
        Bam finalBam = MergeBams.finalBam
        Array[File] QcFiles = QcMetrics.QcFiles
    }
}
