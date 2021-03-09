version 1.0

import "align_fastq_wkf.wdl" as alignFastq
import "merge_bams_wkf.wdl" as mergeBams
import "../wdl_structs.wdl"
import "../qc/qc_wkf.wdl" as qc

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
        File hsMetricsIntervals
        File randomIntervals
        File chromLengths
        File gnomadBiallelic

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
            randomIntervals = randomIntervals,
            qcDir = "Sample_~{sampleId}/qc",
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
            chromLengths = chromLengths,
            gnomadBiallelic = gnomadBiallelic,
            outputDir = "Sample_~{sampleId}/qc"
    }

    output {
        Bam finalBam = MergeBams.finalBam
        Array[File] QcFiles = QcMetrics.QcFiles
        # Dedup metrics.
        File collectWgsMetricsPreBqsr = MergeBams.collectWgsMetricsPreBqsr
        File qualityDistributionPdPreBqsr = MergeBams.qualityDistributionPdPreBqsr
        File qualityByCycleMetricsPreBqsr = MergeBams.qualityByCycleMetricsPreBqsr
        File qualityByCyclePdfPreBqsr = MergeBams.qualityByCyclePdfPreBqsr
        File qualityDistributionMetricsPreBqsr = MergeBams.qualityDistributionMetricsPreBqsr
    }

}
