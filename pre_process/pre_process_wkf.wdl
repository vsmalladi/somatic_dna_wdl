version 1.0

import "align_fastq_wkf.wdl" as alignFastq
import "merge_bams_wkf.wdl" as mergeBams
import "../wdl_structs.wdl"
import "qc_wkf.wdl" as qc
#import "../tasks/bam_cram_conversion.wdl" as cramConversion

workflow Preprocess {
    # command
    #   align FASTQ files
    #   and
    #   merge lane level BAMs
    input {
        Boolean external = false

        Array[Fastqs] listOfFastqPairs
        Boolean trim = true
        BwaReference bwaReference
        File adaptersFa
        #    command merge flowcell
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        #IndexedTable callRegions
        File callRegions
        IndexedReference referenceFa
        File hsMetricsIntervals
        File randomIntervals
        File chromLengths
        IndexedVcf gnomadBiallelic

        # resources
        #    prep flowcell
        Int bwaMem = 86
        Int novosortMem = 32
        Int threads = 16
        Int bwaThreads = 96
    }

    call alignFastq.AlignFastq {
        input:
            listOfFastqPairs = listOfFastqPairs,
            trim = trim,
            adaptersFa = adaptersFa,
            bwaReference = bwaReference,
            bwaMem = bwaMem,
            threads = threads,
            bwaThreads = bwaThreads
    }

    call mergeBams.MergeBams {
        input:
            external = external,
            laneFixmateBams = AlignFastq.laneFixmateBam,
            sampleId = sampleId,
            MillsAnd1000G = MillsAnd1000G,
            Indels = Indels,
            dbsnp = dbsnp,
            callRegions = callRegions,
            referenceFa = referenceFa,
            randomIntervals = randomIntervals,
            qcDir = "Sample_~{sampleId}/qc",
            novosortMem = novosortMem,
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

#    call cramConversion.SamtoolsBamToCram as bamToCram {
#        input:
#            inputBam = MergeBams.finalBam,
#            referenceFa = referenceFa,
#            sampleId = sampleId,
#            diskSize = (ceil(size(MergeBams.finalBam.bam, "GB")) * 2) + 50
#    }

    output {
        Bam finalBam = MergeBams.finalBam
        Cram finalCram = MergeBams.finalCram
#        Cram finalCram = bamToCram.finalCram
        File alignmentSummaryMetrics = QcMetrics.alignmentSummaryMetrics
        File qualityByCyclePdf = QcMetrics.qualityByCyclePdf
        File baseDistributionByCycleMetrics = QcMetrics.baseDistributionByCycleMetrics
        File qualityByCycleMetrics = QcMetrics.qualityByCycleMetrics
        File baseDistributionByCyclePdf = QcMetrics.baseDistributionByCyclePdf
        File qualityDistributionPdf = QcMetrics.qualityDistributionPdf
        File qualityDistributionMetrics = QcMetrics.qualityDistributionMetrics
        File insertSizeHistogramPdf = QcMetrics.insertSizeHistogramPdf
        File insertSizeMetrics = QcMetrics.insertSizeMetrics
        File gcBiasMetrics = QcMetrics.gcBiasMetrics
        File gcBiasSummary = QcMetrics.gcBiasSummary
        File gcBiasPdf = QcMetrics.gcBiasPdf
        File flagStat = QcMetrics.flagStat
        File hsMetrics = QcMetrics.hsMetrics
        File hsMetricsPerTargetCoverage = QcMetrics.hsMetricsPerTargetCoverage
        File hsMetricsPerTargetCoverageAutocorr = QcMetrics.hsMetricsPerTargetCoverageAutocorr
        File autocorroutput1100 = QcMetrics.autocorroutput1100
        File collectOxoGMetrics = QcMetrics.collectOxoGMetrics
        File collectWgsMetrics = QcMetrics.collectWgsMetrics
        File binestCov = QcMetrics.binestCov
        File normCoverageByChrPng = QcMetrics.normCoverageByChrPng
        # Dedup metrics.
        File collectWgsMetricsPreBqsr = MergeBams.collectWgsMetricsPreBqsr
        File qualityDistributionPdfPreBqsr = MergeBams.qualityDistributionPdfPreBqsr
        File qualityByCycleMetricsPreBqsr = MergeBams.qualityByCycleMetricsPreBqsr
        File qualityByCyclePdfPreBqsr = MergeBams.qualityByCyclePdfPreBqsr
        File qualityDistributionMetricsPreBqsr = MergeBams.qualityDistributionMetricsPreBqsr
        File dedupLog = MergeBams.dedupLog
    }

}
