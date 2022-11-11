version 1.0

import "align_fastq_wkf.wdl" as alignFastq
import "merge_bams_wkf.wdl" as mergeBams
import "../wdl_structs.wdl"
import "qc_wkf.wdl" as qc
import "../tasks/bam_cram_conversion.wdl" as cramConversion

workflow Preprocess {
    # command
    #   align FASTQ files
    #   and
    #   merge lane level BAMs
    input {
        Boolean external = false
        Boolean highMem = false

        Array[Fastqs] listOfFastqPairs
        Boolean trim = true
        BwaMem2Reference bwamem2Reference
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
        File markerBedFile

        # resources
        #    prep flowcell
        Int threads = 16
        Int samtoolsThreads = 4
        Int bwamem2Threads = 32
    }

    Int bwamem2MemLow = 32
    Int novosortMemLow = 32

    if (highMem) {
        Int novosortMemHigh = 80
        Int bwamem2MemHigh = 48
    }

    Int novosortMem = select_first([novosortMemHigh, novosortMemLow])
    Int bwamem2Mem = select_first([bwamem2MemHigh, bwamem2MemLow])

    call alignFastq.AlignFastq {
        input:
            listOfFastqPairs = listOfFastqPairs,
            trim = trim,
            adaptersFa = adaptersFa,
            bwamem2Reference = bwamem2Reference,
            bwamem2Mem = bwamem2Mem,
            threads = samtoolsThreads,
            bwamem2Threads = bwamem2Threads
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
            markerBedFile = markerBedFile,
            outputDir = "Sample_~{sampleId}/qc"
    }

    Int bamSize = ceil(size(MergeBams.finalBam.bam, "GB"))
    Int bamToCramMemLow = 8
    Int bamToCramThreadsLow = 8
    if (bamSize > 200) {
        Int bamToCramMemHigh = 16
        Int bamToCramThreadsHigh = 12
    }
    Int bamToCramMem = select_first([bamToCramMemHigh, bamToCramMemLow])
    Int bamToCramThreads = select_first([bamToCramThreadsHigh, bamToCramThreadsLow])
    call cramConversion.SamtoolsBamToCram as bamToCram {
        input:
            inputBam = MergeBams.finalBam,
            referenceFa = referenceFa,
            sampleId = sampleId,
            threads = bamToCramThreads,
            memoryGb = bamToCramMem,
            diskSize = (ceil(size(MergeBams.finalBam.bam, "GB") * 2)) + 20 # 0.7 is estimated cram size
    }

    output {
        Bam finalBam = MergeBams.finalBam
        Cram finalCram = bamToCram.cramInfo.finalCram
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
        File binestSex = QcMetrics.binestSex
        File normCoverageByChrPng = QcMetrics.normCoverageByChrPng
        # Dedup metrics.
        File collectWgsMetricsPreBqsr = MergeBams.collectWgsMetricsPreBqsr
        File qualityDistributionPdfPreBqsr = MergeBams.qualityDistributionPdfPreBqsr
        File qualityByCycleMetricsPreBqsr = MergeBams.qualityByCycleMetricsPreBqsr
        File qualityByCyclePdfPreBqsr = MergeBams.qualityByCyclePdfPreBqsr
        File qualityDistributionMetricsPreBqsr = MergeBams.qualityDistributionMetricsPreBqsr
        File dedupLog = MergeBams.dedupLog
        File pileupsConpair = QcMetrics.pileupsConpair
    }

}
