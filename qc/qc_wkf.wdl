version 1.0


import "qc.wdl" as qc
import "../wdl_structs.wdl"

workflow QcMetrics {
    # command
    input {
        Bam finalBam
        Bam mergedDedupBam
        String MultipleMetricsBase
        IndexedReference referenceFa
        String sampleId
        String MultipleMetricsBasePreBqsrBasename
        File hsMetricsIntervals
        File randomIntervals
        File chromLengths
    }

    Int additionalDiskSize = 50
    Int diskSize = ceil(size(finalBam.bam, "GB") + size(finalBam.bamIndex, "GB")) + additionalDiskSize

    call qc.MultipleMetrics {
        input:
            MultipleMetricsBase = MultipleMetricsBase,
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            diskSize = diskSize
    }

    call qc.MultipleMetricsPreBqsr {
        input:
            MultipleMetricsBasePreBqsrBasename = MultipleMetricsBasePreBqsrBasename,
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
            sampleId = sampleId,
            diskSize = diskSize
    }

    call qc.CollectGcBiasMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            diskSize = diskSize,
    }

    call qc.Flagstat {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            diskSize = diskSize

    }

    call qc.HsMetrics {
        input:
            referenceFa = referenceFa,
            hsMetricsIntervals = hsMetricsIntervals,
            finalBam = finalBam,
            sampleId = sampleId,
            diskSize = diskSize
    }

    call qc.FormatHsMetrics {
        input:
            HsMetricsPerTargetCoverage = HsMetrics.HsMetricsPerTargetCoverage,
            sampleId = sampleId
    }

    call qc.Autocorrelations {
        input:
            HsMetricsPerTargetCoverageAutocorr = FormatHsMetrics.HsMetricsPerTargetCoverageAutocorr,
            sampleId = sampleId,
    }

    call qc.CollectOxoGMetricsWgs {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            diskSize = diskSize
    }

    call qc.CollectWgsMetricsWgsDecoy {
        input:
            referenceFa = referenceFa,
            randomIntervals = randomIntervals,
            finalBam = finalBam,
            sampleId = sampleId,
            diskSize = diskSize
    }

    call qc.Binest as BinestCov {
        input:
            finalBam = finalBam,
            sampleId = sampleId,
            diskSize = diskSize
    }

    call qc.PlotBinCov {
        input:
            chromLengths = chromLengths,
            sampleId = sampleId,
            binestCov = BinestCov.binestCov
    }

    output {
        Array[File] QcFiles = [
            MultipleMetrics.alignmentSummaryMetrics,
            MultipleMetrics.qualityByCyclePdf,
            MultipleMetrics.baseDistributionByCycleMetrics,
            MultipleMetrics.qualityByCycleMetrics,
            MultipleMetrics.baseDistributionByCyclePdf,
            MultipleMetrics.qualityDistributionPdf,
            MultipleMetrics.qualityDistributionMetrics,
            MultipleMetrics.insertSizeHistogramPdf,
            MultipleMetrics.insertSizeMetrics,
            MultipleMetricsPreBqsr.qualityDistributionPdfPreBqsr,
            MultipleMetricsPreBqsr.qualityByCycleMetricsPreBqsr,
            MultipleMetricsPreBqsr.qualityByCyclePdfPreBqsr,
            MultipleMetricsPreBqsr.qualityDistributionMetricsPreBqsr,
            CollectGcBiasMetrics.gcBiasMetrics,
            CollectGcBiasMetrics.gcBiasSummary,
            CollectGcBiasMetrics.gcBiasPdf,
            Flagstat.FlagStat,
            HsMetrics.HsMetrics,
            HsMetrics.HsMetricsPerTargetCoverage,
            FormatHsMetrics.HsMetricsPerTargetCoverageAutocorr,
            Autocorrelations.autocorroutput1100,
            CollectOxoGMetricsWgs.CollectOxoGMetrics,
            CollectWgsMetricsWgsDecoy.CollectWgsMetrics,
            BinestCov.binestCov,
            PlotBinCov.normCoverageByChrPng
            ]
    }

        # File alignmentSummaryMetrics = MultipleMetrics.alignmentSummaryMetrics
        # File qualityByCyclePdf = MultipleMetrics.qualityByCyclePdf
        # File baseDistributionByCycleMetrics = MultipleMetrics.baseDistributionByCycleMetrics
        # File qualityByCycleMetrics = MultipleMetrics.qualityByCycleMetrics
        # File baseDistributionByCyclePdf = MultipleMetrics.baseDistributionByCyclePdf
        # File qualityDistributionPdf = MultipleMetrics.qualityDistributionPdf
        # File qualityDistributionMetrics = MultipleMetrics.qualityDistributionMetrics
        # File insertSizeHistogramPdf = MultipleMetrics.insertSizeHistogramPdf
        # File insertSizeMetrics = MultipleMetrics.insertSizeMetrics

        # File qualityDistributionPdfPreBqsr = MultipleMetricsPreBqsr.qualityDistributionPdfPreBqsr
        # File qualityByCycleMetricsPreBqsr = MultipleMetricsPreBqsr.qualityByCycleMetricsPreBqsr
        # File qualityByCyclePdfPreBqsr = MultipleMetricsPreBqsr.qualityByCyclePdfPreBqsr
        # File qualityDistributionMetricsPreBqsr = MultipleMetricsPreBqsr.qualityDistributionMetricsPreBqsr



}
