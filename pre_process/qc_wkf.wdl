version 1.0


import "qc.wdl" as qc
import "../wdl_structs.wdl"

workflow QcMetrics {
    # command
    input {
        Bam finalBam
        IndexedReference referenceFa
        String sampleId
        File hsMetricsIntervals
        File randomIntervals
        File chromLengths
        IndexedVcf gnomadBiallelic
        String outputDir = "."
    }

    Int additionalDiskSize = 10
    Int diskSize = ceil((size(finalBam.bam, "GB") + size(finalBam.bamIndex, "GB"))) +
                      additionalDiskSize

    call qc.MultipleMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.CollectGcBiasMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize,
    }

    call qc.Flagstat {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize

    }

    call qc.HsMetrics {
        input:
            referenceFa = referenceFa,
            hsMetricsIntervals = hsMetricsIntervals,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.FormatHsMetrics {
        input:
            hsMetricsPerTargetCoverage = HsMetrics.hsMetricsPerTargetCoverage,
            sampleId = sampleId,
            outputDir = outputDir
    }

    call qc.Autocorrelations {
        input:
            hsMetricsPerTargetCoverageAutocorr = FormatHsMetrics.hsMetricsPerTargetCoverageAutocorr,
            sampleId = sampleId,
            outputDir = outputDir
    }

    call qc.CollectOxoGMetricsWgs {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.CollectWgsMetrics {
        input:
            referenceFa = referenceFa,
            randomIntervals = randomIntervals,
            inputBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.Binest as Binest {
        input:
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.PlotBinCov {
        input:
            chromLengths = chromLengths,
            sampleId = sampleId,
            outputDir = outputDir,
            binestCov = Binest.binestCov
    }

    call qc.Pileup {
        input:
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            gnomadBiallelic = gnomadBiallelic,
            diskSize = diskSize
    }

    call qc.CalculateContamination {
        input:
            sampleId = sampleId,
            pileupsTable = Pileup.pileupsTable,
            outputDir = outputDir,
            diskSize = diskSize
    }

    output {
        File alignmentSummaryMetrics = MultipleMetrics.alignmentSummaryMetrics
        File qualityByCyclePdf = MultipleMetrics.qualityByCyclePdf
        File baseDistributionByCycleMetrics = MultipleMetrics.baseDistributionByCycleMetrics
        File qualityByCycleMetrics = MultipleMetrics.qualityByCycleMetrics
        File baseDistributionByCyclePdf = MultipleMetrics.baseDistributionByCyclePdf
        File qualityDistributionPdf = MultipleMetrics.qualityDistributionPdf
        File qualityDistributionMetrics = MultipleMetrics.qualityDistributionMetrics
        File insertSizeHistogramPdf = MultipleMetrics.insertSizeHistogramPdf
        File insertSizeMetrics = MultipleMetrics.insertSizeMetrics
        File gcBiasMetrics = CollectGcBiasMetrics.gcBiasMetrics
        File gcBiasSummary = CollectGcBiasMetrics.gcBiasSummary
        File gcBiasPdf = CollectGcBiasMetrics.gcBiasPdf
        File flagStat = Flagstat.flagStat
        File hsMetrics = HsMetrics.hsMetrics
        File hsMetricsPerTargetCoverage = HsMetrics.hsMetricsPerTargetCoverage
        File hsMetricsPerTargetCoverageAutocorr = FormatHsMetrics.hsMetricsPerTargetCoverageAutocorr
        File autocorroutput1100 = Autocorrelations.autocorroutput1100
        File collectOxoGMetrics = CollectOxoGMetricsWgs.collectOxoGMetrics
        File collectWgsMetrics = CollectWgsMetrics.collectWgsMetrics
        File binestCov = Binest.binestCov
        File binestSex = Binest.binestSex
        File normCoverageByChrPng = PlotBinCov.normCoverageByChrPng
    }

}
