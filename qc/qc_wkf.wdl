version 1.0


import "qc.wdl" as qc
import "../wdl_structs.wdl"

workflow QcMetrics {
    # command
    input {
        Bam finalBam
        Bam mergedDedupBam
        IndexedReference referenceFa
        String sampleId
        File hsMetricsIntervals
        File randomIntervals
        File chromLengths
        File gnomadBiallelic
        String outputDir = "."
    }

    Int additionalDiskSize = 50
    Int diskSize = ceil(size(finalBam.bam, "GB") + size(finalBam.bamIndex, "GB")) + additionalDiskSize

    call qc.MultipleMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            outputDir = outputDir,
            diskSize = diskSize
    }

    call qc.MultipleMetricsPreBqsr {
        input:
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
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
            HsMetricsPerTargetCoverage = HsMetrics.HsMetricsPerTargetCoverage,
            sampleId = sampleId,
            outputDir = outputDir
    }

    call qc.Autocorrelations {
        input:
            HsMetricsPerTargetCoverageAutocorr = FormatHsMetrics.HsMetricsPerTargetCoverageAutocorr,
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

    call qc.Binest as BinestCov {
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
            binestCov = BinestCov.binestCov
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
            CollectWgsMetrics.CollectWgsMetrics,
            BinestCov.binestCov,
            PlotBinCov.normCoverageByChrPng,
            Pileup.pileupsTable,
            CalculateContamination.contaminationTable
            ]
    }

}
