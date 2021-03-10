version 1.0
import "./wdl_structs.wdl"
import "pre_process/pre_process.wdl" as preprocess
import "pre_process/qc_wkf.wdl" as qc

workflow SomaticWorkflow {
    input {
        BwaReference bwaReference
        IndexedReference indexedReference
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf DbSnp
        File chromFile
        File chromLengths  # Is this different from above?
        File hsMetricsIntervals
        File randomIntervals
        Array[Array[Fastqs]]+ tumorFastqs
        # normalFastqs can be missing for some.
        Array[Array[Fastqs]] normalFastqs
        Array[String] tumorIds
        Array[String] normalIds

    }

    scatter (i in range(length(tumorFastqs))) {
        call preprocess.Preprocess as tumorPrep {
            input:
                listOfFastqPairs = tumorFastqs[i],
                sampleId = tumorIds[i],
                bwaReference = bwaReference,
                indexedReference = indexedReference,
                MillsAnd1000G = MillsAnd1000G,
                Indels = Indels,
                DbSnp = DbSnp,
                chromFile = chromFile
        }
    }

    scatter (i in range(length(normalFastqs))) {
        call preprocess.Preprocess as normalPrep {
            input:
                listOfFastqPairs = normalFastqs[i],
                sampleId = normalIds[i],
                bwaReference = bwaReference,
                indexedReference = indexedReference,
                MillsAnd1000G = MillsAnd1000G,
                Indels = Indels,
                DbSnp = DbSnp,
                chromFile = chromFile
        }

    }

    output {
        Array[Bam] tumorBams = tumorPrep.finalBam
        Array[File] tumorAlignmentSummaryMetrics = tumorPrep.alignmentSummaryMetrics
        Array[File] tumorQualityByCyclePdf = tumorPrep.qualityByCyclePdf
        Array[File] tumorBaseDistributionByCycleMetrics = tumorPrep.baseDistributionByCycleMetrics
        Array[File] tumorQualityByCycleMetrics = tumorPrep.qualityByCycleMetrics
        Array[File] tumorBaseDistributionByCyclePdf = tumorPrep.baseDistributionByCyclePdf
        Array[File] tumorQualityDistributionPdf = tumorPrep.qualityDistributionPdf
        Array[File] tumorQualityDistributionMetrics = tumorPrep.qualityDistributionMetrics
        Array[File] tumorInsertSizeHistogramPdf = tumorPrep.insertSizeHistogramPdf
        Array[File] tumorInsertSizeMetrics = tumorPrep.insertSizeMetrics
        Array[File] tumorQualityDistributionPdfPreBqsr = tumorPrep.qualityDistributionPdfPreBqsr
        Array[File] tumorQualityByCycleMetricsPreBqsr = tumorPrep.qualityByCycleMetricsPreBqsr
        Array[File] tumorQualityByCyclePdfPreBqsr = tumorPrep.qualityByCyclePdfPreBqsr
        Array[File] tumorQualityDistributionMetricsPreBqsr = tumorPrep.qualityDistributionMetricsPreBqsr
        Array[File] tumorGcBiasMetrics = tumorPrep.gcBiasMetrics
        Array[File] tumorGcBiasSummary = tumorPrep.gcBiasSummary
        Array[File] tumorGcBiasPdf = tumorPrep.gcBiasPdf
        Array[File] tumorFlagStat = tumorPrep.FlagStat
        Array[File] tumorHsMetrics = tumorPrep.HsMetrics
        Array[File] tumorHsMetricsPerTargetCoverage = tumorPrep.HsMetricsPerTargetCoverage
        Array[File] tumorHsMetricsPerTargetCoverageAutocorr = tumorPrep.HsMetricsPerTargetCoverageAutocorr
        Array[File] tumorAutocorroutput1100 = tumorPrep.autocorroutput1100
        Array[File] tumorCollectOxoGMetrics = tumorPrep.CollectOxoGMetrics
        Array[File] tumorCollectWgsMetrics = tumorPrep.CollectWgsMetrics
        Array[File] tumorBinestCov = tumorPrep.binestCov
        Array[File] tumorNormCoverageByChrPng = tumorPrep.normCoverageByChrPng
        
        # QC files for samples.
        Array[Bam] normalBams = normalPrep.finalBam
        Array[File] normalAlignmentSummaryMetrics = normalPrep.alignmentSummaryMetrics
        Array[File] normalQualityByCyclePdf = normalPrep.qualityByCyclePdf
        Array[File] normalBaseDistributionByCycleMetrics = normalPrep.baseDistributionByCycleMetrics
        Array[File] normalQualityByCycleMetrics = normalPrep.qualityByCycleMetrics
        Array[File] normalBaseDistributionByCyclePdf = normalPrep.baseDistributionByCyclePdf
        Array[File] normalQualityDistributionPdf = normalPrep.qualityDistributionPdf
        Array[File] normalQualityDistributionMetrics = normalPrep.qualityDistributionMetrics
        Array[File] normalInsertSizeHistogramPdf = normalPrep.insertSizeHistogramPdf
        Array[File] normalInsertSizeMetrics = normalPrep.insertSizeMetrics
        Array[File] normalQualityDistributionPdfPreBqsr = normalPrep.qualityDistributionPdfPreBqsr
        Array[File] normalQualityByCycleMetricsPreBqsr = normalPrep.qualityByCycleMetricsPreBqsr
        Array[File] normalQualityByCyclePdfPreBqsr = normalPrep.qualityByCyclePdfPreBqsr
        Array[File] normalQualityDistributionMetricsPreBqsr = normalPrep.qualityDistributionMetricsPreBqsr
        Array[File] normalGcBiasMetrics = normalPrep.gcBiasMetrics
        Array[File] normalGcBiasSummary = normalPrep.gcBiasSummary
        Array[File] normalGcBiasPdf = normalPrep.gcBiasPdf
        Array[File] normalFlagStat = normalPrep.FlagStat
        Array[File] normalHsMetrics = normalPrep.HsMetrics
        Array[File] normalHsMetricsPerTargetCoverage = normalPrep.HsMetricsPerTargetCoverage
        Array[File] normalHsMetricsPerTargetCoverageAutocorr = normalPrep.HsMetricsPerTargetCoverageAutocorr
        Array[File] normalAutocorroutput1100 = normalPrep.autocorroutput1100
        Array[File] normalCollectOxoGMetrics = normalPrep.CollectOxoGMetrics
        Array[File] normalCollectWgsMetrics = normalPrep.CollectWgsMetrics
        Array[File] normalBinestCov = normalPrep.binestCov
        Array[File] normaNormCoverageByChrPng = normalPrep.normCoverageByChrPng
    }
}
