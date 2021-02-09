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
        #File randomIntervals
        #File chromLengths

        Int threads
        Int memoryGb
    }

    call qc.MultipleMetrics {
        input:
            MultipleMetricsBase = MultipleMetricsBase,
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
    }

    call qc.MultipleMetricsPreBqsr {
        input:
            MultipleMetricsBasePreBqsrBasename = MultipleMetricsBasePreBqsrBasename,
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
    }

    call qc.CollectGcBiasMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
    }

    call qc.Flagstat {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
    }

    call qc.HsMetrics {
        input:
            referenceFa = referenceFa,
            hsMetricsIntervals = hsMetricsIntervals,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = 40,
            threads = threads,
    }

    # call qc.FormatHsMetrics {
    #     input:
    #         HsMetricsPerTargetCoverage = HsMetrics.HsMetricsPerTargetCoverage,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    # }

    # call qc.Autocorrelations {
    #     input:
    #         HsMetricsPerTargetCoverageAutocorr = FormatHsMetrics.HsMetricsPerTargetCoverageAutocorr,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    # }

    call qc.CollectOxoGMetricsWgs {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
    }

    # call qc.CollectWgsMetricsWgsDecoy {
    #     input:
    #         referenceFa = referenceFa,
    #         randomIntervals = randomIntervals,
    #         finalBam = finalBam,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    # }

    call qc.Binest {
        input:
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
    }

    # call qc.PlotBinCov {
    #     input:
    #         chromLengths = chromLengths,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    #}
}
