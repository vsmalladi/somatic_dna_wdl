version 1.0


import "qc/qc.wdl" as qc
import "./wdl_structs.wdl"

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
        #File genomeTemplates

        Int threads
        Int memoryGb
        String gatkDockerImage
        String MultipleMetricsBase
        String binestDockerImage
        String rDockerImage
        #String perlDockerImage
    }

    call qc.MultipleMetrics {
        input:
            MultipleMetricsBase = MultipleMetricsBase,
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }

    call qc.MultipleMetricsPreBqsr {
        input:
            MultipleMetricsBasePreBqsrBasename = MultipleMetricsBasePreBqsrBasename,
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }

    call qc.CollectGcBiasMetrics {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }

    call qc.Flagstat {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }

    call qc.HsMetrics {
        input:
            referenceFa = referenceFa,
            hsMetricsIntervals = hsMetricsIntervals,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = 40,
            threads = threads,
            dockerImage = gatkDockerImage
    }

    # call qc.FormatHsMetrics {
    #     input:
    #         HsMetricsPerTargetCoverage = HsMetrics.HsMetricsPerTargetCoverage,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    #         dockerImage = perlDockerImage
    # }

    # call qc.Autocorrelations {
    #     input:
    #         HsMetricsPerTargetCoverageAutocorr = FormatHsMetrics.HsMetricsPerTargetCoverageAutocorr,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    #         dockerImage = rDockerImage
    # }

    call qc.CollectOxoGMetricsWgs {
        input:
            referenceFa = referenceFa,
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }

    # call qc.CollectWgsMetricsWgsDecoy {
    #     input:
    #         referenceFa = referenceFa,
    #         randomIntervals = randomIntervals,
    #         finalBam = finalBam,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    #         dockerImage = gatkDockerImage
    # }

    call qc.Binest {
        input:
            finalBam = finalBam,
            sampleId = sampleId,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = binestDockerImage
    }

    # call qc.PlotBinCov {
    #     input:
    #         genomeTemplates = genomeTemplates,
    #         sampleId = sampleId,
    #         memoryGb = memoryGb,
    #         threads = threads,
    #         dockerImage = rDockerImage
    #}
}
