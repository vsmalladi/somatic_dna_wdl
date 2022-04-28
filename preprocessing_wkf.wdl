version 1.0

import "./wdl_structs.wdl"
import "pre_process/pre_process_wkf.wdl" as preProcess
import "pre_process/qc_wkf.wdl" as qc

# breaking full pipeline to just handle the preprocessing steps to test minimap
# input is a list of samples and there fastq information, no pairing info needed

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================


# for wdl version 1.0

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
    }

    command {
        python /get_index.py \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }

    output {
        Int index = read_int(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-public/workflow_utils@sha256:b7269061a4620c6565566cbeaf61b1a58d49d26c382fa12f05f41b0e5f2e4807"
    }
}


workflow PreprocessWrapper {
    input {
        Boolean external = false

        BwaReference bwaReference
        IndexedReference referenceFa
        File adaptersFa
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        IndexedVcf gnomadBiallelic
        File bqsrCallRegions
        File chromLengths
        File hsMetricsIntervals
        File randomIntervals
        Array[sampleInfo]+ sampleInfos

        Boolean trim = true
    }

    scatter (sampleInfoObj in sampleInfos) {
        call preProcess.Preprocess {
            input:
                external = external,
                listOfFastqPairs = sampleInfoObj.listOfFastqPairs,
                trim = trim,
                adaptersFa = adaptersFa,
                sampleId = sampleInfoObj.sampleId,
                bwaReference = bwaReference,
                referenceFa = referenceFa,
                MillsAnd1000G = MillsAnd1000G,
                gnomadBiallelic = gnomadBiallelic,
                hsMetricsIntervals = hsMetricsIntervals,
                callRegions = bqsrCallRegions,
                randomIntervals = randomIntervals,
                Indels = Indels,
                dbsnp = dbsnp,
                chromLengths = chromLengths
        }

        String sampleIds = sampleInfoObj.sampleId

    }


    PreprocessingOutput workflowOutput = object {
                                             # Bams
                                             finalBams: Preprocess.finalBam,
                                             # QC
                                             alignmentSummaryMetrics: Preprocess.alignmentSummaryMetrics,
                                             qualityByCyclePdf: Preprocess.qualityByCyclePdf,
                                             baseDistributionByCycleMetrics: Preprocess.baseDistributionByCycleMetrics,
                                             qualityByCycleMetrics: Preprocess.qualityByCycleMetrics,
                                             baseDistributionByCyclePdf: Preprocess.baseDistributionByCyclePdf,
                                             qualityDistributionPdf: Preprocess.qualityDistributionPdf,
                                             qualityDistributionMetrics: Preprocess.qualityDistributionMetrics,
                                             insertSizeHistogramPdf: Preprocess.insertSizeHistogramPdf,
                                             insertSizeMetrics: Preprocess.insertSizeMetrics,
                                             gcBiasMetrics: Preprocess.gcBiasMetrics,
                                             gcBiasSummary: Preprocess.gcBiasSummary,
                                             gcBiasPdf: Preprocess.gcBiasPdf,
                                             flagStat: Preprocess.flagStat,
                                             hsMetrics: Preprocess.hsMetrics,
                                             hsMetricsPerTargetCoverage: Preprocess.hsMetricsPerTargetCoverage,
                                             hsMetricsPerTargetCoverageAutocorr: Preprocess.hsMetricsPerTargetCoverageAutocorr,
                                             autocorroutput1100: Preprocess.autocorroutput1100,
                                             collectOxoGMetrics: Preprocess.collectOxoGMetrics,
                                             collectWgsMetrics: Preprocess.collectWgsMetrics,
                                             binestCov: Preprocess.binestCov,
                                             normCoverageByChrPng: Preprocess.normCoverageByChrPng,
                                             # Dedup metrics,
                                             collectWgsMetricsPreBqsr: Preprocess.collectWgsMetricsPreBqsr,
                                             qualityDistributionPdfPreBqsr: Preprocess.qualityDistributionPdfPreBqsr,
                                             qualityByCycleMetricsPreBqsr: Preprocess.qualityByCycleMetricsPreBqsr,
                                             qualityByCyclePdfPreBqsr: Preprocess.qualityByCyclePdfPreBqsr,
                                             qualityDistributionMetricsPreBqsr: Preprocess.qualityDistributionMetricsPreBqsr,
                                         }

    output {
        PreprocessingOutput finalOutput = workflowOutput
    }

}


