version 1.0

import "./wdl_structs.wdl"
import "pre_process/pre_process_wkf.wdl" as preProcess
import "pre_process/qc_wkf.wdl" as qc

# breaking full pipeline to just handle the preprocessing steps to test aligners
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
        docker: "gcr.io/nygc-public/workflow_utils@sha256:40fa18ac3f9d9f3b9f037ec091cb0c2c26ad6c7cb5c32fb16c1c0cf2a5c9caea"
    }
}


workflow PreprocessWrapper {
    input {
        Boolean external = false
        Boolean highMem = false

        BwaMem2Reference bwamem2Reference
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
        Array[SampleInfo]+ sampleInfos

        Boolean trim = true
    }

    scatter (sampleInfoObj in sampleInfos) {
        call preProcess.Preprocess {
            input:
                external = external,
                highMem = highMem,
                listOfFastqPairs = sampleInfoObj.listOfFastqPairs,
                trim = trim,
                adaptersFa = adaptersFa,
                sampleId = sampleInfoObj.sampleAnalysisId,
                bwamem2Reference = bwamem2Reference,
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

        String sampleIds = sampleInfoObj.sampleAnalysisId

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
                                             binestSex: Preprocess.binestSex,
                                             normCoverageByChrPng: Preprocess.normCoverageByChrPng,
                                             # Dedup metrics,
                                             collectWgsMetricsPreBqsr: Preprocess.collectWgsMetricsPreBqsr,
                                             qualityDistributionPdfPreBqsr: Preprocess.qualityDistributionPdfPreBqsr,
                                             qualityByCycleMetricsPreBqsr: Preprocess.qualityByCycleMetricsPreBqsr,
                                             qualityByCyclePdfPreBqsr: Preprocess.qualityByCyclePdfPreBqsr,
                                             qualityDistributionMetricsPreBqsr: Preprocess.qualityDistributionMetricsPreBqsr,
                                             pileupsConpair: Preprocess.pileupsConpair
                                         }

    output {
        PreprocessingOutput finalOutput = workflowOutput
    }

}
