version 1.0

import "./wdl_structs.wdl"
import "pre_process/pre_process_wkf.wdl" as preProcess
import "pre_process/qc_wkf.wdl" as qc
import "calling/calling_wkf.wdl" as calling
import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "alignment_analysis/kourami_wfk.wdl" as kourami
import "alignment_analysis/msi_wkf.wdl" as msi

struct BamMapLike {
    Array[String] sampleId
    Array[Bam] bam
}

task CoerceMap {
    input {
        BamMapLike bamMapLike
        File zipMap = "gs://nygc-comp-s-fd4e-input/zip_map.sh"
    }

    command {
        set -e -o pipefail

        bash ~{zipMap}
    }

    output {
        Map[String, Bam] bamMaps = read_json(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
    }

    meta {
        doc: "temp function works with temp struct until wdl v1.1 is supported"
    }
}

workflow SomaticWorkflow {
    input {
        BwaReference bwaReference
        IndexedReference referenceFa
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        File bqsrCallRegions
        File chromLengths
        File hsMetricsIntervals
        File randomIntervals
        Array[sampleInfo]+ sampleInfos
        Array[PairRelationship]+ listOfPairRelationships

        # calling
        Array[String]+ listOfChroms
        IndexedTable callRegions
        File dbsnpIndels
        Map[String, File] chromBedsWgs

        # merge callers
        File intervalListBed

        String library
        File ponWGSFile
        File ponExomeFile
        IndexedVcf gnomadBiallelic

        IndexedVcf germFile

        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem3Index

        # mantis
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa
    }

    scatter (sampleInfoObj in sampleInfos) {
        call preProcess.Preprocess {
            input:
                listOfFastqPairs = sampleInfoObj.listOfFastqPairs,
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

        call kourami.Kourami {
            input:
                sampleId=sampleInfoObj.sampleId,
                kouramiReference=kouramiReference,
                finalBam=Preprocess.finalBam,
                kouramiFastaGem3Index=kouramiFastaGem3Index
        }

        # for wdl version 1.0
        String sampleIds = sampleInfoObj.sampleId

        # for wdl version 1.1
        # Pair[String, Bam] bamPairs = (sampleInfo.sampleId, Preprocess.finalBam)

    }

    # for wdl version 1.0
    BamMapLike bamMapLike = object {
                                    sampleId: sampleIds,
                                    bam: Preprocess.finalBam
                                    }
    # for wdl version 1.1
    # Map[String, Bam] bamMaps = as_map(bamPairs)

    call CoerceMap {
        input:
            bamMapLike = bamMapLike
    }

    scatter (pairRelationship in listOfPairRelationships) {

        pairInfo pairInfoObject = object {
            pairId : pairRelationship.pairId,
            tumorFinalBam : CoerceMap.bamMaps[pairRelationship.tumor],
            normalFinalBam : CoerceMap.bamMaps[pairRelationship.normal],
            tumor : pairRelationship.tumor,
            normal : pairRelationship.normal
        }

        call calling.Calling {
            input:
                pairInfo = pairInfoObject,
                listOfChroms = listOfChroms,
                referenceFa = referenceFa,
                callRegions = callRegions,
                bwaReference = bwaReference,
                dbsnpIndels = dbsnpIndels,
                chromBedsWgs = chromBedsWgs

        }

        call msi.Msi {
            input:
                normal=pairRelationship.normal,
                pairName=pairRelationship.pairId,
                mantisBed=mantisBed,
                intervalListBed=intervalListBed,
                referenceFa=referenceFa,
                tumorFinalBam=CoerceMap.bamMaps[pairRelationship.tumor],
                normalFinalBam=CoerceMap.bamMaps[pairRelationship.normal]
        }

        PairRawVcfInfo pairRawVcfInfo = object {
            pairId : pairRelationship.pairId,
            filteredMantaSV : Calling.filteredMantaSV,
            strelka2Snv : Calling.strelka2Snv,
            strelka2Indel : Calling.strelka2Indel,
            mutect2 : Calling.mutect2,
            lancet : Calling.lancet,
            svabaSv : Calling.svabaSv,
            svabaIndel : Calling.svabaIndel,
            tumor : pairRelationship.tumor,
            normal : pairRelationship.normal,
            tumorFinalBam : CoerceMap.bamMaps[pairRelationship.tumor],
            normalFinalBam : CoerceMap.bamMaps[pairRelationship.normal]

        }

        if (library == 'WGS') {
            call mergeVcf.MergeVcf as wgsMergeVcf {
                input:
                    pairRawVcfInfo = pairRawVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponWGSFile,
                    germFile = germFile

            }
        }

        if (library == 'Exome') {
            call mergeVcf.MergeVcf as exomeMergeVcf {
                input:
                    pairRawVcfInfo = pairRawVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponExomeFile,
                    germFile = germFile

            }
        }

        File mergedVcf = select_first([wgsMergeVcf.mergedVcf, exomeMergeVcf.mergedVcf])
    }

    output {
        Array[File] mergedVcfs = mergedVcf
        Map[String, Bam] bamMaps = CoerceMap.bamMaps
        Array[PairRawVcfInfo] pairRawVcfInfos = pairRawVcfInfo
        Array[File] kouramiResult = Kourami.result
        Array[File] mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        Array[File] mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        Array[File] mantisExomeTxt = Msi.mantisExomeTxt
        Array[File] mantisStatusFinal = Msi.mantisStatusFinal
        # QC
        Array[File] alignmentSummaryMetrics = Preprocess.alignmentSummaryMetrics
        Array[File] qualityByCyclePdf = Preprocess.qualityByCyclePdf
        Array[File] baseDistributionByCycleMetrics = Preprocess.baseDistributionByCycleMetrics
        Array[File] qualityByCycleMetrics = Preprocess.qualityByCycleMetrics
        Array[File] baseDistributionByCyclePdf = Preprocess.baseDistributionByCyclePdf
        Array[File] qualityDistributionPdf = Preprocess.qualityDistributionPdf
        Array[File] qualityDistributionMetrics = Preprocess.qualityDistributionMetrics
        Array[File] insertSizeHistogramPdf = Preprocess.insertSizeHistogramPdf
        Array[File] insertSizeMetrics = Preprocess.insertSizeMetrics
        Array[File] gcBiasMetrics = Preprocess.gcBiasMetrics
        Array[File] gcBiasSummary = Preprocess.gcBiasSummary
        Array[File] gcBiasPdf = Preprocess.gcBiasPdf
        Array[File] flagStat = Preprocess.flagStat
        Array[File] hsMetrics = Preprocess.hsMetrics
        Array[File] hsMetricsPerTargetCoverage = Preprocess.hsMetricsPerTargetCoverage
        Array[File] hsMetricsPerTargetCoverageAutocorr = Preprocess.hsMetricsPerTargetCoverageAutocorr
        Array[File] autocorroutput1100 = Preprocess.autocorroutput1100
        Array[File] collectOxoGMetrics = Preprocess.collectOxoGMetrics
        Array[File] collectWgsMetrics = Preprocess.collectWgsMetrics
        Array[File] binestCov = Preprocess.binestCov
        Array[File] normCoverageByChrPng = Preprocess.normCoverageByChrPng
        # Dedup metrics.
        Array[File] collectWgsMetricsPreBqsr = Preprocess.collectWgsMetricsPreBqsr
        Array[File] qualityDistributionPdfPreBqsr = Preprocess.qualityDistributionPdfPreBqsr
        Array[File] qualityByCycleMetricsPreBqsr = Preprocess.qualityByCycleMetricsPreBqsr
        Array[File] qualityByCyclePdfPreBqsr = Preprocess.qualityByCyclePdfPreBqsr
        Array[File] qualityDistributionMetricsPreBqsr = Preprocess.qualityDistributionMetricsPreBqsr
    }
}
