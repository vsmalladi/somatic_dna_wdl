version 1.0

import "./wdl_structs.wdl"
import "pre_process/pre_process_wkf.wdl" as preProcess
import "pre_process/qc_wkf.wdl" as qc
import "calling/calling_wkf.wdl" as calling
import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "alignment_analysis/kourami_wfk.wdl" as kourami
import "alignment_analysis/msi_wkf.wdl" as msi

# for wdl version 1.0

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
        File getIndex = "gs://nygc-comp-s-fd4e-input/get_index.py"
    }

    command {
        python ~{getIndex} \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }

    output {
        Int index = read_int(stdout())
    }

    runtime {
        docker: "python:3"
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
        File kouramiFastaGem1Index
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
                sampleId = sampleInfoObj.sampleId,
                kouramiReference = kouramiReference,
                finalBam = Preprocess.finalBam,
                kouramiFastaGem1Index = kouramiFastaGem1Index
        }

        # for wdl version 1.0
        String sampleIds = sampleInfoObj.sampleId

        # for wdl version 1.1
        # Pair[String, Bam] bamPairs = (sampleInfo.sampleId, Preprocess.finalBam)

    }
    
    # for wdl version 1.1
    # Map[String, Bam] bamMaps = as_map(bamPairs)

    
    scatter (pairRelationship in listOfPairRelationships) {
        # for wdl version 1.0
        call GetIndex as tumorGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.tumor 
        }
        
        call GetIndex as normalGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.normal 
        }

        pairInfo pairInfoObject = object {
            pairId : pairRelationship.pairId,
            tumorFinalBam : Preprocess.finalBam[tumorGetIndex.index],
            normalFinalBam : Preprocess.finalBam[normalGetIndex.index],
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
                tumorFinalBam=Preprocess.finalBam[tumorGetIndex.index],
                normalFinalBam=Preprocess.finalBam[normalGetIndex.index]
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
            tumorFinalBam : Preprocess.finalBam[tumorGetIndex.index],
            normalFinalBam : Preprocess.finalBam[normalGetIndex.index]

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
        Array[Bam] finalBams = Preprocess.finalBam
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
