version 1.0

import "./wdl_structs.wdl"
import "pre_process/pre_process_wkf.wdl" as preProcess
import "pre_process/qc_wkf.wdl" as qc
import "calling/calling_wkf.wdl" as calling
import "calling/calling.wdl" as callingTasks
import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "alignment_analysis/kourami_wfk.wdl" as kourami
import "alignment_analysis/msi_wkf.wdl" as msi
import "pre_process/conpair_wkf.wdl" as conpair
import "annotate/annotate_wkf.wdl" as annotate
import "annotate/annotate_cnv_sv_wkf.wdl" as annotate_cnv_sv
import "germline/germline_wkf.wdl" as germline
import "annotate/germline_annotate_wkf.wdl" as germlineAnnotate
import "baf/baf_wkf.wdl" as baf
import "variant_analysis/deconstruct_sigs_wkf.wdl" as deconstructSigs

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
        docker: "gcr.io/nygc-internal-tools/workflow_utils:2.0"
    }
}


workflow SomaticWorkflow {
    input {
        BwaReference bwaReference
        IndexedReference referenceFa
        File adaptersFa
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        File bqsrCallRegions
        File chromLengths
        File hsMetricsIntervals
        File randomIntervals
        Array[sampleInfo]+ normalSampleInfos
        Array[sampleInfo]+ sampleInfos
        Array[PairRelationship]+ listOfPairRelationships

        Boolean trim = true
        Boolean production = true
        Boolean bypassQcCheck = false

        # For Tumor-Normal QC
        File markerBedFile
        File markerTxtFile
        Boolean bypassQcCheck = false

        # calling
        Array[String]+ listOfChromsFull
        Array[String]+ listOfChroms
        IndexedTable callRegions
        File dbsnpIndels
        Map[String, File] chromBedsWgs
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File svabaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow

        #   BicSeq2
        Int readLength
        Int coordReadLength
        Map[Int, Map[String, File]] uniqCoords
        File bicseq2ConfigFile
        File bicseq2SegConfigFile
        Map[String, File] chromFastas
        Int lambda = 4

        # Gridss
        String bsGenome
        File ponTarGz
        Array[File] gridssAdditionalReference

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

        # annotation:
        String vepGenomeBuild
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding

        # Public
        File cancerResistanceMutations
        File vepCache
        File annotations
        File plugins
        String vepGenomeBuild
        IndexedReference vepFastaReference

        # NYGC-only
        IndexedVcf hgmdGene
        IndexedVcf hgmdUd10
        IndexedVcf hgmdPro
        IndexedVcf omimVcf

        # Public
        IndexedVcf chdGenesVcf
        IndexedVcf chdEvolvingGenesVcf
        IndexedVcf chdWhitelistVcf
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf masterMind

        # annotate cnv
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed

        # annotate sv
        String vepGenomeBuild
        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe


        # post annotation
        File cosmicCensus

        File ensemblEntrez

        # germline

        File excludeIntervalList
        Array[File] scatterIntervalsHcs

        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG

        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf
        
        # signatures
        File cosmicSigs
        
        Boolean highMem = false

    }

    scatter (sampleInfoObj in sampleInfos) {
        call preProcess.Preprocess {
            input:
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

        # for wdl version 1.0
        String sampleIds = sampleInfoObj.sampleId
        # for wdl version 1.1
        # Pair[String, Bam] bamPairs = (sampleInfo.sampleId, Preprocess.finalBam)

    }

    scatter (sampleInfoObj in normalSampleInfos) {
        String normalSampleIds = sampleInfoObj.sampleId

        call GetIndex as germlineRunGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = sampleInfoObj.sampleId
        }

        if (!bypassQcCheck) {
            call BamQcCheck {
                input:
                    wgsMetricsFile = Preprocess.collectWgsMetrics[germlineRunGetIndex.index],
                    expectedCoverage = sampleInfoObj.expectedCoverage
            }
        }

        Boolean coveragePass = select_first([BamQcCheck.coveragePass, false])

        if (bypassQcCheck || coveragePass ) {
            call kourami.Kourami {
                input:
                    sampleId = sampleInfoObj.sampleId,
                    kouramiReference = kouramiReference,
                    finalBam = Preprocess.finalBam[germlineRunGetIndex.index],
                    kouramiFastaGem1Index = kouramiFastaGem1Index
            }

            call germline.Germline {
                input:
                    finalBam = Preprocess.finalBam[germlineRunGetIndex.index],
                    normal = sampleInfoObj.sampleId,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    MillsAnd1000G = MillsAnd1000G,
                    hapmap = hapmap,
                    omni = omni,
                    onekG = onekG,
                    dbsnp = dbsnp,
                    nygcAf = nygcAf,
                    excludeIntervalList=excludeIntervalList,
                    scatterIntervalsHcs=scatterIntervalsHcs,
                    pgx=pgx,
                    rwgsPgxBed=rwgsPgxBed,
                    whitelist=whitelist,
                    chdWhitelistVcf=chdWhitelistVcf,
                    deepIntronicsVcf=deepIntronicsVcf,
                    clinvarIntronicsVcf=clinvarIntronicsVcf,
                    highMem=highMem
            }

            call germlineAnnotate.GermlineAnnotate as filteredGermlineAnnotate {
                input:
                    unannotatedVcf = Germline.haplotypecallerFinalFiltered,
                    production = production,
                    referenceFa = referenceFa,
                    normal = sampleInfoObj.sampleId,
                    listOfChroms = listOfChroms,
                    vepGenomeBuild = vepGenomeBuild,
                    cosmicCoding = cosmicCoding,
                    cosmicNoncoding = cosmicNoncoding,
                    # Public
                    cancerResistanceMutations = cancerResistanceMutations,
                    vepCache = vepCache,
                    annotations = annotations,
                    plugins = plugins,
                    vepFastaReference = vepFastaReference,
                    # NYGC-only
                    hgmdGene = hgmdGene,
                    hgmdUd10 = hgmdUd10,
                    hgmdPro = hgmdPro,
                    omimVcf = omimVcf,
                    # Public
                    chdGenesVcf = chdGenesVcf,
                    chdEvolvingGenesVcf = chdEvolvingGenesVcf,
                    chdWhitelistVcf = chdWhitelistVcf,
                    deepIntronicsVcf = deepIntronicsVcf,
                    clinvarIntronicsVcf = clinvarIntronicsVcf,
                    masterMind = masterMind,
                    # post annotation
                    cosmicCensus = cosmicCensus,
                    ensemblEntrez = ensemblEntrez,
                    library = library
            }

            call germlineAnnotate.GermlineAnnotate as unFilteredGermlineAnnotate {
                input:
                    unannotatedVcf = Germline.haplotypecallerVcf,
                    haplotypecallerAnnotatedVcfPath = "~{sampleInfoObj.sampleId}.haplotypecaller.gatk.v4.1.8.0.annotated.unfiltered.vcf",
                    production = production,
                    referenceFa = referenceFa,
                    normal = sampleInfoObj.sampleId,
                    listOfChroms = listOfChroms,
                    vepGenomeBuild = vepGenomeBuild,
                    cosmicCoding = cosmicCoding,
                    cosmicNoncoding = cosmicNoncoding,
                    # Public
                    cancerResistanceMutations = cancerResistanceMutations,
                    vepCache = vepCache,
                    annotations = annotations,
                    plugins = plugins,
                    vepFastaReference = vepFastaReference,
                    # NYGC-only
                    hgmdGene = hgmdGene,
                    hgmdUd10 = hgmdUd10,
                    hgmdPro = hgmdPro,
                    omimVcf = omimVcf,
                    # Public
                    chdGenesVcf = chdGenesVcf,
                    chdEvolvingGenesVcf = chdEvolvingGenesVcf,
                    chdWhitelistVcf = chdWhitelistVcf,
                    deepIntronicsVcf = deepIntronicsVcf,
                    clinvarIntronicsVcf = clinvarIntronicsVcf,
                    masterMind = masterMind,
                    # post annotation
                    cosmicCensus = cosmicCensus,
                    ensemblEntrez = ensemblEntrez,
                    library = library
            }
        }
    }

    # for wdl version 1.1
    # Map[String, Bam] bamMaps = as_map(bamPairs)

    scatter (pairRelationship in listOfPairRelationships) {
        # for wdl version 1.0
        call GetIndex as tumorBafGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.tumor
        }

        call GetIndex as normalBafGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.normal
        }
        call GetIndex as germlineBafGetIndex {
            input:
                sampleIds = normalSampleIds,
                sampleId = pairRelationship.normal
        }
        
        if ( size(unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf[germlineBafGetIndex.index]) > 0 ) {
            call baf.Baf {
                input:
                    referenceFa = referenceFa,
                    pairName = pairRelationship.pairId,
                    sampleId = pairRelationship.normal,
                    tumorFinalBam = Preprocess.finalBam[tumorBafGetIndex.index],
                    normalFinalBam = Preprocess.finalBam[normalBafGetIndex.index],
                    germlineVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf[germlineBafGetIndex.index],
                    listOfChroms = listOfChroms
                }
        }
    }

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

        call conpair.Conpair {
            input:
                finalTumorBam = Preprocess.finalBam[tumorGetIndex.index],
                finalNormalBam = Preprocess.finalBam[normalGetIndex.index],
                tumor = pairRelationship.tumor,
                normal = pairRelationship.normal,
                pairName = pairRelationship.pairId,
                referenceFa = referenceFa,
                markerBedFile = markerBedFile,
                markerTxtFile = markerTxtFile
        }

        if (!bypassQcCheck) {
            call SomaticQcCheck {
                input:
                    tumorWgsMetricsFile = Preprocess.collectWgsMetrics[tumorGetIndex.index],
                    tumorExpectedCoverage = sampleInfos[tumorGetIndex.index].expectedCoverage,
                    normalWgsMetricsFile = Preprocess.collectWgsMetrics[normalGetIndex.index],
                    normalExpectedCoverage = sampleInfos[normalGetIndex.index].expectedCoverage,
                    concordanceFile = Conpair.concordanceAll,
                    contaminationFile = Conpair.contamination
            }
            Boolean qcPass = SomaticQcCheck.qcPass
        }

        if (bypassQcCheck || select_first([SomaticQcCheck.qcPass, false])) {
            call callingTasks.GetInsertSize as tumorGetInsertSize {
                input:
                    insertSizeMetrics = Preprocess.insertSizeMetrics[tumorGetIndex.index]
            }
            
            call callingTasks.GetInsertSize as normalGetInsertSize {
                input:
                    insertSizeMetrics = Preprocess.insertSizeMetrics[normalGetIndex.index]
            }

            call calling.Calling {
                input:
                    mantaJsonLog = mantaJsonLog,
                    lancetJsonLog = lancetJsonLog,
                    mutectJsonLog = mutectJsonLog,
                    mutectJsonLogFilter = mutectJsonLogFilter,
                    svabaJsonLog = svabaJsonLog,
                    strelkaJsonLog = strelkaJsonLog,
                    configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                    pairInfo = pairInfoObject,
                    listOfChroms = listOfChroms,
                    listOfChromsFull = listOfChromsFull,
                    referenceFa = referenceFa,
                    callRegions = callRegions,
                    bwaReference = bwaReference,
                    dbsnpIndels = dbsnpIndels,
                    chromBedsWgs = chromBedsWgs,
                    readLength = readLength,
                    coordReadLength = coordReadLength,
                    uniqCoords = uniqCoords,
                    bicseq2ConfigFile = bicseq2ConfigFile,
                    bicseq2SegConfigFile = bicseq2SegConfigFile,
                    tumorMedianInsertSize = tumorGetInsertSize.insertSize,
                    normalMedianInsertSize = normalGetInsertSize.insertSize,
                    chromFastas = chromFastas,
                    bsGenome = bsGenome,
                    ponTarGz = ponTarGz,
                    gridssAdditionalReference = gridssAdditionalReference,
                    highMem = highMem
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
                gridssVcf : Calling.gridssVcf,
                bicseq2Png : Calling.bicseq2Png,
                bicseq2 : Calling.bicseq2,
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

            call annotate.Annotate {
                input:
                    unannotatedVcf = mergedVcf,
                    referenceFa = referenceFa,
                    production = production,
                    tumor = pairRelationship.tumor,
                    normal = pairRelationship.normal,
                    pairName = pairRelationship.pairId,
                    vepGenomeBuild = vepGenomeBuild,
                    cosmicCoding = cosmicCoding,
                    cosmicNoncoding = cosmicNoncoding,
                    # Public
                    cancerResistanceMutations = cancerResistanceMutations,
                    vepCache = vepCache,
                    annotations = annotations,
                    plugins = plugins,
                    vepFastaReference = vepFastaReference,
                    # NYGC-only
                    hgmdGene = hgmdGene,
                    hgmdUd10 = hgmdUd10,
                    hgmdPro = hgmdPro,
                    omimVcf = omimVcf,
                    # Public
                    chdGenesVcf = chdGenesVcf,
                    chdEvolvingGenesVcf = chdEvolvingGenesVcf,
                    chdWhitelistVcf = chdWhitelistVcf,
                    deepIntronicsVcf = deepIntronicsVcf,
                    clinvarIntronicsVcf = clinvarIntronicsVcf,
                    masterMind = masterMind,
                    # post annotation
                    cosmicCensus = cosmicCensus,
                    ensemblEntrez = ensemblEntrez,
                    library = library
            }

            call annotate_cnv_sv.AnnotateCnvSv {
                input:
                    tumor=pairRawVcfInfo.tumor,
                    normal=pairRawVcfInfo.normal,
                    pairName=pairRawVcfInfo.pairId,
                    listOfChroms=listOfChroms,
                    bicseq2=pairRawVcfInfo.bicseq2,
                    cytoBand=cytoBand,
                    dgv=dgv,
                    thousandG=thousandG,
                    cosmicUniqueBed=cosmicUniqueBed,
                    cancerCensusBed=cancerCensusBed,
                    ensemblUniqueBed=ensemblUniqueBed,

                    filteredMantaSV=pairRawVcfInfo.filteredMantaSV,
                    svabaSv=pairRawVcfInfo.svabaSv,
                    gridssVcf=pairRawVcfInfo.gridssVcf,
                    vepGenomeBuild=vepGenomeBuild,
                    gap=gap,
                    dgvBedpe=dgvBedpe,
                    thousandGVcf=thousandGVcf,
                    svPon=svPon,
                    cosmicBedPe=cosmicBedPe

          }
          
          call deconstructSigs.DeconstructSig {
              input: 
                    pairId = pairRelationship.pairId,
                    mainVcf = mergedVcf,
                    cosmicSigs = cosmicSigs,
                    vepGenomeBuild = vepGenomeBuild
          }
          
          FinalPairInfo finalPairInfos = object {
                pairId : pairRelationship.pairId,
                tumor : pairRelationship.tumor,
                normal : pairRelationship.normal,
                mainVcf : Annotate.pairVcfInfo.mainVcf,
                supplementalVcf : Annotate.pairVcfInfo.supplementalVcf,
                vcfAnnotatedTxt : Annotate.pairVcfInfo.vcfAnnotatedTxt,
                maf : Annotate.pairVcfInfo.maf,
                filteredMantaSV : Calling.filteredMantaSV,
                strelka2Snv : Calling.strelka2Snv,
                strelka2Indel : Calling.strelka2Indel,
                mutect2 : Calling.mutect2,
                lancet : Calling.lancet,
                svabaSv : Calling.svabaSv,
                svabaIndel : Calling.svabaIndel,
                gridssVcf : Calling.gridssVcf,
                bicseq2Png : Calling.bicseq2Png,
                bicseq2 : Calling.bicseq2,
                cnvAnnotatedFinalBed : AnnotateCnvSv.cnvAnnotatedFinalBed,
                cnvAnnotatedSupplementalBed : AnnotateCnvSv.cnvAnnotatedSupplementalBed,
                svFinalBedPe : AnnotateCnvSv.svFinalBedPe,
                svHighConfidenceFinalBedPe : AnnotateCnvSv.svHighConfidenceFinalBedPe,
                svSupplementalBedPe : AnnotateCnvSv.svSupplementalBedPe,
                svHighConfidenceSupplementalBedPe : AnnotateCnvSv.svHighConfidenceSupplementalBedPe,
                tumorFinalBam : Preprocess.finalBam[tumorGetIndex.index],
                normalFinalBam : Preprocess.finalBam[normalGetIndex.index]
        }
          
      }

   }

    output {
        # alignment and calling results (calling results may not exist if qc failed)
        # SNV INDELs CNV SV and BAM output
        Array[FinalPairInfo?] finalPairInfo = finalPairInfos
        
        # MSI
        Array[File?] mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        Array[File?] mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        Array[File?] mantisExomeTxt = Msi.mantisExomeTxt
        Array[File?] mantisStatusFinal = Msi.mantisStatusFinal
        # SIGs
        Array[File?] sigs = DeconstructSig.sigs
        Array[File?] counts = DeconstructSig.counts
        Array[File?] sig_input = DeconstructSig.sigInput
        Array[File?] reconstructed = DeconstructSig.reconstructed
        Array[File?] diff = DeconstructSig.diff
        
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
        # Dedup metrics
        Array[File] collectWgsMetricsPreBqsr = Preprocess.collectWgsMetricsPreBqsr
        Array[File] qualityDistributionPdfPreBqsr = Preprocess.qualityDistributionPdfPreBqsr
        Array[File] qualityByCycleMetricsPreBqsr = Preprocess.qualityByCycleMetricsPreBqsr
        Array[File] qualityByCyclePdfPreBqsr = Preprocess.qualityByCyclePdfPreBqsr
        Array[File] qualityDistributionMetricsPreBqsr = Preprocess.qualityDistributionMetricsPreBqsr

        # Conpair
        Array[File] concordanceAll = Conpair.concordanceAll
        Array[File] concordanceHomoz = Conpair.concordanceHomoz
        Array[File] contamination = Conpair.contamination

        # Germline
        Array[File?] kouramiResult = Kourami.result
        Array[IndexedVcf?] haplotypecallerVcf = Germline.haplotypecallerVcf
        Array[IndexedVcf?] haplotypecallerFinalFiltered = Germline.haplotypecallerFinalFiltered
        Array[File?] filteredHaplotypecallerAnnotatedVcf = filteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        Array[File?] haplotypecallerAnnotatedVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        Array[File?] alleleCountsTxt = Baf.alleleCountsTxt
    }
}

task BamQcCheck {
    input {
        File wgsMetricsFile
        Float expectedCoverage
    }

    command {
        python /check_wgs_coverage.py -m ~{wgsMetricsFile} -e ~{expectedCoverage}
    }

    output {
        Boolean coveragePass = read_boolean(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-internal-tools/workflow_utils:3.0"
    }
}

task SomaticQcCheck {
    # Check coverage, contamination and concordance in one go.
    input {
        File tumorWgsMetricsFile
        File normalWgsMetricsFile
        Float tumorExpectedCoverage
        Float normalExpectedCoverage
        File concordanceFile
        File contaminationFile
        Float minConcordance = 95.0
        Float maxContamination = 0.99
    }

    command {
        python /check_somatic_qc.py \
           --tumor_metrics_file ~{tumorWgsMetricsFile} \
           --tumor_expected_coverage ~{tumorExpectedCoverage} \
           --normal_metrics_file ~{normalWgsMetricsFile} \
           --normal_expected_coverage ~{normalExpectedCoverage} \
           --concordance_file ~{concordanceFile} \
           --contamination_file ~{contaminationFile} \
           --min_concordance ~{minConcordance} \
           --max_contamination ~{maxContamination}
    }

    output {
        Boolean qcPass = read_boolean(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-internal-tools/workflow_utils:3.0"
    }
}
