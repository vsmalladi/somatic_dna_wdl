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
import "alignment_analysis/fastngsadmix_wkf.wdl" as fastNgsAdmix
import "tasks/utils.wdl" as utils

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    James Roche (jroche@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#    Minita Shah
#
# ================== /COPYRIGHT ===============================================



workflow SomaticDNA {
    input {
        Boolean external = false

        BwaReference bwaReference
        BwaMem2Reference bwamem2Reference
        IndexedReference referenceFa
        File adaptersFa
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        File bqsrCallRegions
        File chromLengths
        File hsMetricsIntervals
        File randomIntervals
        Array[SampleInfo]+ normalSampleInfos
        Array[SampleInfo]+ sampleInfos
        Array[PairRelationship]+ listOfPairRelationships

        Boolean trim = true
        Boolean production = true

        # For Tumor-Normal QC
        File markerBedFile
        File markerTxtFile

        Boolean bypassQcCheck = false

        # calling
        Array[String]+ listOfChromsFull
        Array[String]+ listOfChroms
        IndexedTable callRegions
        Map[String, File] chromBedsWgs
        Map[String, File] chromBeds
        File intervalList
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
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

        #fastNgsAdmix
        File fastNgsAdmixChroms

        File fastNgsAdmixContinentalSites
        File fastNgsAdmixContinentalSitesBin
        File fastNgsAdmixContinentalSitesIdx
        File fastNgsAdmixContinentalRef
        File fastNgsAdmixContinentalNind

        File fastNgsAdmixPopulationSites
        File fastNgsAdmixPopulationSitesBin
        File fastNgsAdmixPopulationSitesIdx
        File fastNgsAdmixPopulationRef
        File fastNgsAdmixPopulationNind

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

        # Specific mem values in GB
        Int annotateBicSeq2CnvMem = 36
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
                markerBedFile = markerBedFile,
                hsMetricsIntervals = hsMetricsIntervals,
                callRegions = bqsrCallRegions,
                randomIntervals = randomIntervals,
                Indels = Indels,
                dbsnp = dbsnp,
                chromLengths = chromLengths
        }

        # for wdl version 1.0
        String sampleIds = sampleInfoObj.sampleAnalysisId

        # for wdl version 1.1
        # Pair[String, Bam] bamPairs = (sampleInfo.sampleId, Preprocess.finalBam)

    }


    scatter (sampleInfoObj in normalSampleInfos) {
        String normalSampleIds = sampleInfoObj.sampleAnalysisId

        call utils.GetIndex as germlineRunGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = sampleInfoObj.sampleAnalysisId
        }

        Boolean skipCoverageCheck = select_first([sampleInfoObj.skipCoverageCheck, false])
        if (!bypassQcCheck && !skipCoverageCheck) {
            call BamQcCheck {
                input:
                    wgsMetricsFile = Preprocess.collectWgsMetrics[germlineRunGetIndex.index],
                    expectedCoverage = sampleInfoObj.expectedCoverage
            }
        }

        Boolean coveragePass = select_first([BamQcCheck.coveragePass, false])

        if (bypassQcCheck || skipCoverageCheck || coveragePass ) {
            call kourami.Kourami {
                input:
                    sampleId = sampleInfoObj.sampleAnalysisId,
                    kouramiReference = kouramiReference,
                    finalBam = Preprocess.finalBam[germlineRunGetIndex.index],
                    kouramiFastaGem1Index = kouramiFastaGem1Index,
                    referenceFa = referenceFa
            }

            call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixContinental{
                input:
                    normalFinalBam = Preprocess.finalBam[germlineRunGetIndex.index],
                    fastNgsAdmixSites = fastNgsAdmixContinentalSites,
                    fastNgsAdmixSitesBin = fastNgsAdmixContinentalSitesBin,
                    fastNgsAdmixSitesIdx = fastNgsAdmixContinentalSitesIdx,
                    fastNgsAdmixChroms = fastNgsAdmixChroms,
                    fastNgsAdmixRef = fastNgsAdmixContinentalRef,
                    fastNgsAdmixNind = fastNgsAdmixContinentalNind,
                    outprefix = normalSampleIds + "_continental"
            }

            call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixPopulation{
                input:
                    normalFinalBam = Preprocess.finalBam[germlineRunGetIndex.index],
                    fastNgsAdmixSites = fastNgsAdmixPopulationSites,
                    fastNgsAdmixSitesBin = fastNgsAdmixPopulationSitesBin,
                    fastNgsAdmixSitesIdx = fastNgsAdmixPopulationSitesIdx,
                    fastNgsAdmixChroms = fastNgsAdmixChroms,
                    fastNgsAdmixRef = fastNgsAdmixPopulationRef,
                    fastNgsAdmixNind = fastNgsAdmixPopulationNind,
                    outprefix = normalSampleIds + "_population"
            }

            call germline.Germline {
                input:
                    finalBam = Preprocess.finalBam[germlineRunGetIndex.index],
                    normal = sampleInfoObj.listOfFastqPairs[0].clientSampleId, # SM tag.
                    outputPrefix = sampleInfoObj.sampleAnalysisId,
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
                    normal = sampleInfoObj.sampleAnalysisId,
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
                    haplotypecallerAnnotatedVcfPath = "~{sampleInfoObj.sampleAnalysisId}.haplotypecaller.gatk.annotated.unfiltered.vcf",
                    production = production,
                    referenceFa = referenceFa,
                    normal = sampleInfoObj.sampleAnalysisId,
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
        call utils.GetIndex as tumorBafGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.tumorPrefix
        }

        call utils.GetIndex as normalBafGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.normalPrefix
        }
        call utils.GetIndex as germlineBafGetIndex {
            input:
                sampleIds = normalSampleIds,
                sampleId = pairRelationship.normalPrefix
        }

        if ( size(unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf[germlineBafGetIndex.index]) > 0 ) {
            call baf.Baf {
                input:
                    referenceFa = referenceFa,
                    pairName = pairRelationship.pairId,
                    sampleId = pairRelationship.normalId,
                    tumorFinalBam = Preprocess.finalBam[tumorBafGetIndex.index],
                    normalFinalBam = Preprocess.finalBam[normalBafGetIndex.index],
                    germlineVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf[germlineBafGetIndex.index],
                    listOfChroms = listOfChroms
                }
        }
    }

    scatter (pairRelationship in listOfPairRelationships) {
        # for wdl version 1.0
        call utils.GetIndex as tumorGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.tumorPrefix
        }

        call utils.GetIndex as normalGetIndex {
            input:
                sampleIds = sampleIds,
                sampleId = pairRelationship.normalPrefix
        }

        PairInfo pairInfoObject = object {
            pairId : pairRelationship.pairId,
            tumorFinalBam : Preprocess.finalBam[tumorGetIndex.index],
            normalFinalBam : Preprocess.finalBam[normalGetIndex.index],
            tumorId : pairRelationship.tumorId,    # SM tag
            normalId : pairRelationship.normalId   # SM tag
        }

        call conpair.Conpair {
            input:
                tumorPileupsConpair = Preprocess.pileupsConpair[tumorGetIndex.index],
                normalPileupsConpair = Preprocess.pileupsConpair[normalGetIndex.index],
                tumor = pairRelationship.tumorPrefix,
                normal = pairRelationship.normalPrefix,
                pairName = pairRelationship.pairId,
                markerTxtFile = markerTxtFile
        }

        if (!bypassQcCheck) {
            Boolean tumorSkipCoverageCheck = select_first([sampleInfos[tumorGetIndex.index].skipCoverageCheck, false])
            Boolean normalSkipCoverageCheck = select_first([sampleInfos[normalGetIndex.index].skipCoverageCheck, false])
            call SomaticQcCheck {
                input:
                    pairId = pairRelationship.pairId,
                    tumorWgsMetricsFile = Preprocess.collectWgsMetrics[tumorGetIndex.index],
                    tumorExpectedCoverage = sampleInfos[tumorGetIndex.index].expectedCoverage,
                    tumorSkipCoverageCheck = tumorSkipCoverageCheck,
                    normalWgsMetricsFile = Preprocess.collectWgsMetrics[normalGetIndex.index],
                    normalExpectedCoverage = sampleInfos[normalGetIndex.index].expectedCoverage,
                    normalSkipCoverageCheck = normalSkipCoverageCheck,
                    concordanceFile = Conpair.concordanceAll,
                    contaminationFile = Conpair.contamination
            }
        }

        Boolean somaticQcPass = select_first([SomaticQcCheck.qcPass, bypassQcCheck])
        if (somaticQcPass) {
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
                    strelkaJsonLog = strelkaJsonLog,
                    mutectJsonLog = mutectJsonLog,
                    mutectJsonLogFilter = mutectJsonLogFilter,
                    configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                    pairInfo = pairInfoObject,
                    listOfChroms = listOfChroms,
                    listOfChromsFull = listOfChromsFull,
                    referenceFa = referenceFa,
                    callRegions = callRegions,
                    bwaReference = bwaReference,
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
                    highMem = highMem,
                    library = library,
                    intervalList = intervalList,
                    chromBeds = chromBeds
            }

            call msi.Msi {
                input:
                    normal=pairRelationship.normalId,
                    pairName=pairRelationship.pairId,
                    mantisBed=mantisBed,
                    intervalListBed=intervalListBed,
                    referenceFa=referenceFa,
                    tumorFinalBam=Preprocess.finalBam[tumorGetIndex.index],
                    normalFinalBam=Preprocess.finalBam[normalGetIndex.index]
            }

            PreMergedPairVcfInfo preMergedPairVcfInfo = object {
                pairId : pairRelationship.pairId,
                filteredMantaSV : Calling.filteredMantaSV,
                strelka2Snv : Calling.strelka2Snv,
                strelka2Indel : Calling.strelka2Indel,
                mutect2 : Calling.mutect2,
                lancet : Calling.lancet,
                tumor : pairRelationship.tumorId,
                normal : pairRelationship.normalId,
                tumorFinalBam : Preprocess.finalBam[tumorGetIndex.index],
                normalFinalBam : Preprocess.finalBam[normalGetIndex.index]

            }

            PairRawVcfInfo pairRawVcfInfo = object {
                pairId : pairRelationship.pairId,
                filteredMantaSV : Calling.filteredMantaSV,
                strelka2Snv : Calling.strelka2Snv,
                strelka2Indel : Calling.strelka2Indel,
                mutect2 : Calling.mutect2,
                lancet : Calling.lancet,
                gridssVcf : Calling.gridssVcf,
                bicseq2Png : Calling.bicseq2Png,
                bicseq2 : Calling.bicseq2,
                tumor : pairRelationship.tumorId,
                normal : pairRelationship.normalId,
                tumorFinalBam : Preprocess.finalBam[tumorGetIndex.index],
                normalFinalBam : Preprocess.finalBam[normalGetIndex.index]

            }

            if (library == 'WGS') {
                call mergeVcf.MergeVcf as wgsMergeVcf {
                    input:
                        external = external,
                        preMergedPairVcfInfo = preMergedPairVcfInfo,
                        referenceFa = referenceFa,
                        listOfChroms = listOfChroms,
                        intervalListBed = intervalListBed,
                        ponFile = ponWGSFile,
                        germFile = germFile,
                        library = library

                }
            }

            if (library == 'Exome') {
                call mergeVcf.MergeVcf as exomeMergeVcf {
                    input:
                        external = external,
                        preMergedPairVcfInfo = preMergedPairVcfInfo,
                        referenceFa = referenceFa,
                        listOfChroms = listOfChroms,
                        intervalListBed = intervalListBed,
                        ponFile = ponExomeFile,
                        germFile = germFile,
                        library = library

                }
            }

            File mergedVcf = select_first([wgsMergeVcf.mergedVcf, exomeMergeVcf.mergedVcf])

            call annotate.Annotate {
                input:
                    unannotatedVcf = mergedVcf,
                    referenceFa = referenceFa,
                    production = production,
                    tumor = pairRelationship.tumorId,
                    normal = pairRelationship.normalId,
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
                    gridssVcf=pairRawVcfInfo.gridssVcf,
                    vepGenomeBuild=vepGenomeBuild,
                    gap=gap,
                    dgvBedpe=dgvBedpe,
                    thousandGVcf=thousandGVcf,
                    svPon=svPon,
                    cosmicBedPe=cosmicBedPe,
                    annotateBicSeq2CnvMem = annotateBicSeq2CnvMem
          }

          call deconstructSigs.DeconstructSig {
              input:
                    pairId = pairRelationship.pairId,
                    mainVcf = mergedVcf,
                    cosmicSigs = cosmicSigs,
                    vepGenomeBuild = vepGenomeBuild
          }

          FinalVcfPairInfo finalPairInfos = object {
                pairId : pairRelationship.pairId,
                tumor : pairRelationship.tumorId,
                normal : pairRelationship.normalId,
                mainVcf : Annotate.pairVcfInfo.mainVcf,
                supplementalVcf : Annotate.pairVcfInfo.supplementalVcf,
                vcfAnnotatedTxt : Annotate.pairVcfInfo.vcfAnnotatedTxt,
                maf : Annotate.pairVcfInfo.maf,
                filteredMantaSV : Calling.filteredMantaSV,
                strelka2Snv : Calling.strelka2Snv,
                strelka2Indel : Calling.strelka2Indel,
                mutect2 : Calling.mutect2,
                lancet : Calling.lancet,
                gridssVcf : Calling.gridssVcf,
                bicseq2Png : Calling.bicseq2Png,
                bicseq2 : Calling.bicseq2,
                cnvAnnotatedFinalBed : AnnotateCnvSv.cnvAnnotatedFinalBed,
                cnvAnnotatedSupplementalBed : AnnotateCnvSv.cnvAnnotatedSupplementalBed,
                svFinalBedPe : AnnotateCnvSv.svFinalBedPe,
                svHighConfidenceFinalBedPe : AnnotateCnvSv.svHighConfidenceFinalBedPe,
                svSupplementalBedPe : AnnotateCnvSv.svSupplementalBedPe,
                svHighConfidenceSupplementalBedPe : AnnotateCnvSv.svHighConfidenceSupplementalBedPe,

        }

      }
    }

    SomaticWorkflowOutput workflowOutput = object {
        # alignment and calling results (calling results may not exist if qc failed)
        # SNV INDELs CNV SV and BAM output
        finalPairInfo: finalPairInfos,

        # MSI
        mantisWxsKmerCountsFinal: Msi.mantisWxsKmerCountsFinal,
        mantisWxsKmerCountsFiltered: Msi.mantisWxsKmerCountsFiltered,
        mantisExomeTxt: Msi.mantisExomeTxt,
        mantisStatusFinal: Msi.mantisStatusFinal,
        # SIGs
        deconstructSigSigs: DeconstructSig.sigs,
        deconstructSigCounts: DeconstructSig.counts,
        deconstructSigSigInput: DeconstructSig.sigInput,
        deconstructSigReconstructed: DeconstructSig.reconstructed,
        deconstructSigDiff: DeconstructSig.diff,

        #ancestry
        beagleFileContinental: fastNgsAdmixContinental.beagleFile,
        fastNgsAdmixQoptContinental: fastNgsAdmixContinental.fastNgsAdmixQopt,
        beagleFilePopulation: fastNgsAdmixPopulation.beagleFile,
        fastNgsAdmixQoptPopulation: fastNgsAdmixPopulation.fastNgsAdmixQopt,

        # Bams/Crams
        finalBams: Preprocess.finalBam,
        finalCrams: Preprocess.finalCram,

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
        dedupLog: Preprocess.dedupLog,

        # Conpair
        pileupsConpair: Preprocess.pileupsConpair,
        concordanceAll: Conpair.concordanceAll,
        concordanceHomoz: Conpair.concordanceHomoz,
        contamination: Conpair.contamination,

        # Germline
        kouramiResult: Kourami.result,
        haplotypecallerVcf: Germline.haplotypecallerVcf,
        haplotypecallerFinalFiltered: Germline.haplotypecallerFinalFiltered,
        filteredHaplotypecallerAnnotatedVcf: filteredGermlineAnnotate.haplotypecallerAnnotatedVcf,
        haplotypecallerAnnotatedVcf: unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf,
        alleleCountsTxt: Baf.alleleCountsTxt
    }

    output {
       SomaticWorkflowOutput finalOutput = workflowOutput
       Array[Boolean] allQcPass = somaticQcPass
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
        docker: "gcr.io/nygc-public/workflow_utils@sha256:40fa18ac3f9d9f3b9f037ec091cb0c2c26ad6c7cb5c32fb16c1c0cf2a5c9caea"
    }
}

task SomaticQcCheck {
    # Check coverage, contamination and concordance in one go.
    input {
        String pairId
        File tumorWgsMetricsFile
        File normalWgsMetricsFile
        Float tumorExpectedCoverage
        Float normalExpectedCoverage
        File concordanceFile
        File contaminationFile
        Float minConcordance = 95.0
        Float maxContamination = 0.99
        Boolean tumorSkipCoverageCheck = false
        Boolean normalSkipCoverageCheck = false
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
           --max_contamination ~{maxContamination} \
           ${if tumorSkipCoverageCheck then "--skip_tumor_coverage" else " "} \
           ${if normalSkipCoverageCheck then "--skip_normal_coverage" else " "}

        mv qc_summary.txt ~{pairId}.qc_summary.txt

    }

    output {
        Boolean qcPass = read_boolean(stdout())
        File qcCheckReport = "~{pairId}.qc_summary.txt"
    }

    runtime {
        docker: "gcr.io/nygc-public/workflow_utils@sha256:40fa18ac3f9d9f3b9f037ec091cb0c2c26ad6c7cb5c32fb16c1c0cf2a5c9caea"
    }
}
