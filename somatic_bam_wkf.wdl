version 1.0

import "./wdl_structs.wdl"
import "calling/calling_wkf.wdl" as calling
import "calling/calling.wdl" as callingTasks
import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "alignment_analysis/kourami_wfk.wdl" as kourami
import "alignment_analysis/msi_wkf.wdl" as msi
import "alignment_analysis/fastngsadmix_wkf.wdl" as fastNgsAdmix
import "pre_process/conpair_wkf.wdl" as conpair
import "pre_process/qc.wdl" as qc
import "annotate/annotate_wkf.wdl" as annotate
import "annotate/annotate_cnv_sv_wkf.wdl" as annotate_cnv_sv
import "germline/germline_wkf.wdl" as germline
import "annotate/germline_annotate_wkf.wdl" as germlineAnnotate
import "baf/baf_wkf.wdl" as baf
import "variant_analysis/deconstruct_sigs_wkf.wdl" as deconstructSigs
import "tasks/bam_cram_conversion.wdl" as cramConversion
import "tasks/reheader_bam_wkf.wdl" as reheaderBam

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


workflow SomaticBamWorkflow {
    input {
        Boolean local = false

        BwaReference bwaReference
        IndexedReference referenceFa

        Boolean external = false
        Boolean production = true

        Array[PairInfo]+ pairInfos
        Array[SampleBamInfo]+ normalSampleBamInfos

        # For Tumor-Normal QC
        File markerBedFile
        File markerTxtFile

        # calling
        Array[String]+ listOfChromsFull
        Array[String]+ listOfChroms
        IndexedTable callRegions
        Map[String, File] chromBedsWgs
        Map[String, File] chromBeds
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow
        File intervalList

        #   BicSeq2
        Int readLength
        Int coordReadLength
        Map[Int, Map[String, File]] uniqCoords
        File bicseq2ConfigFile
        File bicseq2SegConfigFile
        Map[String, File] chromFastas
        Int tumorMedianInsertSize = 400
        Int normalMedianInsertSize = 400
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

        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG
        IndexedVcf dbsnp

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
        Boolean createCramBasedObjects = false
    }

    # need to find UNIQUE bams (don't convert if part of more than one pair)
    call cramConversion.UniqueBams as uniqueBams {
        input:
            pairInfosJson = write_json(pairInfos)
    }

    scatter(bamInfo in uniqueBams.uniqueBams) {        
        
        call reheaderBam.Reheader {
            input:
                finalBam = bamInfo.finalBam,
                sampleId = bamInfo.sampleId
        }
        
        call cramConversion.SamtoolsBamToCram as bamToCram {
            input:
                inputBam = Reheader.sampleBamMatched,
                referenceFa = referenceFa,
                sampleId = bamInfo.sampleId,
                diskSize = (ceil(size(Reheader.sampleBamMatched.bam, "GB") * 1.7)) + 20 # 0.7 is estimated cram size
        }

        call qc.ConpairPileup {
            input:
                markerBedFile = markerBedFile,
                referenceFa = referenceFa,
                finalBam = Reheader.sampleBamMatched,
                sampleId = bamInfo.sampleId,
                memoryGb = 4,
                threads = 1,
                diskSize = ceil(size(Reheader.sampleBamMatched.bam, "GB"))  + 20
       }

      String uniqueSampleIds = bamInfo.sampleId
    }
    if (createCramBasedObjects) {
        call cramConversion.UpdateCramInfos as updateCramInfo {
            input:
                pairInfosJson = write_json(pairInfos),
                normalInfosJson = write_json(normalSampleBamInfos),
                cramInfosJson = write_json(bamToCram.cramInfo)
        }
    }

    scatter (normalSampleBamInfo in normalSampleBamInfos) {
        String normalSampleIds = normalSampleBamInfo.sampleId
        
        call utils.GetIndex as normalGetIndex {
            input:
                sampleIds = uniqueSampleIds,
                sampleId = normalSampleBamInfo.sampleId
        }
        
        call kourami.Kourami {
            input:
                sampleId = normalSampleBamInfo.sampleId,
                kouramiReference = kouramiReference,
                finalBam = Reheader.sampleBamMatched[normalGetIndex.index],
                kouramiFastaGem1Index = kouramiFastaGem1Index,
                referenceFa = referenceFa
        }

        call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixContinental{
            input:
                normalFinalBam = Reheader.sampleBamMatched[normalGetIndex.index],
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
                normalFinalBam = Reheader.sampleBamMatched[normalGetIndex.index],
                fastNgsAdmixSites = fastNgsAdmixPopulationSites,
                fastNgsAdmixSitesBin = fastNgsAdmixPopulationSitesBin,
                fastNgsAdmixSitesIdx = fastNgsAdmixPopulationSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixRef = fastNgsAdmixPopulationRef,
                fastNgsAdmixNind = fastNgsAdmixPopulationNind,
                outprefix = normalSampleIds + "population"
        }

        call germline.Germline {
            input:
                finalBam = Reheader.sampleBamMatched[normalGetIndex.index],
                normal = normalSampleBamInfo.sampleId,
                outputPrefix = normalSampleBamInfo.sampleId,
                referenceFa = referenceFa,
                listOfChroms = listOfChroms,
                MillsAnd1000G = MillsAnd1000G,
                hapmap = hapmap,
                omni = omni,
                onekG = onekG,
                dbsnp = dbsnp,
                nygcAf = nygcAf,
                excludeIntervalList = excludeIntervalList,
                scatterIntervalsHcs = scatterIntervalsHcs,
                pgx = pgx,
                rwgsPgxBed = rwgsPgxBed,
                whitelist = whitelist,
                chdWhitelistVcf = chdWhitelistVcf,
                deepIntronicsVcf = deepIntronicsVcf,
                clinvarIntronicsVcf = clinvarIntronicsVcf,
                highMem = highMem
        }

        call germlineAnnotate.GermlineAnnotate as filteredGermlineAnnotate {
            input:
                unannotatedVcf = Germline.haplotypecallerFinalFiltered,
                production = production,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
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
                haplotypecallerAnnotatedVcfPath = "~{normalSampleBamInfo.sampleId}.haplotypecaller.gatk.annotated.unfiltered.vcf",
                production = production,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
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

    scatter(pairInfo in pairInfos) {
        call utils.GetIndex as germlineGetIndex {
            input:
                sampleIds = normalSampleIds,
                sampleId = pairInfo.normalId
        }
        
        call utils.GetIndex as tumorGetIndex {
            input:
                sampleIds = uniqueSampleIds,
                sampleId = pairInfo.tumorId
        }

        call baf.Baf {
            input:
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                sampleId = pairInfo.normalId,
                tumorFinalBam = Reheader.sampleBamMatched[tumorGetIndex.index],
                normalFinalBam = Reheader.sampleBamMatched[germlineGetIndex.index],
                germlineVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf[germlineGetIndex.index],
                listOfChroms = listOfChroms
        }
    }

    scatter(pairInfo in pairInfos) {
        call utils.GetIndex as tumorCallingGetIndex {
            input:
                sampleIds = uniqueSampleIds,
                sampleId = pairInfo.tumorId
        }

        call utils.GetIndex as normalCallingGetIndex {
            input:
                sampleIds = uniqueSampleIds,
                sampleId = pairInfo.normalId
        }

        # tumor insert size
        Int tumorDiskSize = ceil(size(Reheader.sampleBamMatched[tumorCallingGetIndex.index].bam, "GB")) + 30

        call qc.MultipleMetrics as tumorMultipleMetrics {
            input:
                referenceFa = referenceFa,
                finalBam = Reheader.sampleBamMatched[tumorCallingGetIndex.index],
                sampleId = pairInfo.tumorId,
                diskSize = tumorDiskSize
        }

        call callingTasks.GetInsertSize as tumorGetInsertSize {
            input:
                insertSizeMetrics = tumorMultipleMetrics.insertSizeMetrics
        }

        # normal insert size
        Int normalDiskSize = ceil(size(Reheader.sampleBamMatched[normalCallingGetIndex.index].bam, "GB")) + 30

        call qc.MultipleMetrics as normalMultipleMetrics {
            input:
                referenceFa = referenceFa,
                finalBam = Reheader.sampleBamMatched[normalCallingGetIndex.index],
                sampleId = pairInfo.normalId,
                diskSize = normalDiskSize
        }

        call callingTasks.GetInsertSize as normalGetInsertSize {
            input:
                insertSizeMetrics = normalMultipleMetrics.insertSizeMetrics
        }

        call conpair.Conpair {
            input:
                tumorPileupsConpair = ConpairPileup.pileupsConpair[tumorCallingGetIndex.index],
                normalPileupsConpair = ConpairPileup.pileupsConpair[normalCallingGetIndex.index],
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                pairName = pairInfo.pairId,
                markerTxtFile = markerTxtFile
        }
        
        PairInfo callingPairInfo = object {
                pairId : pairInfo.pairId,
                tumorFinalBam : Reheader.sampleBamMatched[tumorCallingGetIndex.index],
                normalFinalBam : Reheader.sampleBamMatched[normalCallingGetIndex.index],
                tumorId : pairInfo.tumorId,
                normalId : pairInfo.normalId
            }

        call calling.Calling {
            input:
                mantaJsonLog = mantaJsonLog,
                lancetJsonLog = lancetJsonLog,
                mutectJsonLog = mutectJsonLog,
                mutectJsonLogFilter = mutectJsonLogFilter,
                strelkaJsonLog = strelkaJsonLog,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                intervalList = intervalList,
                pairInfo = callingPairInfo,
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
                chromBeds = chromBeds,
                highMem = highMem,
                library = library
        }

        call msi.Msi {
            input:
                normal = pairInfo.normalId,
                pairName = pairInfo.pairId,
                mantisBed = mantisBed,
                intervalListBed = intervalListBed,
                referenceFa = referenceFa,
                tumorFinalBam = Reheader.sampleBamMatched[tumorCallingGetIndex.index],
                normalFinalBam = Reheader.sampleBamMatched[normalCallingGetIndex.index]
        }

        PreMergedPairVcfInfo preMergedPairVcfInfo = object {
            pairId : pairInfo.pairId,
            filteredMantaSV : Calling.filteredMantaSV,
            strelka2Snv : Calling.strelka2Snv,
            strelka2Indel : Calling.strelka2Indel,
            mutect2 : Calling.mutect2,
            lancet : Calling.lancet,
            tumor : pairInfo.tumorId,
            normal : pairInfo.normalId,
            tumorFinalBam : Reheader.sampleBamMatched[tumorCallingGetIndex.index],
            normalFinalBam : Reheader.sampleBamMatched[normalCallingGetIndex.index]

        }

        PairRawVcfInfo pairRawVcfInfo = object {
            pairId : pairInfo.pairId,
            filteredMantaSV : Calling.filteredMantaSV,
            strelka2Snv : Calling.strelka2Snv,
            strelka2Indel : Calling.strelka2Indel,
            mutect2 : Calling.mutect2,
            lancet : Calling.lancet,
            gridssVcf : Calling.gridssVcf,
            bicseq2Png : Calling.bicseq2Png,
            bicseq2 : Calling.bicseq2,
            tumor : pairInfo.tumorId,
            normal : pairInfo.normalId,
            tumorFinalBam : Reheader.sampleBamMatched[tumorCallingGetIndex.index],
            normalFinalBam : Reheader.sampleBamMatched[normalCallingGetIndex.index]

        }

        if (library == 'WGS') {
            call mergeVcf.MergeVcf as wgsMergeVcf {
                input:
                    external = external,
                    library = library,
                    preMergedPairVcfInfo = preMergedPairVcfInfo,
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
                    external = external,
                    preMergedPairVcfInfo = preMergedPairVcfInfo,
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
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                pairName = pairInfo.pairId,
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
                tumor = pairRawVcfInfo.tumor,
                normal = pairRawVcfInfo.normal,
                pairName = pairRawVcfInfo.pairId,
                listOfChroms =listOfChroms,
                bicseq2 = pairRawVcfInfo.bicseq2,
                cytoBand = cytoBand,
                dgv = dgv,
                thousandG = thousandG,
                cosmicUniqueBed = cosmicUniqueBed,
                cancerCensusBed = cancerCensusBed,
                ensemblUniqueBed = ensemblUniqueBed,
                filteredMantaSV = pairRawVcfInfo.filteredMantaSV,
                gridssVcf = pairRawVcfInfo.gridssVcf,
                vepGenomeBuild = vepGenomeBuild,
                gap = gap,
                dgvBedpe = dgvBedpe,
                thousandGVcf = thousandGVcf,
                svPon = svPon,
                cosmicBedPe = cosmicBedPe
        }

        call deconstructSigs.DeconstructSig {
            input:
                pairId = pairInfo.pairId,
                mainVcf = mergedVcf,
                cosmicSigs = cosmicSigs,
                vepGenomeBuild = vepGenomeBuild
        }

        FinalVcfPairInfo finalVcfPairInfo = object {
            pairId : pairInfo.pairId,
            tumor : pairInfo.tumorId,
            normal : pairInfo.normalId,
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
            svHighConfidenceSupplementalBedPe : AnnotateCnvSv.svHighConfidenceSupplementalBedPe
        }
   }

    output {
        # Germline
        Array[IndexedVcf] haplotypecallerFinalFiltered = Germline.haplotypecallerFinalFiltered
        Array[File] filteredHaplotypecallerAnnotatedVcf = filteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        Array[File] haplotypecallerAnnotatedVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        Array[File?] alleleCountsTxt = Baf.alleleCountsTxt
        Array[File] kouramiResult = Kourami.result
        # CNV SV output and SNV INDELs
        Array[FinalVcfPairInfo] finalVcfPairInfos = finalVcfPairInfo
        # MSI
        Array[File] mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        Array[File] mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        Array[File] mantisExomeTxt = Msi.mantisExomeTxt
        Array[File] mantisStatusFinal = Msi.mantisStatusFinal
        # ancestry
        Array[File] beagleFileContinental = fastNgsAdmixContinental.beagleFile
        Array[File] fastNgsAdmixQoptContinental = fastNgsAdmixContinental.fastNgsAdmixQopt
        Array[File] beagleFilePopulation = fastNgsAdmixPopulation.beagleFile
        Array[File] fastNgsAdmixQoptPopulation = fastNgsAdmixPopulation.fastNgsAdmixQopt
        # sigs
        Array[File] sigs = DeconstructSig.sigs
        Array[File] counts = DeconstructSig.counts
        Array[File] sig_input = DeconstructSig.sigInput
        Array[File] reconstructed = DeconstructSig.reconstructed
        Array[File] diff = DeconstructSig.diff

        # Conpair
        Array[File] concordanceAll = Conpair.concordanceAll
        Array[File] concordanceHomoz = Conpair.concordanceHomoz
        Array[File] contamination = Conpair.contamination
        Array[File] pileupsConpair = ConpairPileup.pileupsConpair

        # Cram
         Array[SampleCramInfo] crams = bamToCram.cramInfo
    }
}
