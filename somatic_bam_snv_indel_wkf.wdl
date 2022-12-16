version 1.0

import "./wdl_structs.wdl"
import "calling/snv_indel_wkf.wdl" as calling
import "alignment_analysis/msi_wkf.wdl" as msi
import "calling/calling.wdl" as callingTasks
import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "pre_process/qc.wdl" as qc
import "annotate/annotate_wkf.wdl" as annotate
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


# for wdl version 1.0


workflow SomaticBamWorkflow {
    input {
        Boolean external = false

        BwaReference bwaReference
        IndexedReference referenceFa

        Boolean production = true

        Array[PairInfo]+ pairInfos

        # For Tumor-Normal QC
        File markerBedFile
        File markerTxtFile

        # calling
        Array[String]+ listOfChromsFull
        Array[String]+ listOfChroms
        IndexedTable callRegions
        Map[String, File] chromBedsWgs
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow


        # merge callers
        File intervalListBed

        String library
        File ponWGSFile
        File ponExomeFile
        IndexedVcf gnomadBiallelic

        # mantis
        File mantisBed

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

        IndexedVcf germFile

        # signatures
        File cosmicSigs

        Boolean highMem = false
    }

    scatter(pairInfo in pairInfos) {

        call calling.Calling {
            input:
                mantaJsonLog = mantaJsonLog,
                lancetJsonLog = lancetJsonLog,
                mutectJsonLog = mutectJsonLog,
                mutectJsonLogFilter = mutectJsonLogFilter,
                strelkaJsonLog = strelkaJsonLog,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                pairInfo = pairInfo,
                listOfChroms = listOfChroms,
                listOfChromsFull = listOfChromsFull,
                referenceFa = referenceFa,
                callRegions = callRegions,
                bwaReference = bwaReference,
                chromBedsWgs = chromBedsWgs,
                highMem = highMem
        }

        call msi.Msi {
            input:
                normal = pairInfo.normalId,
                pairName = pairInfo.pairId,
                mantisBed = mantisBed,
                intervalListBed = intervalListBed,
                referenceFa = referenceFa,
                tumorFinalBam = pairInfo.tumorFinalBam,
                normalFinalBam = pairInfo.normalFinalBam
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
            tumorFinalBam : pairInfo.tumorFinalBam,
            normalFinalBam : pairInfo.normalFinalBam

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

        File runMainVcf = Annotate.pairVcfInfo.mainVcf
        File runSupplementalVcf = Annotate.pairVcfInfo.supplementalVcf
        File runVcfAnnotatedTxt = Annotate.pairVcfInfo.vcfAnnotatedTxt

        call deconstructSigs.DeconstructSig {
            input:
                pairId = pairInfo.pairId,
                mainVcf = mergedVcf,
                cosmicSigs = cosmicSigs,
                vepGenomeBuild = vepGenomeBuild
        }
   }

    output {

        # CNV SV output and SNV INDELs
        Array[File] mainVcf = runMainVcf
        Array[File] supplementalVcf = runSupplementalVcf
        Array[File] vcfAnnotatedTxt = runVcfAnnotatedTxt
        Array[File] filteredMantaSV = Calling.filteredMantaSV
        Array[File] strelka2Snv = Calling.strelka2Snv
        Array[File] strelka2Indel = Calling.strelka2Indel
        Array[File] mutect2 = Calling.mutect2
        Array[File] lancet = Calling.lancet
        # MSI
        Array[File] mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        Array[File] mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        Array[File] mantisExomeTxt = Msi.mantisExomeTxt
        Array[File] mantisStatusFinal = Msi.mantisStatusFinal
        # sigs
        Array[File] sigs = DeconstructSig.sigs
        Array[File] counts = DeconstructSig.counts
        Array[File] sig_input = DeconstructSig.sigInput
        Array[File] reconstructed = DeconstructSig.reconstructed
        Array[File] diff = DeconstructSig.diff

    }
}
