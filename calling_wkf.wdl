version 1.0

import "calling/gridss_wkf.wdl" as gridss
import "calling/bicseq2_wkf.wdl" as bicseq2
import "calling/mutect2_wkf.wdl" as mutect2
import "calling/strelka2_wkf.wdl" as strelka2
import "calling/manta_wkf.wdl" as manta
import "calling/lancet_wkf.wdl" as lancet
import "calling/gridss_wkf.wdl" as gridss
import "calling/bicseq2_wkf.wdl" as bicseq2

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

import "wdl_structs.wdl"

workflow Calling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[PairInfo]+ pairInfos
        # strelka2
        File strelkaJsonLog
        File configureStrelkaSomaticWorkflow
        #   mutect2
        File mutectJsonLog
        Array[String]+ listOfChroms
        IndexedReference referenceFa
        #   Manta
        IndexedTable callRegions
        File mantaJsonLog
        File mutectJsonLogFilter
        BwaReference bwaReference
        #   Lancet
        Map[String, File] chromBedsWgs
        File lancetJsonLog
        #   BicSeq2
        Array[String]+ listOfChromsFull
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
        Boolean highMem = false
    }
    scatter(pairInfo in pairInfos) {
        call gridss.Gridss {
            input:
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                pairName = pairInfo.pairId,
                bwaReference = bwaReference,
                gridssAdditionalReference = gridssAdditionalReference,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                bsGenome = bsGenome,
                ponTarGz = ponTarGz,
                highMem = highMem
        }

        call bicseq2.BicSeq2 {
            input:
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                readLength = readLength,
                coordReadLength = coordReadLength,
                uniqCoords = uniqCoords,
                bicseq2ConfigFile = bicseq2ConfigFile,
                bicseq2SegConfigFile = bicseq2SegConfigFile,
                chromFastas = chromFastas,
                listOfChromsFull = listOfChromsFull,
                pairName = pairInfo.pairId,
                referenceFa = referenceFa,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                tumorMedianInsertSize = tumorMedianInsertSize,
                normalMedianInsertSize = normalMedianInsertSize,
                lambda = lambda
        }

        call mutect2.Mutect2 {
            input:
                mutectJsonLog = mutectJsonLog,
                mutectJsonLogFilter = mutectJsonLogFilter,
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                listOfChroms = listOfChroms,
                pairName = pairInfo.pairId,
                referenceFa = referenceFa,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                highMem = highMem
        }

        call manta.Manta {
            input:
                mantaJsonLog = mantaJsonLog,
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                callRegions = callRegions,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                highMem = highMem
        }

        call strelka2.Strelka2 {
            input:
                strelkaJsonLog = strelkaJsonLog,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                callRegions = callRegions,
                candidateSmallIndels = Manta.candidateSmallIndels,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }

        call lancet.Lancet {
            input:
                lancetJsonLog = lancetJsonLog,
                tumor = pairInfo.tumorId,
                normal = pairInfo.normalId,
                listOfChroms = listOfChroms,
                chromBedsWgs = chromBedsWgs,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
    }

    output {
        # Gridss
        Array[IndexedVcf] gridssVcf = Gridss.gridssVcf
        # Bicseq2
        Array[File] bicseq2Png = BicSeq2.bicseq2Png
        Array[File] bicseq2 = BicSeq2.bicseq2
        # Mutect2
        Array[File] mutect2 = Mutect2.mutect2
        Array[File] mutect2Unfiltered = Mutect2.mutect2_unfiltered
        # Manta
        Array[IndexedVcf] candidateSmallIndels = Manta.candidateSmallIndels
        Array[IndexedVcf] diploidSV = Manta.diploidSV
        Array[IndexedVcf] somaticSV = Manta.somaticSV
        Array[IndexedVcf] candidateSV = Manta.candidateSV
        Array[File] unfilteredMantaSV = Manta.unfilteredMantaSV
        Array[File] filteredMantaSV = Manta.filteredMantaSV
        # Strelka2
        Array[IndexedVcf] strelka2Snvs = Strelka2.strelka2Snvs
        Array[IndexedVcf] strelka2Indels = Strelka2.strelka2Indels
        Array[File] strelka2Snv = Strelka2.strelka2Snv
        Array[File] strelka2Indel = Strelka2.strelka2Indel
        # Lancet
        Array[File] lancet = Lancet.lancet
    }
}
