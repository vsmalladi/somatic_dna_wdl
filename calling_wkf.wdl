version 1.0

import "calling/calling_wkf.wdl" as calling
import "tasks/bam_cram_conversion.wdl" as cramConversion
import "tasks/reheader_bam_wkf.wdl" as reheaderBam
import "tasks/utils.wdl" as utils

import "wdl_structs.wdl"

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
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

workflow Calling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean local = false
        String library
        
        Array[PairInfo]+ pairInfos
        # strelka2
        File strelkaJsonLog
        File configureStrelkaSomaticWorkflow
        File intervalListBed
        #   mutect2
        File mutectJsonLog
        Array[String]+ listOfChroms
        Array[String]+ callerIntervals
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
        Map[String, File] chromBeds
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
        
        String uniqueSampleIds = bamInfo.sampleId
    }
    
    scatter(pairInfo in pairInfos) {
        call utils.GetIndex as normalGetIndex {
            input:
                sampleIds = uniqueSampleIds,
                sampleId = pairInfo.normalId
        }
        
        call utils.GetIndex as tumorGetIndex {
            input:
                sampleIds = uniqueSampleIds,
                sampleId = pairInfo.tumorId
        }
        
        PairInfo callingPairInfo = object {
                pairId : pairInfo.pairId,
                tumorFinalBam : Reheader.sampleBamMatched[tumorGetIndex.index],
                normalFinalBam : Reheader.sampleBamMatched[normalGetIndex.index],
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
                pairInfo = callingPairInfo,
                listOfChroms = listOfChroms,
                listOfChromsFull = listOfChromsFull,
                callerIntervals = callerIntervals,
                chromBeds = chromBeds,
                referenceFa = referenceFa,
                callRegions = callRegions,
                intervalListBed = intervalListBed,
                bwaReference = bwaReference,
                chromBedsWgs = chromBedsWgs,
                highMem = highMem,
                library = library
        }
    }

    output {
        # Gridss
        Array[IndexedVcf] gridssVcf = Calling.gridssVcf
        # Bicseq2
        Array[File?] bicseq2Png = Calling.bicseq2Png
        Array[File?] bicseq2 = Calling.bicseq2
        # Mutect2
        Array[File] mutect2 = Calling.mutect2
        Array[File] mutect2Unfiltered = Calling.mutect2Unfiltered
        # Manta
        Array[IndexedVcf?] candidateSmallIndels = Calling.candidateSmallIndels
        Array[IndexedVcf?] diploidSV = Calling.diploidSV
        Array[IndexedVcf?] somaticSV = Calling.somaticSV
        Array[IndexedVcf?] candidateSV = Calling.candidateSV
        Array[File?] unfilteredMantaSV = Calling.unfilteredMantaSV
        Array[File?] filteredMantaSV = Calling.filteredMantaSV
        # Strelka2
        Array[IndexedVcf] strelka2Snvs = Calling.strelka2Snvs
        Array[IndexedVcf] strelka2Indels = Calling.strelka2Indels
        Array[File] strelka2Snv = Calling.strelka2Snv
        Array[File] strelka2Indel = Calling.strelka2Indel
        # Lancet
        Array[File] lancet = Calling.lancet
    }
}
