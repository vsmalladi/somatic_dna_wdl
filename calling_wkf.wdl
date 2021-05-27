version 1.0

import "calling/gridss_wkf.wdl" as gridss
import "calling/bicseq2_wkf.wdl" as bicseq2
import "calling/mutect2_wkf.wdl" as mutect2
import "calling/strelka2_wkf.wdl" as strelka2
import "calling/manta_wkf.wdl" as manta
import "calling/svaba_wkf.wdl" as svaba
import "calling/lancet_wkf.wdl" as lancet
import "calling/gridss_wkf.wdl" as gridss
import "calling/bicseq2_wkf.wdl" as bicseq2


import "wdl_structs.wdl"

workflow Calling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[pairInfo]+ pairInfos
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
        #   Svaba
        File dbsnpIndels
        File svabaJsonLog
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
        Map[String, Map[String, File]] bicseq2ConfigMaps
        Map[String, File] chromFastas
        Int tumorMedianInsertSize = 400
        Int normalMedianInsertSize = 400
        Int lambda = 4
        # Gridss
        String bsGenome
        File ponTarGz
        Array[File] gridssAdditionalReference
    }
    scatter(pairInfo in pairInfos) {
        call gridss.Gridss {
            input:
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                pairName = pairInfo.pairId,
                bwaReference = bwaReference,
                gridssAdditionalReference = gridssAdditionalReference,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam,
                bsGenome = bsGenome,
                ponTarGz = ponTarGz
        }    
    
        call bicseq2.BicSeq2 {
            input:
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                readLength = readLength,
                coordReadLength = coordReadLength,
                uniqCoords = uniqCoords,
                tumorConfigFile = bicseq2ConfigMaps[pairInfo.pairId]["tumorConfigFile"],
                normalConfigFile = bicseq2ConfigMaps[pairInfo.pairId]["normalConfigFile"],
                segConfigFile = bicseq2ConfigMaps[pairInfo.pairId]["segConfigFile"],
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
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                listOfChroms = listOfChroms,
                pairName = pairInfo.pairId,
                referenceFa = referenceFa,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
        
        call manta.Manta {
            input:
                mantaJsonLog = mantaJsonLog,
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                callRegions = callRegions,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
        
        call strelka2.Strelka2 {
            input:
                strelkaJsonLog = strelkaJsonLog,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                callRegions = callRegions,
                candidateSmallIndels = Manta.candidateSmallIndels,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
        
        call svaba.Svaba {
            input:
                svabaJsonLog = svabaJsonLog,
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                dbsnpIndels = dbsnpIndels,
                bwaReference = bwaReference,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
        
        call lancet.Lancet {
            input:
                lancetJsonLog = lancetJsonLog,
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
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
        Array[File] gridssVcf = Gridss.gridssVcf
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
        # Svaba
        Array[File] svabaRawGermlineIndel = Svaba.svabaRawGermlineIndel
        Array[File] svabaRawGermlineSv = Svaba.svabaRawGermlineSv
        Array[File] svabaSv = Svaba.svabaSv
        Array[File] svabaIndel = Svaba.svabaIndel
        # Lancet
        Array[File] lancet = Lancet.lancet
    }
}