version 1.0

import "mutect2_wkf.wdl" as mutect2
import "strelka2_wkf.wdl" as strelka2
import "manta_wkf.wdl" as manta
import "svaba_wkf.wdl" as svaba
import "lancet_wkf.wdl" as lancet
import "gridss_wkf.wdl" as gridss
import "bicseq2_wkf.wdl" as bicseq2

import "../wdl_structs.wdl"

workflow Calling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean local = false
        
        pairInfo pairInfo
        #   mutect2
        Array[String]+ listOfChroms
        Array[String]+ listOfChromsFull
        IndexedReference referenceFa
        #   Manta
        IndexedTable callRegions
        #   Svaba
        File dbsnpIndels
        BwaReference bwaReference
        #   Lancet
        Map[String, File] chromBedsWgs

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

        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File svabaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow
        
        Boolean highMem = false
    }
    call mutect2.Mutect2 {
        input:
            local = local,
            mutectJsonLogFilter = mutectJsonLogFilter,
            mutectJsonLog = mutectJsonLog,
            tumor = pairInfo.tumor,
            normal = pairInfo.normal,
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
            tumor = pairInfo.tumor,
            normal = pairInfo.normal,
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
            callRegions = callRegions,
            pairName = pairInfo.pairId,
            normalFinalBam = pairInfo.normalFinalBam,
            tumorFinalBam = pairInfo.tumorFinalBam,
            highMem = highMem
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

    call gridss.Gridss {
        input:
            tumor = pairInfo.tumor,
            normal = pairInfo.normal,
            pairName = pairInfo.pairId,
            bwaReference = bwaReference,
            gridssAdditionalReference = gridssAdditionalReference,
            listOfChroms = listOfChroms,
            normalFinalBam = pairInfo.normalFinalBam,
            tumorFinalBam = pairInfo.tumorFinalBam,
            bsGenome = bsGenome,
            ponTarGz = ponTarGz,
            highMem = highMem
    }

    call bicseq2.BicSeq2 {
        input:
            tumor = pairInfo.tumor,
            normal = pairInfo.normal,
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

    output {
        # Mutect2
        File mutect2 = Mutect2.mutect2
        File mutect2Unfiltered = Mutect2.mutect2_unfiltered
        # Manta
        IndexedVcf candidateSmallIndels = Manta.candidateSmallIndels
        IndexedVcf diploidSV = Manta.diploidSV
        IndexedVcf  somaticSV = Manta.somaticSV
        IndexedVcf  candidateSV = Manta.candidateSV
        File unfilteredMantaSV = Manta.unfilteredMantaSV
        File filteredMantaSV = Manta.filteredMantaSV
        # Strelka2
        IndexedVcf strelka2Snvs = Strelka2.strelka2Snvs
        IndexedVcf strelka2Indels = Strelka2.strelka2Indels
        File strelka2Snv = Strelka2.strelka2Snv
        File strelka2Indel = Strelka2.strelka2Indel
        # Svaba
        File svabaRawGermlineIndel = Svaba.svabaRawGermlineIndel
        File svabaRawGermlineSv = Svaba.svabaRawGermlineSv
        File svabaSv = Svaba.svabaSv
        File svabaIndel = Svaba.svabaIndel
        # Lancet
        File lancet = Lancet.lancet
        # Gridss
        IndexedVcf gridssVcf = Gridss.gridssVcf
        # Bicseq2
        File bicseq2Png = BicSeq2.bicseq2Png
        File bicseq2 = BicSeq2.bicseq2
    }
}
