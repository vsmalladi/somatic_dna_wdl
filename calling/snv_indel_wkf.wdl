version 1.0

import "mutect2_wkf.wdl" as mutect2
import "strelka2_wkf.wdl" as strelka2
import "manta_wkf.wdl" as manta
import "svaba_wkf.wdl" as svaba
import "lancet_wkf.wdl" as lancet

import "../wdl_structs.wdl"

workflow Calling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
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
            mutectJsonLogFilter = mutectJsonLogFilter,
            mutectJsonLog = mutectJsonLog,
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

    call svaba.Svaba {
        input:
            svabaJsonLog = svabaJsonLog,
            tumor = pairInfo.tumorId,
            normal = pairInfo.normalId,
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
            tumor = pairInfo.tumorId,
            normal = pairInfo.normalId,
            listOfChroms = listOfChroms,
            chromBedsWgs = chromBedsWgs,
            referenceFa = referenceFa,
            pairName = pairInfo.pairId,
            normalFinalBam = pairInfo.normalFinalBam,
            tumorFinalBam = pairInfo.tumorFinalBam
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
    }
}
