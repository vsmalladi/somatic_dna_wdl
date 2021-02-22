version 1.0

import "calling/mutect2_wkf.wdl" as mutect2
import "calling/strelka2_wkf.wdl" as strelka2
import "calling/manta_wkf.wdl" as manta
import "calling/svaba_wkf.wdl" as svaba
import "calling/lancet_wkf.wdl" as lancet
import "wdl_structs.wdl"

workflow Calling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[pairInfo]+ pairInfos
        #   mutect2
        Array[String]+ listOfChroms
        IndexedReference referenceFa
        #   Manta
        IndexedTable callRegions
        #   Svaba
        IndexedVcf dbsnp
        #   Lancet
        Map[String, File] chromBeds
    }
    scatter(pairInfo in pairInfos) {
        call mutect2.Mutect2 {
            input:
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
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                dbsnp = dbsnp,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
        
        call lancet.Lancet {
            input:
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                listOfChroms = listOfChroms,
                chromBeds = chromBeds,
                referenceFa = referenceFa,
                pairName = pairInfo.pairId,
                normalFinalBam = pairInfo.normalFinalBam,
                tumorFinalBam = pairInfo.tumorFinalBam
        }
    }

    output {
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
        Array[Array[File]] svabaInternalInput = Svaba.svabaInternalInput
        Array[File] svabaSv = Svaba.svabaSv
        Array[File] svabaIndel = Svaba.svabaIndel
        # Lancet
        Array[File] lancet = Lancet.lancet
    }
}