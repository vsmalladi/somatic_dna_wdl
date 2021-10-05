version 1.0


import "calling/mutect2_pon_wkf.wdl" as mutect2Pon
import "calling/manta_pon_wkf.wdl" as mantaPon
import "calling/svaba_pon_wkf.wdl" as svabaPon
import "merge_vcf/merge_vcf_pon_wkf.wdl" as MergeVcfPon

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

workflow CallingPon {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[SampleBamInfo]+ tumorInfos
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
        Boolean highMem = false
    }
    scatter(tumorInfo in tumorInfos) {
        
        call mutect2Pon.Mutect2Pon {
            input:
                mutectJsonLogFilter = mutectJsonLogFilter,
                tumor = tumorInfo.sampleId,
                listOfChroms = listOfChroms,
                referenceFa = referenceFa,
                tumorFinalBam = tumorInfo.finalBam,
                highMem = highMem
        }
        
        call mantaPon.MantaPon {
            input:
                mantaJsonLog = mantaJsonLog,
                tumor = tumorInfo.sampleId,
                callRegions = callRegions,
                referenceFa = referenceFa,
                tumorFinalBam = tumorInfo.finalBam,
                highMem = highMem
        }
        
        call svabaPon.SvabaPon {
            input:
                svabaJsonLog = svabaJsonLog,
                tumor = tumorInfo.sampleId,
                dbsnpIndels = dbsnpIndels,
                referenceFa = referenceFa,
                bwaReference = bwaReference,
                tumorFinalBam = tumorInfo.finalBam,
                highMem = highMem
        }
        call MergeVcfPon.MergeVcfPon as wgsMergeVcfPon {
                    input:
                        filteredMantaSV = MantaPon.filteredMantaSV,
                        mutect2 = Mutect2Pon.mutect2,
                        svabaIndel = SvabaPon.svabaIndel,
                        
                        referenceFa = referenceFa,
                        listOfChroms = listOfChroms
        }
    }
 
    output {
        # Mutect2
        Array[File] mutect2 = Mutect2Pon.mutect2
        # Manta 
        Array[File] filteredMantaSV = MantaPon.filteredMantaSV
        # Strelka2

        # Svaba
        Array[File] svabaIndel = SvabaPon.svabaIndel
        
        Array[File] mergedVcf = wgsMergeVcfPon.mergedVcf
    }
}