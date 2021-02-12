version 1.0

import "calling/mutect2_wkf.wdl" as mutect2
import "wdl_structs.wdl"

workflow Calling {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[pairInfo]+ pairInfos
        #    command mutect2
        Array[String]+ listOfChroms
        IndexedReference referenceFa
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
    }

    output {
        Array[File] mutect2 = Mutect2.mutect2
        Array[File] mutect2Unfiltered = Mutect2.mutect2_unfiltered
    }
}