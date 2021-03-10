version 1.0

import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "wdl_structs.wdl"

workflow MergeVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[PairRawVcfInfo]+ pairRawVcfInfos
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        # merge callers
        File knownGeneBed
        File ponFile
        IndexedVcf germFile
    }
    
    scatter(pairRawVcfInfo in pairRawVcfInfos) {
        call mergeVcf.MergeVcf {
            input:
                pairRawVcfInfo = pairRawVcfInfo,
                referenceFa = referenceFa,
                listOfChroms = listOfChroms,
                knownGeneBed = knownGeneBed,
                ponFile = ponFile,
                germFile = germFile
                
        }
    }
    
    output {
        Array[File] mergedVcfs = MergeVcf.mergedVcf
    }
}