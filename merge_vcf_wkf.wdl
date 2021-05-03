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
        File intervalListBed
        
        String library
        File ponWGSFile
        File ponExomeFile

        IndexedVcf germFile
    }
    
    scatter(pairRawVcfInfo in pairRawVcfInfos) {
        if (library == 'WGS') {
            call mergeVcf.MergeVcf as wgsMergeVcf {
                input:
                    pairRawVcfInfo = pairRawVcfInfo,
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
                    pairRawVcfInfo = pairRawVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponExomeFile,
                    germFile = germFile
                    
            }
        }
        
        File mergedVcf = select_first([wgsMergeVcf.mergedVcf, exomeMergeVcf.mergedVcf])
    }
    
    output {
        Array[File] mergedVcf = mergedVcf
    }
}