version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

workflow MergeChroms {
    input {
        String tumor
        String normal
        String pairName
        String sortedVcfPath = "~{pairName}.sorted.merged.v6.vcf"
        String orderedVcfPath = "~{pairName}.merged.v6.vcf"
        Array[File]+ finalChromVcf
        IndexedReference referenceFa
    }
    
    call merge_vcf.Gatk4MergeSortVcf {
            input:
                tempVcfs = finalChromVcf,
                sortedVcfPath = sortedVcfPath,
                referenceFa = referenceFa,
                memoryGb = 8,
                diskSize = 10
        }
    
    call calling.ReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            orderedVcfPath = orderedVcfPath,
            memoryGb = 2,
            diskSize = 1
    }
    
    output {
        File unannotatedVcf = ReorderVcfColumns.orderedVcf
    }

}