version 1.0

import "merge_vcf.wdl" as merge_vcf
import "../calling/calling.wdl" as calling
import "../wdl_structs.wdl"

workflow MergeChroms {
    input {
        String tumor
        String normal
        String pairName
        String sortedVcfPath = "~{pairName}.merged.v6.vcf"
        String orderedVcfPath = "~{pairName}.merged.v6.vcf"
        Array[File]+ finalChromVcf
        IndexedReference referenceFa
        Int threads
        Int memory_gb
        String gatkDockerImage
        String pysamDockerImage
    }
    
    call merge_vcf.Gatk4MergeSortVcf {
            input:
                tempVcfs = finalChromVcf,
                sortedVcfPath = sortedVcfPath,
                referenceFa = referenceFa,
                memory_gb = memory_gb,
                threads = threads,
                dockerImage = gatkDockerImage 
        }
    
    call calling.ReorderVcfColumns {
        input:
            tumor = tumor,
            normal = normal,
            rawVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            orderedVcfPath = orderedVcfPath,
            memory_gb = memory_gb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    output {
        File unannotatedVcf = ReorderVcfColumns.orderedVcf
    }

}