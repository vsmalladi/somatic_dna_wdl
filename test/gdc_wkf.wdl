version 1.0

import "tests.wdl" as tests
import "../wdl_structs.wdl"
import "../merge_vcf/merge_vcf.wdl" as mergeVcf


workflow GdcComparePair {
    input {
        String name = "gdc_comparison"
        
        String pairId
        File nygcVcf
        File gdcSnvVcf
        File gdcIndelVcf
    
        IndexedReference referenceFa
    }
 
        
    call mergeVcf.Gatk4MergeSortVcf {
        input:
            tempVcfs = [gdcSnvVcf, gdcIndelVcf],
            sortedVcfPath = "~{pairId}.gdc.vcf",
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }
    
    call tests.FilterHighConfidence {
        input:
            vcf = nygcVcf,
            filteredVcfPath = "~{pairId}.nygc.filtered.vcf",
    }
    
    call tests.SomPy {
        input:
            referenceFa = referenceFa,
            pairId = pairId,
            oldVcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            newVcf = FilterHighConfidence.filteredVcf
    }
    
    call tests.SummarizeVcf as finalSummarizeVcf {
        input:
            pairId = pairId,
            name = "final.vcf",
            oldOnlyVcf = SomPy.oldOnlyVcf,
            newOnlyVcf = SomPy.newOnlyVcf,
            concordantVcf = SomPy.concordantVcf
    }

    
    
    output {
        File summaryOutputTable = finalSummarizeVcf.summaryOutputTable
        File detailedOutputTable = finalSummarizeVcf.detailedOutputTable

    }
}

