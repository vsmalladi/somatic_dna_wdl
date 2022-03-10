version 1.0

import "tests.wdl" as tests
import "../wdl_structs.wdl"

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
    }
    
    command {
        python /get_index.py \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }
    
    output {
        Int index = read_int(stdout())
    }
    
    runtime {
        docker: "gcr.io/nygc-public/workflow_utils@sha256:b7269061a4620c6565566cbeaf61b1a58d49d26c382fa12f05f41b0e5f2e4807"
    }
}

workflow SampleReport {
    input {
        String pipeline = "v7"
        String pairId
        File chromLengths
        File cosmicCensus
        Array[String] listOfChroms
        # All calls
        FinalVcfPairInfo finalVcfPairInfo
        IndexedReference referenceFa
        
        Int bedPeDiskSize = ceil( size(finalVcfPairInfo.svFinalBedPe, "GB") * 2) + 5
        Int vcfDiskSize = ceil( size(finalVcfPairInfo.mainVcf, "GB") * 2) + 5
    }
    
    call tests.SummarizeFinalVcf {
            input:
                pairId = pairId,
                vcf = finalVcfPairInfo.mainVcf,
                diskSize = vcfDiskSize
        }
    
    call tests.DescribeBedPe as highConfidenceDescribeBedPe {
        input:
            pairId = finalVcfPairInfo.pairId,
            name = ".highConfidence",
            bedpe = finalVcfPairInfo.svHighConfidenceFinalBedPe,
            chromLengths = chromLengths,
            listOfChroms = listOfChroms,
            diskSize = bedPeDiskSize
    }
    
    call tests.DescribeBedPe as allSomaticDescribeBedPe {
        input:
            pairId = finalVcfPairInfo.pairId,
            name = ".AllSomatic",
            bedpe = finalVcfPairInfo.svFinalBedPe,
            chromLengths = chromLengths,
            listOfChroms = listOfChroms,
            diskSize = bedPeDiskSize
    }
    
    
    
    
    output {
        
        File highConfidenceSvTable = highConfidenceDescribeBedPe.svTable
        File allSomaticSvTable = allSomaticDescribeBedPe.svTable
        
        File detailedVcfTable = SummarizeFinalVcf.detailedVcfTable
        File summaryVcfTable = SummarizeFinalVcf.summaryVcfTable

    }
}

