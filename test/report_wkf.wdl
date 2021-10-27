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
        docker: "gcr.io/nygc-internal-tools/workflow_utils@sha256:e9ea11001d93b56c200094d9eb2919e5eef8dbeb66e78f9dc9f43447264e8a48"
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
        # BAF
        File alleleCountsTxt
        # HC
        File haplotypecallerAnnotatedVcf
        # MSI
        File mantisStatusFinal
        # baf
        File alleleCountsTxt
        # Hla
        File kouramiResult
        IndexedReference referenceFa
        
        Int bedDiskSize = ceil( size(finalVcfPairInfo.cnvAnnotatedSupplementalBed, "GB") * 2) + 5
        Int bedPeDiskSize = ceil( size(finalVcfPairInfo.svFinalBedPe, "GB") * 2) + 5
        Int vcfDiskSize = ceil( size(finalVcfPairInfo.mainVcf, "GB") * 2) + 5
        Int bafDiskSize = ceil( size(alleleCountsTxt, "GB") * 2) +  ceil( size(haplotypecallerAnnotatedVcf, "GB")) + 5
    }
    
    call tests.DescribeBed {
        input:
            pairId = finalVcfPairInfo.pairId,
            bed = finalVcfPairInfo.cnvAnnotatedSupplementalBed,
            chromLengths = chromLengths,
            listOfChroms = listOfChroms,
            diskSize = bedDiskSize
    }
    
    call tests.DescribeBedGenes {
        input:
            pairId = finalVcfPairInfo.pairId,
            bed = finalVcfPairInfo.cnvAnnotatedSupplementalBed,
            diskSize = bedDiskSize
    }
    
    call tests.SummarizeFinalVcf {
            input:
                pairId = pairId,
                vcf = finalVcfPairInfo.mainVcf,
                diskSize = vcfDiskSize
        }
        
    call tests.DescribeBedPeGenes {
        input:
            pairId = finalVcfPairInfo.pairId,
            bedpe = finalVcfPairInfo.svHighConfidenceFinalBedPe,
            diskSize = bedPeDiskSize
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
    
    #call tests.DescribeBaf {
    #    input:
    #        pairId = finalVcfPairInfo.pairId,
    #        vcf = haplotypecallerAnnotatedVcf,
    #        baf = alleleCountsTxt,
    #        chromLengths = chromLengths,
    #        listOfChroms = listOfChroms
    #}
    
    #call tests.DraftSampleReport {
    #    input:
    #        listOfChroms = listOfChroms,
    #        pairId = finalVcfPairInfo.pairId,
    
    #        cnvTable = DescribeBed.cnvTable,
    #        cnvGenesTable = DescribeBedGenes.cnvGenesTable,
    
    #        svGenesTable = DescribeBedPeGenes.svGenesTable
    #        svTable = DescribeBedPeGenes.svGenesTable
    
    #        detailedVcfTable = SummarizeFinalVcf.detailedVcfTable,
    #        summaryVcfTable = SummarizeFinalVcf.summaryVcfTable,
    #        diskSize = 10
    #}
    
    
    
    output {
        File cnvTable = DescribeBed.cnvTable
        File cnvGenesTable = DescribeBedGenes.cnvGenesTable
        
        File svGenesTable = DescribeBedPeGenes.svGenesTable
        File highConfidenceSvTable = highConfidenceDescribeBedPe.svTable
        File allSomaticSvTable = allSomaticDescribeBedPe.svTable
        
        File detailedVcfTable = SummarizeFinalVcf.detailedVcfTable
        File summaryVcfTable = SummarizeFinalVcf.summaryVcfTable
        
        #File md = DraftSampleReport.md
        #File header = DraftSampleReport.header

    }
}

