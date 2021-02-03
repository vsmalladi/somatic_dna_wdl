version 1.0

import "../wdl_structs.wdl"
import "baf.wdl"


workflow Baf {
    # command 
    input {
        String pairName
        Bam normalFinalBam
        Bam tumorFinalBam
        File finalGermlineVcf
        IndexedReference referenceFa
        Int threads
        Int memoryGb
        String gatkDockerImage
        String bcftoolsDockerImage
        String pysamDockerImage
    }
    
    call baf.MultiMerge {
        input:
            pairName = pairName,
            finalGermlineVcf = finalGermlineVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bcftoolsDockerImage
    }
    
    call baf.FilterForHetSnps {
        input:
            pairName = pairName,
            referenceFa = referenceFa,
            multiallelicVcf = MultiMerge.multiallelicVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call baf.FilterBaf {
        input:
            pairName = pairName,
            hetVcf = FilterForHetSnps.hetVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call baf.AlleleCounts {
        input:
            pairName = pairName,
            referenceFa = referenceFa,
            normalFinalBam = normalFinalBam,
            tumorFinalBam = tumorFinalBam,
            knownHetVcf = FilterBaf.knownHetVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage
    }
    
    call baf.CalcBaf {
        input:
            pairName = pairName,
            alleleCountsTxt = AlleleCounts.alleleCountsTxt,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = pysamDockerImage 
    }
    
    output {
        File alleleCountsTxt = CalcBaf.bafTxt
    }
}