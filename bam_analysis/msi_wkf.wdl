version 1.0

import "bam_analysis.wdl" as bam_analysis
import "../wdl_structs.wdl"

workflow Msi {
    input {
        String normal
        String pairName
        String mantisBedByIntervalListPath
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa
        
        Bam normalFinalBam
        Bam tumorFinalBam
        
        String bedtoolsDockerImage
        String mantisDockerImage
        Int threads
        Int memoryGb
    }
    
    call bam_analysis.BedtoolsIntersect {
        input:
            mantisBedByIntervalListPath = mantisBedByIntervalListPath,
            mantisBed = mantisBed,
            intervalListBed = intervalListBed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bedtoolsDockerImage
    }
    
    call bam_analysis.MantisExome {
        input:
            tumorFinalBam = tumorFinalBam,
            normalFinalBam = normalFinalBam,
            mantisBedByIntervalList = BedtoolsIntersect.mantisBedByIntervalList,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = mantisDockerImage,
            referenceFa = referenceFa
    }
    
    call bam_analysis.MantisRethreshold {
        input:
            normal = normal,
            pairName = pairName,
            mantisWxsStatus = MantisExome.mantisWxsStatus,
            threads = threads,
            dockerImage = mantisDockerImage
    }
}