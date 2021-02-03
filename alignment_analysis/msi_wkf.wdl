version 1.0

import "alignment_analysis.wdl" as alignment_analysis
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
    
    call alignment_analysis.BedtoolsIntersect {
        input:
            mantisBedByIntervalListPath = mantisBedByIntervalListPath,
            mantisBed = mantisBed,
            intervalListBed = intervalListBed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bedtoolsDockerImage
    }
    
    call alignment_analysis.MantisExome {
        input:
            tumorFinalBam = tumorFinalBam,
            normalFinalBam = normalFinalBam,
            mantisBedByIntervalList = BedtoolsIntersect.mantisBedByIntervalList,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = mantisDockerImage,
            referenceFa = referenceFa
    }
    
    call alignment_analysis.MantisRethreshold {
        input:
            normal = normal,
            pairName = pairName,
            mantisWxsStatus = MantisExome.mantisWxsStatus,
            threads = threads,
            dockerImage = mantisDockerImage
    }
}
