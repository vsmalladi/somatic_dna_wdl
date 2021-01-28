version 1.0

import "../wdl_structs.wdl"
import "qc.wdl"

workflow CalculateContamination {
    # command 
    input {
        Bam finalNormalBam
        Bam finalTumorBam
        IndexedReference referenceFa
        String tumor
        String normal
        String pairName
        
        File gnomadBiallelic
    
        Int threads
        Int memoryGb
        String gatkDockerImage
    }
    
    call qc.Pileup as normalPileup {
        input:
            sampleId = normal,
            finalBam = finalNormalBam,
            gnomadBiallelic = gnomadBiallelic,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage     
    }
    
    call qc.Pileup as tumorPileup {
        input:
            sampleId = tumor,
            finalBam = finalTumorBam,
            gnomadBiallelic = gnomadBiallelic,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage     
    }
    
    call qc.CalculateContaminationPaired {
        input:
            pairName = pairName,
            pileupsNormalTable = normalPileup.pileupsTable,
            pileupsTumorTable = tumorPileup.pileupsTable,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage    
    }
}