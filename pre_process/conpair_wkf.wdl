version 1.0

import "../wdl_structs.wdl"
import "qc.wdl"

workflow Conpair {
    # command 
    input {
        Bam finalNormalBam
        Bam finalTumorBam
        IndexedReference referenceFa
        String tumor
        String normal
        String pairName
        
        File markerTxtFile
        File markerBedFile
    
        Int threads
        Int memoryGb
        String gatkDockerImage
        String conpairDockerImage
    }
    
    call qc.ConpairPileup as tumorConpairPileup {
        input:
            markerBedFile = markerBedFile,
            referenceFa = referenceFa,
            finalBam = finalTumorBam,
            sampleId = tumor,
            memoryGb = memoryGb,
            threads = threads,
    }
    
    call qc.ConpairPileup as normalConpairPileup {
        input:
            markerBedFile = markerBedFile,
            referenceFa = referenceFa,
            finalBam = finalNormalBam,
            sampleId = normal,
            memoryGb = memoryGb,
            threads = threads,
    }
    
    call qc.VerifyConcordanceAll {
        input:
            pileupsTumorConpair = tumorConpairPileup.pileupsConpair,
            pileupsNormalConpair = normalConpairPileup.pileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
    }
    
    call qc.VerifyConcordanceHomoz {
        input:
            pileupsTumorConpair = tumorConpairPileup.pileupsConpair,
            pileupsNormalConpair = normalConpairPileup.pileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
    }
    
    call qc.Contamination {
        input:
            pileupsTumorConpair = tumorConpairPileup.pileupsConpair,
            pileupsNormalConpair = normalConpairPileup.pileupsConpair,
            markerTxtFile = markerTxtFile,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
    }
    
    output {
        File concordanceAll = VerifyConcordanceAll.concordanceAll
        File concordanceHomoz = VerifyConcordanceHomoz.concordanceHomoz
        File contamination = Contamination.contamination
    }
}

