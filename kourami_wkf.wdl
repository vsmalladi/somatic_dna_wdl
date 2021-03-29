version 1.0

import "alignment_analysis/kourami_wfk.wdl" as kourami
import "wdl_structs.wdl"

workflow AlignmentAnalysis {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[SampleBamInfo]+ sampleBamInfos
        
        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index

    }
    
    scatter(sampleBamInfo in sampleBamInfos) {
        call kourami.Kourami {
            input:
                sampleId=sampleBamInfo.sampleId,
                kouramiReference=kouramiReference,
                finalBam=sampleBamInfo.finalBam,
                kouramiFastaGem1Index=kouramiFastaGem1Index
        }
    }
    
    output {
        Array[File] kouramiResult = Kourami.result
    }
}