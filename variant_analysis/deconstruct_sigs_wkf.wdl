version 1.0

import "variant_analysis.wdl" as variant_analysis
import "../wdl_structs.wdl"

workflow DeconstructSig {
    # command 
    #   run DeconstructSig on merged VCF
    input {
        String pairName
        File mainVcf
        File header
        
        IndexedReference customFasta
        File deconstructsigsBs
        
        Int threads
        Int memoryGb
        String deconstructSigDockerImage
        String bcftoolsDockerImage
    }
    
    call variant_analysis.DeconstructsigPrep38 {
        input:
            mainVcf = mainVcf,
            header = header,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bcftoolsDockerImage
        
    }
    
    call variant_analysis.Deconstructsig {
        input:
            highconfidence = DeconstructsigPrep38.highconfidence,
            deconstructsigsBs = deconstructsigsBs,
            customFasta = customFasta,
            pairName = pairName,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = deconstructSigDockerImage
    }
    
    output {
        File diff = Deconstructsig.diff
        File trinuc = Deconstructsig.trinuc
        File input_file = Deconstructsig.input_file
        File highconfidencePng = Deconstructsig.highconfidencePng
        File highconfidenceTxt = Deconstructsig.highconfidenceTxt
        File reconstructed = Deconstructsig.reconstructed
    }
}
    
