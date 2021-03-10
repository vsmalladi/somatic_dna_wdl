version 1.0

import "variant_analysis.wdl" as variant_analysis
import "../wdl_structs.wdl"

workflow DeconstructSig {
    # command 
    #   run DeconstructSig on merged VCF
    input {
        String pairId
        File mainVcf
        
        File deconstructsigsFasta
        String deconstructsigsBs
    }
    
    call variant_analysis.DeconstructsigPrep38 {
        input:
            mainVcf = mainVcf,
            pairId = pairId
        
    }
    
    call variant_analysis.Deconstructsig {
        input:
            highconfidence = DeconstructsigPrep38.highconfidence,
            deconstructsigsBs = deconstructsigsBs,
            deconstructsigsFasta = deconstructsigsFasta,
            pairId = pairId
    }
    
    output {
        File diff = Deconstructsig.diff
        File trinuc = Deconstructsig.trinuc
        File input_file = Deconstructsig.input_file
        # File highconfidencePng = Deconstructsig.highconfidencePng
        File highconfidenceTxt = Deconstructsig.highconfidenceTxt
        File reconstructed = Deconstructsig.reconstructed
    }
}
    
