version 1.0

import "variant_analysis/deconstruct_sigs_wkf.wdl" as deconstructSigs
import "wdl_structs.wdl"

workflow VariantAnalysis {
    input {
        Array[PairVcfInfo]+ pairVcfInfos
        String bsGenome
        File deconstructsigsFasta
    }
    
    scatter(pairVcfInfo in pairVcfInfos) {
        call deconstructSigs.DeconstructSig {
            input:
                mainVcf = pairVcfInfo.mainVcf,
                pairId = pairVcfInfo.pairId,
                bsGenome = bsGenome,
                deconstructsigsFasta = deconstructsigsFasta
        }
        
    }

    output {
        Array[File] diff = DeconstructSig.diff
        Array[File] trinuc = DeconstructSig.trinuc
        Array[File] input_file = DeconstructSig.input_file
        # Array[File] highconfidencePng = DeconstructSig.highconfidencePng
        Array[File] highconfidenceTxt = DeconstructSig.highconfidenceTxt
        Array[File] reconstructed = DeconstructSig.reconstructed
    }
}