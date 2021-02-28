version 1.0

import "variant_analysis/deconstruct_sigs_wkf.wdl" as deconstructSigs
import "wdl_structs.wdl"

workflow VariantAnalysis {
    input {
        Array[PairVcfInfo]+ PairVcfInfos
    }
    
    scatter(PairVcfInfo in PairVcfInfos) {
        call deconstructSigs.DeconstructSig {
            input:
                mainVcf = PairVcfInfo.mainVcf,
                pairId = PairVcfInfo.pairId
        }
        
    }

    output {
        Array[File] diff = DeconstructSig.diff
        Array[File] trinuc = DeconstructSig.trinuc
        Array[File] input_file = DeconstructSig.input_file
        Array[File] highconfidencePng = DeconstructSig.highconfidencePng
        Array[File] highconfidenceTxt = DeconstructSig.highconfidenceTxt
        Array[File] reconstructed = DeconstructSig.reconstructed
    }
}