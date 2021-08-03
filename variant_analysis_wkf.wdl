version 1.0

import "variant_analysis/deconstruct_sigs_wkf.wdl" as deconstructSigs
import "wdl_structs.wdl"

workflow VariantAnalysis {
    input {
        Array[MergedPairVcfInfo]+ mergedPairVcfInfos
        String vepGenomeBuild
        File cosmicSigs
    }
    
    scatter(mergedPairVcfInfo in mergedPairVcfInfos) {
        call deconstructSigs.DeconstructSig {
            input:
                mainVcf = mergedPairVcfInfo.unannotatedVcf,
                pairId = mergedPairVcfInfo.pairId,
                vepGenomeBuild = vepGenomeBuild,
                cosmicSigs = cosmicSigs
        }
        
    }

    output {
         Array[File] sigs = DeconstructSig.sigs
         Array[File] counts = DeconstructSig.counts
         Array[File] sig_input = DeconstructSig.sigInput
         Array[File] reconstructed = DeconstructSig.reconstructed
         Array[File] diff = DeconstructSig.diff
    }
}
