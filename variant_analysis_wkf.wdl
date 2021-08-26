version 1.0

import "variant_analysis/deconstruct_sigs_wkf.wdl" as deconstructSigs
import "wdl_structs.wdl"

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================


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
