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