version 1.0

import "../wdl_structs.wdl"
import "../alignment_analysis/alignment_analysis.wdl" as alignmentAnalysis

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2022) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    James Roche (jroche@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

workflow Reheader {
    input {
        
        Bam finalBam
        String sampleId
    
    }

    # using small disk size because the file is not localized (on servers that support this)
    Int basicDiskSize = 4
    
    call alignmentAnalysis.GetSampleName {
        input:
            finalBam = finalBam.bam,
            finalBai = finalBam.bamIndex,
            diskSize = basicDiskSize
    }
    
    if (GetSampleName.bamSampleId != sampleId ) {
        Int renameDiskSize = (ceil( size(finalBam.bam, "GB") )  * 2 ) + 4
        
        call alignmentAnalysis.UpdateBamSampleName {
            input:
                finalBam = finalBam,
                sampleId = sampleId,
                outputPrefix = sampleId,
                diskSize = renameDiskSize
        }
    }
    
    Bam sampleBam = select_first([UpdateBamSampleName.reheaderBam, finalBam])
    
    output {
        Bam sampleBamMatched = sampleBam
    }
}
