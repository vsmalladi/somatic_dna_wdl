version 1.0

import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
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

workflow MergeVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[PairRawVcfInfo]+ pairRawVcfInfos
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        # merge callers
        File intervalListBed
        
        String library
        File ponWGSFile
        File ponExomeFile

        IndexedVcf germFile
    }
    
    scatter(pairRawVcfInfo in pairRawVcfInfos) {
        if (library == 'WGS') {
            call mergeVcf.MergeVcf as wgsMergeVcf {
                input:
                    pairRawVcfInfo = pairRawVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponWGSFile,
                    germFile = germFile
                    
            }
        }
        
        if (library == 'Exome') {
            call mergeVcf.MergeVcf as exomeMergeVcf {
                input:
                    pairRawVcfInfo = pairRawVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponExomeFile,
                    germFile = germFile
                    
            }
        }
        
        File mergedVcfFile = select_first([wgsMergeVcf.mergedVcf, exomeMergeVcf.mergedVcf])
    }
    
    output {
        Array[File] mergedVcf = mergedVcfFile
    }
}