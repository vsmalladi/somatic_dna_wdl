version 1.0

import "annotate/annotate_cnv_sv_wkf.wdl" as annotate_cnv_sv
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

workflow AnnotateCnvSv {
    input {
        Array[PairRawVcfInfo]+ pairRawVcfInfos
        Array[String] listOfChroms
        
        # cnv
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed
        
        # sv
        String vepGenomeBuild
        
        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe
    }
    
    scatter(pairRawVcfInfo in pairRawVcfInfos) {
        call annotate_cnv_sv.AnnotateCnvSv {
            input:
                tumor=pairRawVcfInfo.tumor,
                normal=pairRawVcfInfo.normal,
                pairName=pairRawVcfInfo.pairId,
                listOfChroms=listOfChroms,
                bicseq2=pairRawVcfInfo.bicseq2,
                cytoBand=cytoBand,
                dgv=dgv,
                thousandG=thousandG,
                cosmicUniqueBed=cosmicUniqueBed,
                cancerCensusBed=cancerCensusBed, 
                ensemblUniqueBed=ensemblUniqueBed,
                
                filteredMantaSV=pairRawVcfInfo.filteredMantaSV,
                gridssVcf=pairRawVcfInfo.gridssVcf,
                vepGenomeBuild=vepGenomeBuild,
                gap=gap,
                dgvBedpe=dgvBedpe,
                thousandGVcf=thousandGVcf,
                svPon=svPon,
                cosmicBedPe=cosmicBedPe
                
        }
    }
    
    output {
        Array[File] cnvAnnotatedFinalBed  = AnnotateCnvSv.cnvAnnotatedFinalBed
        Array[File] cnvAnnotatedSupplementalBed  = AnnotateCnvSv.cnvAnnotatedSupplementalBed
        
        Array[File] svFinalBedPe = AnnotateCnvSv.svFinalBedPe
        Array[File] svHighConfidenceFinalBedPe = AnnotateCnvSv.svHighConfidenceFinalBedPe
    
        Array[File] svSupplementalBedPe = AnnotateCnvSv.svSupplementalBedPe
        Array[File] svHighConfidenceSupplementalBedPe = AnnotateCnvSv.svHighConfidenceSupplementalBedPe
    }
}