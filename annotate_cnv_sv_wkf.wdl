version 1.0

import "annotate/annotate_cnv_sv_wkf.wdl" as annotate_cnv_sv
import "wdl_structs.wdl"

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
                svabaSv=pairRawVcfInfo.svabaSv,
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