version 1.0

import "alignment_analysis/kourami_wfk.wdl" as kourami
import "alignment_analysis/msi_wkf.wdl" as msi
import "wdl_structs.wdl"

workflow AlignmentAnalysis {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[pairInfo]+ pairInfos
        
        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index
        
        # mantis
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa
    }
    
    scatter(pairInfo in pairInfos) {
        call kourami.Kourami as kouramiNormal {
            input:
                sampleId=pairInfo.normal,
                kouramiReference=kouramiReference,
                finalBam=pairInfo.normalFinalBam,
                kouramiFastaGem1Index=kouramiFastaGem1Index
        }
        
        call kourami.Kourami as kouramiTumor {
            input:
                sampleId=pairInfo.tumor,
                kouramiReference=kouramiReference,
                finalBam=pairInfo.tumorFinalBam,
                kouramiFastaGem1Index=kouramiFastaGem1Index
        }
        
        call msi.Msi {
            input:
                normal=pairInfo.normal,
                pairName=pairInfo.pairId,
                mantisBed=mantisBed,
                intervalListBed=intervalListBed,
                referenceFa=referenceFa,
                tumorFinalBam=pairInfo.tumorFinalBam,
                normalFinalBam=pairInfo.normalFinalBam
        }
    }
    
    output {
        Array[File] normalKouramiResult = kouramiNormal.result
        Array[File] tumorKouramiResult = kouramiTumor.result
        Array[File] mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        Array[File] mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        Array[File] mantisExomeTxt = Msi.mantisExomeTxt
        Array[File] mantisStatusFinal = Msi.mantisStatusFinal
    }
}