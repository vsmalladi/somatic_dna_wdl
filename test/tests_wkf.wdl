version 1.0

import "tests.wdl" as tests
import "../wdl_structs.wdl"

workflow Tests {
    input {
        Array[PairVcfInfo] oldpairVcfInfos
        Array[PairVcfInfo] newpairVcfInfos
        Array[File] oldcnvAnnotatedFinalBed
        Array[File] newcnvAnnotatedFinalBed
        Array[File] oldsvFinalBedPe
        Array[File] newsvFinalBedPe
        IndexedReference referenceFa
    }
    
    call tests.CompareBedPe {
        input:
            pairId = oldpairVcfInfos[0].pairId,
            oldBedpe = oldsvFinalBedPe[0],
            newBedpe = newsvFinalBedPe[0]
    }
    
    call tests.CompareBed {
        input:
            pairId = oldpairVcfInfos[0].pairId,
            oldBed = oldcnvAnnotatedFinalBed[0],
            newBed = newcnvAnnotatedFinalBed[0]
    }
    
    call tests.SomPy {
        input:
            referenceFa = referenceFa,
            pairId = oldpairVcfInfos[0].pairId,
            oldVcf = oldpairVcfInfos[0].mainVcf,
            newVcf = newpairVcfInfos[0].mainVcf
    }
    
    call tests.SummarizeVcf {
        input:
            pairId = oldpairVcfInfos[0].pairId,
            name = "final.vcf",
            oldOnlyVcf = SomPy.oldOnlyVcf,
            newOnlyVcf = SomPy.newOnlyVcf,
            concordantVcf = SomPy.concordantVcf
    }
    
    output {
        File BedPeSummary = CompareBedPe.outFileSummary
        File outFileBedpe = CompareBedPe.outFileBedpe
        File BedSummary = CompareBed.outFileSummary
        File outFileBed = CompareBed.outFileBed
        File detailedOutputTable = SummarizeVcf.detailedOutputTable
        File summaryOutputTable = SummarizeVcf.summaryOutputTable
    }
}

