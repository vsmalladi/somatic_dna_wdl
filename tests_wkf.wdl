version 1.0

import "test/tests_wkf.wdl" as tests
import "wdl_structs.wdl"

workflow PipelineTests {
    input {
        Array[PairVcfInfo] oldpairVcfInfos
        Array[PairVcfInfo] newpairVcfInfos
        Array[File] oldcnvAnnotatedFinalBed
        Array[File] newcnvAnnotatedFinalBed
        Array[File] oldsvFinalBedPe
        Array[File] newsvFinalBedPe
        IndexedReference referenceFa
    }
    
    call tests.Tests {
        input:
            oldpairVcfInfos = oldpairVcfInfos,
            newpairVcfInfos = newpairVcfInfos,
            oldcnvAnnotatedFinalBed = oldcnvAnnotatedFinalBed,
            newcnvAnnotatedFinalBed = newcnvAnnotatedFinalBed,
            oldsvFinalBedPe = oldsvFinalBedPe,
            newsvFinalBedPe = newsvFinalBedPe,
            referenceFa = referenceFa
    }

    output {

        File detailedOutputTable = Tests.detailedOutputTable
        File summaryOutputTable = Tests.summaryOutputTable
        File BedPeSummary = Tests.BedPeSummary
        File outFileBedpe = Tests.outFileBedpe
        File BedSummary = Tests.BedSummary
        File outFileBed = Tests.outFileBed
    }
}