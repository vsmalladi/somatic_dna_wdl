version 1.0

import "test/tests_wkf.wdl" as tests
import "wdl_structs.wdl"


workflow PipelineTests {
    input {
        # BAMs
        Array[File] newflagStat
        Array[File] oldflagStat
        # SNV/INDEL calls
        Array[PairVcfInfo] oldpairVcfInfos
        Array[PairVcfInfo] newpairVcfInfos
        # CNV/SV calls
        Array[File] oldcnvAnnotatedFinalBed
        Array[File] newcnvAnnotatedFinalBed
        Array[File] oldsvFinalBedPe
        Array[File] newsvFinalBedPe
        # MSI
        Array[File] oldmantisStatusFinal
        Array[File] newmantisStatusFinal
        IndexedReference referenceFa
    }
    
    scatter (pairVcfInfo in oldpairVcfInfos) {
        String pairIds = pairVcfInfo.pairId
    }
    
    call tests.Tests {
        input:
            newflagStat = newflagStat,
            oldflagStat = oldflagStat,
            
            oldmantisStatusFinal = oldmantisStatusFinal,
            newmantisStatusFinal = newmantisStatusFinal,
            
            pairIds = pairIds,
            
            oldpairVcfInfos = oldpairVcfInfos,
            newpairVcfInfos = newpairVcfInfos,
            oldcnvAnnotatedFinalBed = oldcnvAnnotatedFinalBed,
            newcnvAnnotatedFinalBed = newcnvAnnotatedFinalBed,
            oldsvFinalBedPe = oldsvFinalBedPe,
            newsvFinalBedPe = newsvFinalBedPe,
            referenceFa = referenceFa
    }

    output {
        File finalCnvGenesTable = Tests.finalCnvGenesTable
        File finalSvGenesTable = Tests.finalSvGenesTable
        File flagstatOutputTable = Tests.flagstatOutputTable
        File msiOutputTable = Tests.msiOutputTable
        File finalBedTable = Tests.finalBedTable
        File finalBedPeTable = Tests.finalBedPeTable
        File finalVcfTable = Tests.finalVcfTable
        Array[File] finalBedpeJoined = Tests.finalBedpeJoined
        Array[File] finalBedJoined = Tests.finalBedJoined
        Array[File] detailedOutputVcfTable = Tests.detailedOutputVcfTable
    }
    
    meta {
        usage : "Run on output from v7 all lists of output must be in the same order (based on the sampleId/pairId scatter used to create the output)"
        notes : "Append new or old to the begining of the WDL name found in the cromwell log or the custom outputInfo.json. For the moment sampleId order is assumed from filenames (needs updating)"
    }
}