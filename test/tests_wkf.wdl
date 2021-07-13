version 1.0

import "tests.wdl" as tests
import "../wdl_structs.wdl"

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
    }
    
    command {
        python /get_index.py \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }
    
    output {
        Int index = read_int(stdout())
    }
    
    runtime {
        docker: "gcr.io/nygc-internal-tools/workflow_utils:2.0"
    }
}

workflow Tests {
    input {
        String name = "comparison"
        # BAM
        Array[File] newflagStat
        Array[File] oldflagStat
        
        Array[String] pairIds
        
        Array[PairVcfInfo] oldpairVcfInfos
        Array[PairVcfInfo] newpairVcfInfos
        Array[File] oldcnvAnnotatedFinalBed
        Array[File] newcnvAnnotatedFinalBed
        Array[File] oldsvFinalBedPe
        Array[File] newsvFinalBedPe
        # MSI
        Array[File] oldmantisStatusFinal
        Array[File] newmantisStatusFinal
        IndexedReference referenceFa
    }
    
    scatter (flagStat in zip(newflagStat, oldflagStat)) {
        String sampleIdFlagstat = basename(flagStat.left, ".FlagStat.txt") 
        call tests.SummarizeFlagstat {
            input:
                sampleId = sampleIdFlagstat,
                newflagStat = flagStat.left,
                oldflagStat = flagStat.right
        }
    }
    
    scatter (mantisStatusFinal in zip(newmantisStatusFinal, oldmantisStatusFinal)) {
        String sampleId = basename(mantisStatusFinal.left, ".mantis.v1.0.4.WGS-targeted.status.final.tsv") 
        call tests.SummarizeMsi {
            input:
                sampleId = sampleId,
                newmantisStatusFinal = mantisStatusFinal.left,
                oldmantisStatusFinal = mantisStatusFinal.right
        }
    }
    
    scatter (pairId in pairIds) {
        call GetIndex {
            input:
                sampleIds = pairIds,
                sampleId = pairId
        }
        
        call tests.CompareBedPe as finalCompareBedPe {
            input:
                pairId = pairId,
                oldBedpe = oldsvFinalBedPe[GetIndex.index],
                newBedpe = newsvFinalBedPe[GetIndex.index]
        }
        
        call tests.CompareSvGenes as finalCompareSvGenes {
            input:
                pairId = pairId,
                name = name,
                concordanceBedPe = finalCompareBedPe.outFileBedpe
        }
        
        call tests.CompareBed as finalCompareBed {
            input:
                pairId = pairId,
                oldBed = oldcnvAnnotatedFinalBed[GetIndex.index],
                newBed = newcnvAnnotatedFinalBed[GetIndex.index]
        }
        
        call tests.CompareCnvGenes as finalCompareCnvGenes {
            input:
                pairId = pairId,
                name = name,
                concordanceBed = finalCompareBed.outFileBed
        }
        
        call tests.SomPy {
            input:
                referenceFa = referenceFa,
                pairId = pairId,
                oldVcf = oldpairVcfInfos[GetIndex.index].mainVcf,
                newVcf = newpairVcfInfos[GetIndex.index].mainVcf
        }
        
        call tests.SummarizeVcf as finalSummarizeVcf {
            input:
                pairId = pairId,
                name = "final.vcf",
                oldOnlyVcf = SomPy.oldOnlyVcf,
                newOnlyVcf = SomPy.newOnlyVcf,
                concordantVcf = SomPy.concordantVcf
        }
    }
    
    call tests.ConcateTables as finalCnvGenesConcateTables {
        input:
            tables = finalCompareCnvGenes.cnvGenesTable,
            outputTablePath = "~{name}.finalCnvGenes.csv"
    }
    
    call tests.ConcateTables as finalSvGenesConcateTables {
        input:
            tables = finalCompareSvGenes.svGenesTable,
            outputTablePath = "~{name}.finalSvGenes.csv"
    }
    
    call tests.ConcateTables as flagstatConcateTables {
        input:
            tables = SummarizeFlagstat.outputTable,
            outputTablePath = "~{name}.flagstat.csv"
    }
    
    call tests.ConcateTables as msiConcateTables {
        input:
            tables = SummarizeMsi.outputTable,
            outputTablePath = "~{name}.msi.csv"
    }
    
    call tests.ConcateTables as finalBedConcateTables {
        input:
            tables = finalCompareBed.outFileSummary,
            outputTablePath = "~{name}.finalBed.csv"
    }
    
    call tests.ConcateTables as finalBedPeConcateTables {
        input:
            tables = finalCompareBedPe.outFileSummary,
            outputTablePath = "~{name}.finalBedPe.csv"
    }
    
    call tests.ConcateTables as finalVcfConcateTables {
        input:
            tables = finalSummarizeVcf.summaryOutputTable,
            outputTablePath = "~{name}.finalVcf.csv"
    }
    
    output {
        File finalCnvGenesTable = finalCnvGenesConcateTables.outputTable
        File finalSvGenesTable = finalSvGenesConcateTables.outputTable
        File flagstatOutputTable = flagstatConcateTables.outputTable
        File msiOutputTable = msiConcateTables.outputTable
        File finalBedTable = finalBedConcateTables.outputTable
        File finalBedPeTable = finalBedPeConcateTables.outputTable
        File finalVcfTable = finalVcfConcateTables.outputTable
        Array[File] finalBedpeJoined = finalCompareBedPe.outFileBedpe
        Array[File] finalBedJoined = finalCompareBed.outFileBed
        Array[File] detailedOutputVcfTable = finalSummarizeVcf.detailedOutputTable

    }
}

