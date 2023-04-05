version 1.0

import "test/report_mini_wkf.wdl" as reports
import "test/tests.wdl" as tests
import "tasks/utils.wdl" as utils
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
#    James Roche (jroche@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#    Minita Shah
#
# ================== /COPYRIGHT ===============================================

# for wdl version 1.0


workflow PipelineReports {
    input {
        String pipeline = "v7"
        String name = "project"
        File chromLengths
        File cosmicCensus
        Array[String] listOfChroms
        # normal order
        Array[String] normalSampleIds
        # All calls
        Array[FinalVcfPairInfo] finalVcfPairInfos
        # HC
        Array[File] haplotypecallerAnnotatedVcf
        # MSI
        Array[File] mantisStatusFinal
        # baf
        Array[File] alleleCountsTxt
        # Hla
        Array[File] kouramiResult
        # flagStat
        Array[File] flagStat
        # collectWgsMetrics
        Array[File] collectWgsMetrics
        # insertSizeMetrics
        Array[File] insertSizeMetrics
        # qualityByCycleMetrics
        Array[File] qualityByCycleMetrics
        IndexedReference referenceFa
    }
    scatter (finalVcfPairInfo in finalVcfPairInfos) {
        String pairIds = finalVcfPairInfo.pairId
    }
    
    scatter (finalVcfPairInfo in finalVcfPairInfos) {
        
        call reports.SampleReport {
            input:
                chromLengths = chromLengths,
                cosmicCensus = cosmicCensus,
                listOfChroms = listOfChroms,
                referenceFa = referenceFa,
                pairId = finalVcfPairInfo.pairId,
                finalVcfPairInfo = finalVcfPairInfo
                
        }
    }
    
    call tests.SummarizeMantis {
        input:
            mantisStatusFinal = mantisStatusFinal,
            name = name
    }
    
    call tests.SummarizeQualityByCycle {
        input:
            qualityByCycleMetrics = qualityByCycleMetrics,
            name = name
    }
    
    call tests.SummarizeFlagStat {
        input:
            flagStats = flagStat,
            name = name
    }
    
    call tests.SummarizeCollectWgsMetrics {
        input:
            collectWgsMetrics = collectWgsMetrics,
            name = name
    }
    
    call tests.SummarizeInsertSizeMetrics {
        input:
            insertSizeMetrics = insertSizeMetrics,
            name = name
    }
    
    call tests.ConcateTables as summaryVcfTableConcateTables {
        input:
            tables = SampleReport.summaryVcfTable,
            outputTablePath = "~{name}.~{pipeline}.summarySnvIndel.csv"
    }
    
    call tests.SummarizeSvs {
        input:
            highConfidenceSvTables = SampleReport.highConfidenceSvTable,
            allSomaticSvTables = SampleReport.allSomaticSvTable,
            name = name
    }
    
    call tests.ConcateTables as mantisTableConcateTables {
        input:
            tables = mantisStatusFinal,
            outputTablePath = "~{name}.~{pipeline}.summarySnvIndel.csv"
    }

    output {
        File mantisStatusTable = SummarizeMantis.mantisStatusTable
        File qualityByCycleTable = SummarizeQualityByCycle.qualityByCycleTable
        File insertSizeTable = SummarizeInsertSizeMetrics.insertSizeTable
        File coverageTable = SummarizeCollectWgsMetrics.coverageTable
        File summarySvTable = SummarizeSvs.summarySvTable
        File flagStatTable = SummarizeFlagStat.flagStatTable
        File summaryVcfTable = summaryVcfTableConcateTables.outputTable
        
        Array[File] highConfidenceSvTable = SampleReport.highConfidenceSvTable
        Array[File] allSomaticSvTable = SampleReport.allSomaticSvTable
        
        Array[File] detailedVcfTable = SampleReport.detailedVcfTable
        Array[File] summaryVcfTables = SampleReport.summaryVcfTable
        
    }
    
    meta {
        usage : "Run on output from v7 all lists of output must be in the same order (based on the sampleId/pairId scatter used to create the output)"
        notes : "For the moment sampleId order is assumed from filenames (needs updating)"
    }
}