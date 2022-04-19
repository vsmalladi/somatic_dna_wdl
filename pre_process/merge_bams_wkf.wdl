version 1.0

import "./merge_bams.wdl" as mergeBams
import "./qc.wdl" as qc
import "../wdl_structs.wdl"

workflow MergeBams {
    # command
    #     merge lane level BAMs
    input {
        Boolean external = false
        #    command merge flowcell
        Array[File] laneFixmateBams
        String sampleId
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        File callRegions
        #IndexedTable callRegions
        IndexedReference referenceFa
        File randomIntervals
        String qcDir = '.'       # To pass to the two QC tasks.
        # resources
        Int memoryGb
        Int threads

    }

    Int laneFixmateBamsSize = ceil(size(laneFixmateBams, "GB"))

    if (!external) { 
        call mergeBams.NovosortMarkDup as novosort {
            input:
                laneBams = laneFixmateBams,
                sampleId = sampleId,
                memoryGb = 20,
                threads = threads,
                # novosort uses a lot of memory and a lot of disk.
                diskSize = (3 * laneFixmateBamsSize)
        }
    }
    
    if (external) { 
        call mergeBams.NovosortMarkDupExternal as novosortExternal {
            input:
                laneBams = laneFixmateBams,
                sampleId = sampleId,
                memoryGb = 20,
                threads = 1,
                # novosort uses a lot of memory and a lot of disk.
                diskSize = (3 * laneFixmateBamsSize)
        }
    }
    
    Bam mergedDedupBam = select_first([novosort.mergedDedupBam, novosortExternal.mergedDedupBam])
    

    # This task runs in parallel with Bqsr38. We are missing coverage check
    # tasks. The idea is that we check coverage and if it's lower than expected
    # coverage, the check task fails thereby stopping the workflow.
    call qc.CollectWgsMetrics {
        input:
            inputBam = mergedDedupBam,
            sampleId = sampleId,
            outputDir = qcDir,
            collectWgsMetricsPath = "~{qcDir}/~{sampleId}.CollectWgsMetrics.dedup.txt",
            referenceFa = referenceFa,
            randomIntervals = randomIntervals,
            diskSize = ceil(size(mergedDedupBam.bam, "GB")) + 10
    }

    call qc.MultipleMetricsPreBqsr {
        input:
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
            outputDir = qcDir,
            sampleId = sampleId,
            diskSize = ceil(size(mergedDedupBam.bam, "GB")) + 10
    }

    call mergeBams.Downsample {
        input:
            mergedDedupBam = mergedDedupBam,
            sampleId = sampleId,
            diskSize = ceil(size(mergedDedupBam.bam, "GB") * 1.5)
    }
    call mergeBams.Bqsr38 {
        input:
            mergedDedupBam = Downsample.downsampleMergedDedupBam,
            MillsAnd1000G = MillsAnd1000G,
            referenceFa = referenceFa,
            Indels = Indels,
            dbsnp = dbsnp,
            callRegions = callRegions,
            sampleId = sampleId,
            diskSize = ceil(size(Downsample.downsampleMergedDedupBam.bam, "GB")) + 20
    }

    call mergeBams.PrintReads {
        input:
            referenceFa = referenceFa,
            mergedDedupBam = mergedDedupBam,
            recalGrp = Bqsr38.recalGrp,
            sampleId = sampleId,
            diskSize = (3 * laneFixmateBamsSize)
    }

    output {
        Bam finalBam = PrintReads.finalBam
        File collectWgsMetricsPreBqsr = CollectWgsMetrics.collectWgsMetrics
        File qualityDistributionPdfPreBqsr = MultipleMetricsPreBqsr.qualityDistributionPdfPreBqsr
        File qualityByCycleMetricsPreBqsr = MultipleMetricsPreBqsr.qualityByCycleMetricsPreBqsr
        File qualityByCyclePdfPreBqsr = MultipleMetricsPreBqsr.qualityByCyclePdfPreBqsr
        File qualityDistributionMetricsPreBqsr = MultipleMetricsPreBqsr.qualityDistributionMetricsPreBqsr
  }
}
