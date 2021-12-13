version 1.0

import "./merge_bams.wdl" as mergeBams
import "./qc.wdl" as qc
import "../wdl_structs.wdl"

workflow MergeBams {
    # command
    #     merge lane level BAMs
    input {
        #    command merge flowcell
        Array[File] laneFixmateBams
        Array[Int] sample_bam_sizes
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
        Int novosortMem
        Int threads

    }

    # There has to be a better way to add over a range of array. But I haven't found it.
    call SumFloats {
        input:
            sizes = sample_bam_sizes
    }

    Int diskSize = ceil(SumFloats.total_size * 2)

    call mergeBams.NovosortMarkDup as novosort {
        input:
            laneBams = laneFixmateBams,
            sampleId = sampleId,
            logDir = qcDir,      # Send novosort log to qc dir.
            mem = novosortMem,
            threads = threads,
            # novosort uses a lot of memory and a lot of disk.
            diskSize = ceil((SumFloats.total_size * 3))
    }

    # This task runs in parallel with Bqsr38. We are missing coverage check
    # tasks. The idea is that we check coverage and if it's lower than expected
    # coverage, the check task fails thereby stopping the workflow.
    call qc.CollectWgsMetrics {
        input:
            inputBam = novosort.mergedDedupBam,
            sampleId = sampleId,
            outputDir = qcDir,
            collectWgsMetricsPath = "~{qcDir}/~{sampleId}.CollectWgsMetrics.dedup.txt",
            referenceFa = referenceFa,
            randomIntervals = randomIntervals,
            diskSize = diskSize
    }

    call qc.MultipleMetricsPreBqsr {
        input:
            referenceFa = referenceFa,
            mergedDedupBam = novosort.mergedDedupBam,
            outputDir = qcDir,
            sampleId = sampleId,
            diskSize = diskSize
    }

    call mergeBams.Downsample {
        input:
            mergedDedupBam = novosort.mergedDedupBam,
            sampleId = sampleId,
            diskSize = diskSize
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
            diskSize = diskSize
    }

    call mergeBams.PrintReads {
        input:
            referenceFa = referenceFa,
            mergedDedupBam = novosort.mergedDedupBam,
            recalGrp = Bqsr38.recalGrp,
            sampleId = sampleId,
            diskSize = ceil((SumFloats.total_size * 3))
    }

    output {
        Bam finalBam = PrintReads.finalBam
        File collectWgsMetricsPreBqsr = CollectWgsMetrics.collectWgsMetrics
        File qualityDistributionPdfPreBqsr = MultipleMetricsPreBqsr.qualityDistributionPdfPreBqsr
        File qualityByCycleMetricsPreBqsr = MultipleMetricsPreBqsr.qualityByCycleMetricsPreBqsr
        File qualityByCyclePdfPreBqsr = MultipleMetricsPreBqsr.qualityByCyclePdfPreBqsr
        File qualityDistributionMetricsPreBqsr = MultipleMetricsPreBqsr.qualityDistributionMetricsPreBqsr
        File dedupLog = novosort.dedupLog
  }
}

# This task should live in some shared utils.
task SumFloats {
    input {
        Array[Float] sizes
    }

    command {
        python -c "print ~{sep="+" sizes}"
    }

    output {
        Float total_size = read_float(stdout())
    }

    runtime {
        docker: "python:2.7"
    }
}
