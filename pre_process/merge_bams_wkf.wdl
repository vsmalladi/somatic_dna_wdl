version 1.0

import "./merge_bams.wdl" as mergeBams
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
        IndexedTable callRegions
        IndexedReference referenceFa
        # resources
        Int mem
        Int threads

    }

    # There has to be a better way to add over a range of array. But I haven't found it.
    call SumFloats {
        input:
            sizes = sample_bam_sizes
    }

    Int diskSize = ceil(SumFloats.total_size) + 50

    call mergeBams.NovosortMarkDup as novosort {
        input:
            laneBams = laneFixmateBams,
            sampleId = sampleId,
            mem = mem,
            threads = threads,
            diskSize = diskSize
    }


    call mergeBams.Bqsr38 {
        input:
            mergedDedupBam = novosort.mergedDedupBam,
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
            diskSize = diskSize
    }

    output {
        Bam mergedDedupBam = novosort.mergedDedupBam
        Bam finalBam = PrintReads.finalBam
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
