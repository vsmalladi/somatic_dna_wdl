version 1.0

import "../wdl_structs.wdl"

task NovosortMarkDup {
   input {
       # command
       Array[File]+ laneBams
       String sampleId
       String mergedDedupBamPath = "~{sampleId}.merged_dedup.bam"
       # resources
       Int mem = 20
       Int threads = 8
       Int diskSize
   }

    command {
        /bin/novosort \
        -c ~{threads} \
        -m 9216M \
        -i \
        -o ~{mergedDedupBamPath} \
        --forcesort \
        --markDuplicates \
        ${sep=' ' laneBams}
    }

    output {
        Bam mergedDedupBam = object {
            bam : mergedDedupBamPath,
            bamIndex : mergedDedupBamPath + ".bai"
        }
    }

    runtime {
        cpu : threads
        memory : mem + " GB"
        docker :  "gcr.io/nygc-compbio/novosort@sha256:ec3b0a6b1293916df21e69ec668f6ee23f33ce643f8c455c21d5f4234be25193"
        # Per clinical team novosort runs significantly faster with SSD
        disks: "local-disk " + diskSize + " LOCAL"
    }
}

task IndexBam {
    input {
        # command
        File bam
        # resources
        Int diskSize
   }

    command {
        samtools \
        index \
        ~{bam}
    }

    output {
        Bam indexedBam = object {
                bam : bam,
                bamIndex : sub(basename(bam), ".bam$", ".bai")
            }
    }

    runtime {
        docker : "gcr.io/nygc-public/samtools@sha256:963b0b2f24908832efab8ddccb7a7f3ba5dca9803bc099be7cf3a455766610fd"
        disks: "local-disk " + diskSize + " HDD"
    }
}

task Bqsr38 {
    input {
        # command
        Bam mergedDedupBam
        IndexedReference referenceFa
        File callRegions
        #IndexedTable callRegions
        String sampleId
        String recalGrpPath = "~{sampleId}.recal_data.grp"
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        # resources
        Int mem = 12
        Int cpu = 2
        Int diskSize
   }

    Int jvmHeap = mem * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        gatk \
        BaseRecalibrator \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -L ~{callRegions} \
        -R ~{referenceFa.fasta} \
        -I ~{mergedDedupBam.bam} \
        -O ~{recalGrpPath} \
        --known-sites ~{MillsAnd1000G.vcf} \
        --known-sites ~{Indels.vcf} \
        --known-sites ~{dbsnp.vcf}
    }

    output {
        File recalGrp = recalGrpPath
    }

    runtime {
        cpu : cpu
        memory : mem + " GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
        disks: "local-disk " + diskSize + " HDD"
    }
}

task Downsample{
    input {
        String sampleId
        String downsampleMergedDedupBamPath = "~{sampleId}.merged_dedup_10_percent.bam"
        Bam mergedDedupBam
        Int mem = 68  #GB
        Int diskSize
    }
    Int jvmHeap = mem * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk DownsampleSam \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --STRATEGY Chained \
        --RANDOM_SEED 1 \
        --CREATE_INDEX \
        --MAX_RECORDS_IN_RAM 2000000 \
        --VALIDATION_STRINGENCY SILENT \
        -P 0.1 \
        -I ~{mergedDedupBam.bam} \
        -O ~{downsampleMergedDedupBamPath} \
    }
    output {
        Bam downsampleMergedDedupBam = object {
            bam : downsampleMergedDedupBamPath,
            bamIndex : sub(downsampleMergedDedupBamPath, ".bam$", ".bai")
        }
    }
    runtime {
        cpu: 2
        memory: mem + "GB"
        docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
        disks: "local-disk " + diskSize + " HDD"
    }
}

task PrintReads {
    input {
        # command
        Bam mergedDedupBam
        File recalGrp
        IndexedReference referenceFa
        String sampleId
        String finalBamPath = "~{sampleId}.final.bam"
        # resources
        Int mem = 16
        Int cpu = 2
        Int diskSize
    }

    Int jvmHeap = mem * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        gatk ApplyBQSR \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -R ~{referenceFa.fasta} \
        -I ~{mergedDedupBam.bam} \
        -O ~{finalBamPath} \
        --bqsr-recal-file ~{recalGrp}
    }

    output {
        Bam finalBam = object {
                bam : finalBamPath,
                bamIndex : sub(basename(finalBamPath), ".bam$", ".bai")
            }
    }

    runtime {
        cpu : cpu
        memory : mem + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
        disks: "local-disk " + diskSize + " HDD"
        preemptible: 1
     }
}
