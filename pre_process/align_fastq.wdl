version 1.0

import "../wdl_structs.wdl"

task Skewer {
    input {
        # command
        Fastqs fastqs
        File adaptersFa
        # replace .fastq.gz or .fq.gz with "-trimmed-pair2.fastq.gz"
        String fastqPrefix = sub(basename(fastqs.fastqR1), ".\\.f.*q\\.gz$", "")
        String fastqOutR1Path = sub(basename(fastqs.fastqR1), ".\\.f.*q\\.gz$", "-trimmed-pair1.fastq.gz")
        String fastqOutR2Path = sub(basename(fastqs.fastqR1), ".\\.f.*q\\.gz$", "-trimmed-pair2.fastq.gz")
        # resources
        Int threads = 16
        Int mem = 24
        Int diskSize
    }

    command {
        skewer \
        --compress \
        -x adaptersFa \
        -f "sanger" \
        -t ~{threads} \
        -m "pe" \
        -o ~{fastqPrefix} \
        ~{fastqs.fastqR1} \
        ~{fastqs.fastqR2}
    }

    output {
        File fastqOutR1 = "~{fastqOutR1Path}"
        File fastqOutR2 = "~{fastqOutR2Path}"
    }

    runtime {
        cpu : threads
        memory : mem + " GB"
        docker : "gcr.io/nygc-public/skewer:v0.2.2"
        disks: "local-disk " + diskSize + " HDD"
    }
}

task AlignBwaMem {
    input {
        # command
        Fastqs fastqs
        File fastqR1
        File fastqR2
        BwaReference bwaReference
        String laneBamPath = "~{fastqs.readgroupId}.readgroup.bam"
        # resources
        Int mem
        Int threads
        Int bwaThreads
        Int totalThreads = 80
        Int diskSize

        # Values used in RG tags. These are overridden for external fastqs or if we start using
        # other sequencing platforms.
        String platform = "illumina"
        String machineType = "NovaSeq"
        String center = "NYGenome"

    }

    command {
        set -e -o pipefail
        bwa mem \
        -Y \
        -K 100000000 \
        -t ~{bwaThreads} \
        -R '@RG\tID:~{fastqs.readgroupId}\tPL:~{platform}\tPM:~{machineType}\tLB:~{fastqs.sampleId}\tDS:hg38\tSM:~{fastqs.sampleId}\tCN:~{center}\tPU:${fastqs.rgpu}' \
        ~{bwaReference.fasta} \
        ~{fastqR1} \
        ~{fastqR2} \
        | samtools view \
        -@ ~{threads} \
        -Shb \
        -o ~{laneBamPath} \
        -
    }

    output {
        File laneBam = laneBamPath
    }

    runtime {
        cpu : totalThreads
        memory : mem + " GB"
        disks: "local-disk " + diskSize + " LOCAL"
        docker : "gcr.io/nygc-public/bwa-kit:0.7.15"
    }
}

task ShortAlignMark {
    input {
        # command
        File laneBam
        String bamBase
        String laneBamMarkPath = "~{bamBase}.readgroup_mark.bam"
        # resources
        Int mem = 16
        Int diskSize
    }

    command {
        set -e -o pipefail
        filter_bam \
        -I ~{laneBam} \
        -A1 30 \
        -A2 30 \
        -o ~{laneBamMarkPath} \
        | samtools view \
        -b \
        -o ~{laneBamMarkPath} \
        -
    }

    output {
        File laneBamMark = laneBamMarkPath
    }

    runtime {
        memory : mem + " GB"
        docker : "gcr.io/nygc-public/nygc-short-alignment-marking:v2.1"
        disks: "local-disk " + diskSize + " HDD"
        preemptible: 1
    }
}

task Fixmate {
    input {
        #command
        File laneBamMark
        String bamBase
        String laneFixmateBamPath = "~{bamBase}.readgroup_fixmate.bam"
        # resources
        Int mem = 8
        Int diskSize
    }

    Int jvmHeap = mem * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        gatk \
        FixMateInformation \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=1" \
        --MAX_RECORDS_IN_RAM 2000000 \
        --VALIDATION_STRINGENCY SILENT \
        --ADD_MATE_CIGAR true \
        --ASSUME_SORTED true \
        -I ~{laneBamMark} \
        -O ~{laneFixmateBamPath}
    }

    output {
        File laneFixmateBam = laneFixmateBamPath
    }

    runtime {
        memory : mem + " GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
        disks: "local-disk " + diskSize + " HDD"
    }
}
