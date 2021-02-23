version 1.0

import "../wdl_structs.wdl"

task AlignBwaMem {
    input {
        # command
        Fastqs fastqs
        BwaReference bwaReference
        String laneBamPath = "~{fastqs.readgroupId}.readgroup.bam"
        # resources
        Int mem
        Int threads
        Int diskSize
    }

    command {
        set -e -o pipefail
        bwa mem \
        -Y \
        -K 100000000 \
        -t ~{threads} \
        -R '@RG\tID:~{fastqs.readgroupId}\tPL:illumina\tPM:NovaSeq\tLB:~{fastqs.sampleId}\tDS:hg38\tSM:~{fastqs.sampleId}\tCN:NYGenome\tPU:${fastqs.flowcell}.${fastqs.lane}.${fastqs.barcode}' \
        ~{bwaReference.fasta} \
        ~{fastqs.fastqR1} \
        ~{fastqs.fastqR2} \
        | samtools view \
        -Shb \
        -o ~{laneBamPath} \
        -

    }

    output {
        File laneBam = laneBamPath
    }

    runtime {
        cpu : threads
        memory : mem + " GB"
        disks: "local-disk " + diskSize + " HDD"
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

    command {
        gatk \
        FixMateInformation \
        --java-options -XX:ParallelGCThreads=1 \
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
