version 1.0

import "../wdl_structs.wdl"

task NovosortMarkDup {
   input {
       # command
       Array[File]+ laneBams
       String sampleId
       String mergedDedupBamPath = "~{sampleId}.merged_dedup.bam"
       # resources
       Int mem = 80
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
        docker :  "gcr.io/nygc-compbio/novosort:v1.03.01"
        # Per clinical team novosort runs significantly faster with SSD
        disks: "local-disk " + diskSize + " SSD"
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
        docker : "gcr.io/nygc-public/samtools:1.9.1"
        disks: "local-disk " + diskSize + " HDD"
    }
}

task Bqsr38 {
    input {
        # command
        Bam mergedDedupBam
        IndexedReference referenceFa
        IndexedTable callRegions
        String sampleId
        String recalGrpPath = "~{sampleId}.recal_data.grp"
        IndexedVcf MillsAnd1000G
        IndexedVcf Indels
        IndexedVcf dbsnp
        # resources
        Int mem = 36
        Int cpu = 2
        Int diskSize
   }

    command {
        gatk \
        BaseRecalibrator \
        --java-options "-XX:ParallelGCThreads=1" \
        -L ~{callRegions.table} \
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

task PrintReads {
    input {
        # command
        Bam mergedDedupBam
        File recalGrp
        IndexedReference referenceFa
        String sampleId
        String finalBamPath = "~{sampleId}.final.bam"
        # resources
        Int mem = 36
        Int cpu = 2
        Int diskSize
    }

    command {
        gatk \
        ApplyBQSR \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
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
        memory : mem
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
        disks: "local-disk " + diskSize + " HDD"
     }
}
