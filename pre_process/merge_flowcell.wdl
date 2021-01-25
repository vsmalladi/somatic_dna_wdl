version 1.0

import "../wdl_structs.wdl"

task NovosortMarkDup {
   input {
    # command
    Array[Bam]+ laneBams
    Array[String]+ laneBamsLists
    String sampleId
    String mergedDedupBamPath = "~{sampleId}.merged_dedup.bam"
    # resources
    Int mem
    Int threads
    String dockerImage
   }

    command {
        novosort \
        -c ~{threads} \
        -m ~{mem} \
        -i \
        -o mergedDedupBamPath \
        --forcesort \
        --markDuplicates \
        ${sep=' ' laneBamsLists}
    }
    
    output {
        File mergedDedupBamOnly = mergedDedupBamPath
    }
    
    runtime {
        cpu : threads
        memory : mem
        docker : dockerImage
    }
}

task IndexBam {
   input {
    # command
    File bam
    # resources
    String dockerImage
   }

    command {
        samtools \
        index \
        ~{bam}
    }
    
    output {
        Bam indexedBam = object {
                bam : bam, 
                bamIndex : sub(bam, ".bam$", ".bai")
            }
    }
    
    runtime {
        docker : dockerImage
    }
}

task Bqsr38 {
   input {
    # command
    Bam mergedDedupBam
    IndexedReference indexedReference
    File chromFile
    String sampleId
    String recalGrpPath = "~{sampleId}.recal_data.grp"
    IndexedVcf MillsAnd1000G
    IndexedVcf Indels
    IndexedVcf DbSnp
    # resources
    Int mem
    Int threads
    String dockerImage
   }

    command {
        gatk \
        BaseRecalibrator \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        -L ~{chromFile} \
        -R ~{indexedReference.fasta} \
        -I ~{mergedDedupBam.bam} \
        -O ~{recalGrpPath} \
        --known-sites ~{MillsAnd1000G} \
        --known-sites ~{Indels} \
        --known-sites ~{DbSnp}
    }
    
    output {
        File recalGrp = recalGrpPath
    }
    
    runtime {
        cpu : threads
        memory : mem
        docker : dockerImage
    }
}

task PrintReads {
   input {
    # command
    Bam mergedDedupBam
    File recalGrp
    IndexedReference indexedReference
    String sampleId
    String finalBamPath = "~{sampleId}.final.bam"
    # resources
    Int mem
    Int threads
    String dockerImage
   }

    command {
        gatk \
        ApplyBQSR \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
        -R ~{indexedReference.fasta} \
        -I ~{mergedDedupBam.bam} \
        -O ~{finalBamPath} \
        --bqsr-recal-file ~{recalGrp}
    }
    
    output {
        Bam finalBam = object {
                bam : finalBamPath, 
                bamIndex : sub(finalBamPath, ".bam$", ".bai")
            }
    }
    
    runtime {
        cpu : threads
        memory : mem
        docker : dockerImage
    }
}
