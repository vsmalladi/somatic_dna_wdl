version 1.0

import "../wdl_structs.wdl"

task BedtoolsIntersect {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String mantisBedByIntervalListPath
        File mantisBed
        File intervalListBed
    }

    command {
        bedtools \
        intersect \
        -a ~{mantisBed} \
        -b ~{intervalListBed} \
        -u \
        > ~{mantisBedByIntervalListPath} \
    }

    output {
        File mantisBedByIntervalList = "~{mantisBedByIntervalListPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


task MantisExome {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String cores
        String pairName
        String mantisExomeTxtPath = "~{pairName}.mantis.v1.0.4.WGS-targeted.txt"
        String mantisWxsKmerCountsPath = "~{pairName}.mantis.v1.0.4.WGS-targeted.kmer_counts.txt"
        String mantisWxsKmerCountsFinalPath = "~{pairName}.mantis.v1.0.4.WGS-targeted.kmer_counts.txt"
        Bam tumorFinalBam
        Bam normalFinalBam
        File mantisBedByIntervalList
        IndexedReference referenceFa
    }

    command {
        python \
        mantis.py \
        --bedfile ~{mantisBedByIntervalList} \
        --genome ~{referenceFa.fasta} \
        -mrq 20.0 \
        -mlq 25.0 \
        -mlc 20 \
        -mrr 1 \
        --threads ~{cores} \
        -n ~{normalFinalBam.bam} \
        -t ~{tumorFinalBam.bam} \
        -o ~{mantisExomeTxtPath} \
        && \
        mv \
        ~{mantisWxsKmerCountsPath} \
        ~{mantisWxsKmerCountsFinalPath}
    }

    output {
        File mantisWxsKmerCountsFinal = "~{mantisWxsKmerCountsFinalPath}"
        File mantisWxsKmerCountsFiltered = "~{pairName}.mantis.v1.0.4.WGS-targeted.kmer_counts_filtered.txt"
        File mantisWxsStatus = "~{pairName}.mantis.v1.0.4.WGS-targeted.txt.status"
        File mantisExomeTxt = "~{mantisExomeTxtPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task MantisRethreshold {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String mantisStatusFinalPath = "~{pairName}.mantis.v1.0.4.WGS-targeted.status.final.tsv"
        String normal
        File mantisWxsStatus
    }

    command {
        python \
        reset_mantis.py \
        ~{mantisWxsStatus} \
        ~{mantisStatusFinalPath} \
        ~{normal}
    }

    output {
        File mantisStatusFinal = "~{mantisStatusFinalPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}
