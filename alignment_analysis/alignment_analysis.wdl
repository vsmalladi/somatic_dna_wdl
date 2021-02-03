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

task GemSelect {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String chr6Contigs
        String sampleId
        Bam finalBam
        File hlaGem
        String r1FilePath = "~{sampleId}.first_pair"
        String r2FilePath = "~{sampleId}.second_pair"
        Float maxMismatches = 0.04
        String alignmentHistoPath = "~{sampleId}.alignment.pdf"
        String r1MappedFastqPath = "~{sampleId}.R1_mapped.fastq"
        String r2MappedFastqPath = "~{sampleId}.R2_mapped.fastq"
    }

    command {
        set -e -o pipefail
        samtools view \
        --threads ~{threads} \
        -h \
        -f 1 \
        ~{finalBam.bam} \
        ~{chr6Contigs} \
        | note_pair.py \
        ~{r1FilePath} \
        ~{r2FilePath} \
        | samtools fastq \
        --threads ~{threads} \
        - \
        | gem-mapper \
        -T ~{threads} \
        --verbose \
        -I ~{hlaGem} \
        -m ~{maxMismatches} \
        -e ~{maxMismatches} \
        --mismatch-alphabet ATCGN \
        --fast-mapping \
        -q ignore \
        | describe_alignments.py \
        ~{alignmentHistoPath} \
        | gem_to_fastq.py \
        ~{r1MappedFastqPath} \
        ~{r2MappedFastqPath}
    }

    output {
        File r2File = "~{r1FilePath}"
        Fastqs r2MappedFastq = "~{r2MappedFastqPath}"
        File r1File = "~{r2FilePath}"
        Fastqs r1MappedFastq = "~{r1MappedFastqPath}"
        File alignmentHisto = "~{sampleId}.alignment.pdf"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task LookUpMates {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String r1UnmappedFilePath = "~{sampleId}.first_pair_unmapped"
        String r2UnmappedFilePath = "~{sampleId}.second_pair_unmapped"
        File r2File
        Fastqs r2MappedFastq
        File r1File
        Fastqs r1MappedFastq
    }

    command {
        look_up_mates.py \
        ~{r1File} \
        ~{r2File} \
        ~{r1MappedFastq} \
        ~{r2MappedFastq} \
        ~{r1UnmappedFilePath} \
        ~{r2UnmappedFilePath}
    }

    output {
        File r1UnmappedFile = "~{r1UnmappedFilePath}"
        File r2UnmappedFile = "~{r2UnmappedFilePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


task GetMates {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String r1UnmappedFastqPath = "~{sampleId}.R1_unmapped.fastq"
        String r2UnmappedFastqPath = "~{sampleId}.R2_unmapped.fastq"
        Bam finalBam
        File r1UnmappedFile
        File r2UnmappedFile
    }

    command {
        set -e -o pipefail
        samtools view \
        --threads ~{threads} \
        -h \
        -f 1 \
        ~{finalBam.bam} \
        | get_mates.py \
        ~{r1UnmappedFile} \
        ~{r2UnmappedFile} \
        | samtools fastq \
        --threads ~{threads} \
        -1 ~{r1UnmappedFastqPath} \
        -2 ~{r2UnmappedFastqPath} \
        -
    }

    output {
        Fastqs r1UnmappedFastq = "~{r1UnmappedFastqPath}"
        Fastqs r2UnmappedFastq = "~{r2UnmappedFastqPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task SortFastqs {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        String fastqPairId
        String sortedFastqPath = "~{sampleId}.~{fastqPairId}_sorted.fastq"
        Fastqs chr6MappedFastq
        Fastqs chr6MappedMatesFastq
    }

    command {
        cat \
        ~{chr6MappedFastq} \
        ~{chr6MappedMatesFastq} \
        | seqkit fx2tab \
        | match_header.py \
        | sort \
        --dictionary-order \
        -k1,1 \
        -T \
        . \
        | seqkit tab2fx \
        > ~{sortedFastqPath}
    }

    output {
        Fastqs sortedFastq = "~{sortedFastqPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task AlignToPanel {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String memoryPerThread
        String sampleId
        String kouramiBamPath = "~{sampleId}.kourami.bam"
        Fastqs r2SortedFastq
        File mergedHlaPanel
        Fastqs r1SortedFastq
    }

    command {
        bwa mem \
        -t ~{threads} \
        ~{mergedHlaPanel} \
        ~{r1SortedFastq} \
        ~{r2SortedFastq} \
        | samtools sort \
        --threads ~{threads} \
        -m ~{memoryPerThread} \
        -O BAM \
        - \
        > ~{kouramiBamPath}
    }

    output {
        Bam kouramiBam = object {
            bam : kouramiBamPath,
            bamIndex : sub(kouramiBamPath, ".bam$", ".bai")
        }
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Kourami {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sampleId
        Bam kouramiBam
    }

    command {
        java \
        -jar kourami.jar \
        -d "./" \
        -o ~{sampleId} \
        ~{kouramiBam.bam}
    }

    output {
        File result = "~{sampleId}.result"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


