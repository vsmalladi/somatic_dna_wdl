version 1.0

import "../wdl_structs.wdl"

task BedtoolsIntersect {
    input {
        String mantisBedByIntervalListPath
        File mantisBed
        File intervalListBed
        Int memoryGb = 1
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
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bedtools@sha256:9e737f5c96c00cf3b813d419d7a7b474c4013c9aa9dfe704eb36417570c6474e"
    }
}


task MantisExome {
    input {
        String pairName
        String mantisExomeTxtPath = "~{pairName}.mantis.WGS-targeted.txt"
        String mantisWxsKmerCountsPath = "~{pairName}.mantis.WGS-targeted.kmer_counts.txt"
        
        Bam tumorFinalBam
        Bam normalFinalBam
        String tumorFinalBamPath = basename(tumorFinalBam.bam)
        String tumorFinalBamIndexPath = basename(tumorFinalBam.bamIndex)
        String normalFinalBamPath = basename(normalFinalBam.bam)
        String normalFinalBamIndexPath = basename(normalFinalBam.bamIndex)

        File mantisBedByIntervalList
        IndexedReference referenceFa
        Int threads = 8
        Int memoryGb = 4
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 30

    }

    command {
        set -e -o pipefail

        # make a .bam.bai index available
        # normal
        ln -s \
        ~{normalFinalBam.bam} \
        ~{normalFinalBamPath}

        ln -s \
        ~{normalFinalBam.bamIndex} \
        ~{normalFinalBamIndexPath}

        ln -s \
        ~{normalFinalBamIndexPath} \
        ~{normalFinalBamPath}.bai

        # tumor
        ln -s \
        ~{tumorFinalBam.bam} \
        ~{tumorFinalBamPath}

        ln -s \
        ~{tumorFinalBam.bamIndex} \
        ~{tumorFinalBamIndexPath}

        ln -s \
        ~{tumorFinalBamIndexPath} \
        ~{tumorFinalBamPath}.bai

        ls -thl

        python \
        /MANTIS-1.0.4/mantis.py \
        --bedfile ~{mantisBedByIntervalList} \
        --genome ~{referenceFa.fasta} \
        -mrq 20.0 \
        -mlq 25.0 \
        -mlc 20 \
        -mrr 1 \
        --threads ~{threads} \
        -n ~{normalFinalBamPath} \
        -t ~{tumorFinalBamPath} \
        -o ~{mantisExomeTxtPath}
    }

    output {
        File mantisWxsKmerCountsFinal = "~{mantisWxsKmerCountsPath}"
        File mantisWxsKmerCountsFiltered = "~{pairName}.mantis.WGS-targeted.kmer_counts_filtered.txt"
        File mantisWxsStatus = "~{pairName}.mantis.WGS-targeted.txt.status"
        File mantisExomeTxt = "~{mantisExomeTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/mantis@sha256:9cf1311c5198b8fa5fecff387a50dfa9408707f7b914a99dc548c6eb14f42c19"
    }
}

task MantisRethreshold {
    input {
        String pairName
        String mantisStatusFinalPath = "~{pairName}.mantis.WGS-targeted.status.final.tsv"
        String normal
        File mantisWxsStatus
    }

    command {
        python \
        /reset_mantis.py \
        ~{mantisWxsStatus} \
        ~{mantisStatusFinalPath} \
        ~{normal}
    }

    output {
        File mantisStatusFinal = "~{mantisStatusFinalPath}"
    }

    runtime {
        docker : "gcr.io/nygc-public/somatic_tools@sha256:9ae77f7d96a3c100319cf0fac2429f8f84301003480b7b7eb72994ca9f358512"
    }
}

task GetChr6Contigs {
    input  {
        IndexedReference referenceFa
        Int diskSize
        Int memoryGb = 2
    }

    command {
        /lookup_contigs.py ~{referenceFa.fasta}
    }

    output {
        String chr6Contigs = read_string(stdout())
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
    }
}

task GemSelect {
    input {
        Int threads = 8
        Int samtoolsThreads = 8
        Int gemThreads = 8
        Int memoryGb = 4
        Int diskSize
        String sampleId
        String chr6Contigs
        Bam finalBam
        File kouramiFastaGem1Index
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
        --threads ~{samtoolsThreads} \
        -h \
        -f 1 \
        ~{finalBam.bam} \
        ~{chr6Contigs} \
        | /note_pair.py \
        ~{r1FilePath} \
        ~{r2FilePath} \
        | samtools fastq \
        --threads ~{samtoolsThreads} \
        - \
        | gem-mapper \
        -T ~{gemThreads} \
        --verbose \
        -I ~{kouramiFastaGem1Index} \
        -m ~{maxMismatches} \
        -e ~{maxMismatches} \
        --mismatch-alphabet ATCGN \
        --fast-mapping \
        -q ignore \
        | /describe_alignments.py \
        ~{alignmentHistoPath} \
        | /gem_to_fastq.py \
        ~{r1MappedFastqPath} \
        ~{r2MappedFastqPath}
    }

    output {
        File r2File = "~{r1FilePath}"
        File r2MappedFastq = "~{r2MappedFastqPath}"
        File r1File = "~{r2FilePath}"
        File r1MappedFastq = "~{r1MappedFastqPath}"
        File alignmentHisto = "~{sampleId}.alignment.pdf"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
    }
}

task LookUpMates {
    input {
        Int memoryGb = 2
        Int diskSize = 4
        String sampleId
        String r1UnmappedFilePath = "~{sampleId}.first_pair_unmapped"
        String r2UnmappedFilePath = "~{sampleId}.second_pair_unmapped"
        File r2File
        File r2MappedFastq
        File r1File
        File r1MappedFastq

    }

    command {
        /look_up_mates.py \
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
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
    }
}


task GetMates {
    input {
        Int threads = 8
        Int samtoolsThreads = 4
        Int memoryGb = 2
        Int diskSize
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
        --threads ~{samtoolsThreads} \
        -h \
        -f 1 \
        ~{finalBam.bam} \
        | /get_mates.py \
        ~{r1UnmappedFile} \
        ~{r2UnmappedFile} \
        | samtools fastq \
        --threads ~{samtoolsThreads} \
        -1 ~{r1UnmappedFastqPath} \
        -2 ~{r2UnmappedFastqPath} \
        -
    }

    output {
        File r1UnmappedFastq = "~{r1UnmappedFastqPath}"
        File r2UnmappedFastq = "~{r2UnmappedFastqPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
    }
}

task SortFastqs {
    input {
        Int memoryGb = 2
        Int diskSize = 4
        String sampleId
        String fastqPairId
        String sortedFastqPath = "~{sampleId}.~{fastqPairId}_sorted.fastq"
        File chr6MappedFastq
        File chr6MappedMatesFastq
    }

    command {
        set -e -o pipefail

        cat \
        ~{chr6MappedFastq} \
        ~{chr6MappedMatesFastq} \
        | seqkit fx2tab \
        | /match_header.py \
        | sort \
        --dictionary-order \
        -k1,1 \
        -T \
        . \
        | seqkit tab2fx \
        > ~{sortedFastqPath}
    }

    output {
        File sortedFastq = "~{sortedFastqPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/hla_prep@sha256:a490cf449eeb98b997f0dd87ff1c23ff77d724c7c2072b6c44f75a713ecc2d36"
    }
}

task AlignToPanel {
    input {
        Int threads = 82
        Int bwaThreads = 80
        Int samtoolsSortThreads = 2
        Int memoryGb = 4
        Int diskSize = 4
        String sampleId
        String kouramiBamPath = "~{sampleId}.kourami.bam"
        File r2SortedFastq
        # mergedHlaPanel
        BwaReference kouramiReference
        File r1SortedFastq
    }

    command {
        set -e -o pipefail

        bwa mem \
        -t ~{bwaThreads} \
        ~{kouramiReference.fasta} \
        ~{r1SortedFastq} \
        ~{r2SortedFastq} \
        | samtools sort \
        --threads ~{samtoolsSortThreads} \
        -m 10G \
        -o ~{kouramiBamPath}
    }

    output {
        File kouramiBam = "~{kouramiBamPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bwa-kit@sha256:0642151a32fe8f90ece70cde3bd61a03c7421314c37c1de2c0ee5e368d2bfc7a"
    }
}

task Kourami {
    input {
        Int threads = 1
        Int memoryGb = 8
        String sampleId
        File kouramiBam
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        java \
        -Xmx~{jvmHeap}m -XX:ParallelGCThreads=4 \
        -jar /Kourami.jar \
        -d /kourami-0.9.6/db/ \
        -o ~{sampleId} \
        ~{kouramiBam}
    }

    output {
        File result = "~{sampleId}.kourami.result"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/kourami@sha256:d4b906b979c24ee4669fdbf7ee1dfbdeb5c89d0e34b4b4aaf21ee070e988d74b"
    }
}
