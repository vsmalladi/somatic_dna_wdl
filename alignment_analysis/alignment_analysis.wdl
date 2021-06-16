version 1.0

import "../wdl_structs.wdl"

task BedtoolsIntersect {
    input {
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
        docker : "gcr.io/nygc-public/bedtools:v2.26.0"
    }
}


task MantisExome {
    input {        
        String pairName
        String mantisExomeTxtPath = "~{pairName}.mantis.v1.0.4.WGS-targeted.txt"
        String mantisWxsKmerCountsPath = "~{pairName}.mantis.v1.0.4.WGS-targeted.kmer_counts.txt"
        
        Bam tumorFinalBam
        Bam normalFinalBam
        String tumorFinalBamPath = basename(tumorFinalBam.bam)
        String tumorFinalBamIndexPath = basename(tumorFinalBam.bamIndex)
        String normalFinalBamPath = basename(normalFinalBam.bam)
        String normalFinalBamIndexPath = basename(normalFinalBam.bamIndex)
        
        File mantisBedByIntervalList
        IndexedReference referenceFa
        Int threads = 16
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
        File mantisWxsKmerCountsFiltered = "~{pairName}.mantis.v1.0.4.WGS-targeted.kmer_counts_filtered.txt"
        File mantisWxsStatus = "~{pairName}.mantis.v1.0.4.WGS-targeted.txt.status"
        File mantisExomeTxt = "~{mantisExomeTxtPath}"
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/mantis:1.0.4"
    }
}

task MantisRethreshold {
    input {
        String pairName
        String mantisStatusFinalPath = "~{pairName}.mantis.v1.0.4.WGS-targeted.status.final.tsv"
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
        docker : "gcr.io/nygc-internal-tools/somatic_tools:1.0.2"
    }
}

task GetChr6Contigs {
    input  {
        Bam finalBam
        Int diskSize
        Int memoryGb = 2
    }
    
    command {
        /lookup_contigs.py ~{finalBam.bam}
    }
    
    output {
        String chr6Contigs = read_string(stdout())
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-internal-tools/hla_prep:3.3.0"
    }
}

task GemSelect {
    input {
        Int threads = 48
        Int samtoolsThreads = 16
        Int gemThreads = 16
        Int memoryGb = 12
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
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-internal-tools/hla_prep:3.3.0"
    }
}

task LookUpMates {
    input {
        Int memoryGb = 12
        Int diskSize = 10
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
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-internal-tools/hla_prep:3.3.0"
    }
}


task GetMates {
    input {
        Int threads = 32
        Int samtoolsThreads = 16
        Int memoryGb = 12
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
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/hla_prep:3.3.0"
    }
}

task SortFastqs {
    input {
        Int memoryGb = 12
        Int diskSize = 10
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
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-internal-tools/hla_prep:3.3.0"
    }
}

task AlignToPanel {
    input {
        Int threads = 82
        Int bwaThreads = 80
        Int samtoolsSortThreads = 2
        Int memoryGb = 48
        Int diskSize = 10
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
        > test.output.sam
        
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
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bwa-kit:0.7.15"
    }
}

task Kourami {
    input {
        Int threads = 16
        Int memoryGb = 8
        String sampleId
        File kouramiBam
    }

    command {
        java \
        -jar /Kourami.jar \
        -d /kourami-0.9.6/db/ \
        -o ~{sampleId} \
        ~{kouramiBam}
    }

    output {
        File result = "~{sampleId}.result"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/kourami:v0.9.6"
    }
}


