version 1.0

import "../wdl_structs.wdl"

task GetSampleName {
    input {
        File finalBam
        File finalBai
        String sampleIdPath = "sampleId.txt"
        Int memoryGb = 1
        Int diskSize
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
            /gatk/gatk \
            --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            GetSampleName \
            -I ~{finalBam} \
            --read-index ~{finalBai} \
            -O ~{sampleIdPath}
    }

    output {
        String bamSampleId = read_string("~{sampleIdPath}")
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
    
    parameter_meta {
        finalBam: {
            localization_optional: true
        }
        
        finalBai: {
            localization_optional: true
        }
    }
}

task UpdateBamSampleName {
    input {
        Bam finalBam
        String sampleId
        String outputPrefix = "~{sampleId}"
        String headerPath = "~{outputPrefix}.reheader.txt"
        String reheaderBamPath = "~{outputPrefix}.reheader.bam"
        String bamIndexPath = "~{outputPrefix}.reheader.bai"
        # resources
        Int diskSize
        Int threads = 4
        Int memoryGb = 8
    }

    command {
    
        set -e -o pipefail
        
        samtools \
        view -H \
        ~{finalBam.bam} \
        | sed "/^@RG/ s/SM:\S+/SM:~{sampleId}/" \
        > ~{headerPath}
        
        samtools \
        reheader \
        ~{headerPath} \
        ~{finalBam.bam} \
        > ~{reheaderBamPath}
        
        samtools \
        index \
        -@ ~{threads} \
        ~{reheaderBamPath} \
        ~{bamIndexPath}
    }

    output {
        Bam reheaderBam = object {
            bam : reheaderBamPath,
            bamIndex : bamIndexPath
        }
    }

    runtime {
        docker : "gcr.io/nygc-public/samtools@sha256:32f29fcd7af01b3941e6f93095e8d899741e81b50bcc838329bd8df43e120cc3"
        disks: "local-disk " + diskSize + " LOCAL"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
    }
}



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
        docker : "gcr.io/nygc-comp-s-fd4e/somatic_dna_tool@sha256:2d61578c9bc8e5ef4a4b25c1e6a32883c67e0fcff6259e43c7da91dbc9c100d7"
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
        Int threads = 8
        Int bwaThreads = 4
        Int samtoolsSortThreads = 4
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
        String resultPrefix = "~{sampleId}.kourami"
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        set -e -o pipefail
        
        java \
        -Xmx~{jvmHeap}m -XX:ParallelGCThreads=4 \
        -jar /Kourami.jar \
        -d /kourami-0.9.6/db/ \
        -o ~{resultPrefix} \
        ~{kouramiBam}
        
    }

    output {
        File result = "~{resultPrefix}.result"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/kourami@sha256:d4b906b979c24ee4669fdbf7ee1dfbdeb5c89d0e34b4b4aaf21ee070e988d74b"
    }
}

task Angsd {
    input {
        Bam normalFinalBam
        File fastNgsAdmixSites
        File fastNgsAdmixSitesBin
        File fastNgsAdmixSitesIdx
        File fastNgsAdmixChroms
        Int threads
        String outprefix
        Int memoryGb = 25
        Int diskSize = ceil(size(normalFinalBam.bam, "GB")) + 30
    }

    command {
        set -exo pipefail

        angsd \
            -i ~{normalFinalBam.bam} \
            -GL 2 \
            -rf ~{fastNgsAdmixChroms} \
            -sites ~{fastNgsAdmixSites} \
            -doMajorMinor 3 \
            -doGlf 2 \
            -minMapQ 30 \
            -minQ 20 \
            -doDepth 1 \
            -doCounts 1 \
            -nThreads ~{threads} \
            -out ~{outprefix}
    }

    runtime {
        memory: memoryGb + "G"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        docker: "gcr.io/nygc-public/angsd@sha256:cd13820de0bc8d400c3e3ff96be6b885b6d3289d53fda56de30bd08508a0bac7"
        disks: "local-disk " + diskSize + " HDD"
    }

    output {
        File beagleFile = "${outprefix}.beagle.gz"
        File beagleLog = "${outprefix}.arg"
        File beagleDepth = "${outprefix}.depthGlobal"
        File beagleSample = "${outprefix}.depthSample"
    }
}

task FastNgsAdmix {
    input {
        File beagleFile
        File fastNgsAdmixRef
        File fastNgsAdmixNind
        String outprefix
        Int memoryGb = 15
        Int threads = 1
        Int diskSize = 30
    }

    command {
        set -exo pipefail

        fastNGSadmix  \
            -likes ~{beagleFile} \
            -fname ~{fastNgsAdmixRef} \
            -Nname ~{fastNgsAdmixNind} \
            -whichPops all \
            -out  ~{outprefix}
    }

    runtime {
        memory: memoryGb + "G"
        mem: memoryGb + "G"
        cpu: threads
        cpus: threads
        docker: "gcr.io/nygc-public/fastngsadmix@sha256:f0a336e9f193ab1b4f1484cbec56e1abef063913102d3126f4b3a6ed7784d7f1"
        disks: "local-disk " + diskSize + " HDD"
    }

    output {
        File fastNgsAdmixQopt = "${outprefix}.qopt"
        File fastNgsAdmixLog = "${outprefix}.log"
    }
}
