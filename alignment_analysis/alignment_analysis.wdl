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
        
        String altTumorIndexPath = sub(basename(tumorFinalBam.bamIndex), ".bai$", ".bam.bai")
        String altNormalIndexPath = sub(basename(normalFinalBam.bamIndex), ".bai$", ".bam.bai")
        File mantisBedByIntervalList
        IndexedReference referenceFa
        Int threads = 16
        Int memoryGb = 4
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 30
        
    }

    command {
        set -e -o pipefail
        
        ln -s \
        ~{normalFinalBam.bamIndex} \
        ~{altNormalIndexPath}
        
        ln -s \
        ~{tumorFinalBam.bamIndex} \
        ~{altTumorIndexPath}
        
        python \
        /MANTIS-1.0.4/mantis.py \
        --bedfile ~{mantisBedByIntervalList} \
        --genome ~{referenceFa.fasta} \
        -mrq 20.0 \
        -mlq 25.0 \
        -mlc 20 \
        -mrr 1 \
        --threads ~{threads} \
        -n ~{normalFinalBam.bam} \
        -t ~{tumorFinalBam.bam} \
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
        File resetMantis = "gs://nygc-comp-s-fd4e-input/reset_mantis.py"
    }

    command {
    
        set -e -o pipefail
    
        chmod 755 ~{resetMantis}
        
        python \
        ~{resetMantis} \
        ~{mantisWxsStatus} \
        ~{mantisStatusFinalPath} \
        ~{normal}
    }

    output {
        File mantisStatusFinal = "~{mantisStatusFinalPath}"
    }

    runtime {
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.2"
    }
}

task GetChr6Contigs {
    input  {
        Bam finalBam
        Int diskSize
        Int memoryGb = 2
    }
    
    command {
        set -e -o pipefail
        
        export CONDA_ALWAYS_YES="true"
        
        conda config --add channels r
        conda config --add channels bioconda    
        
        conda install pysam &> "pysam_install.log"
        
        /lookup_contigs.py ~{finalBam.bam}
    }
    
    output {
        String chr6Contigs = read_string(stdout())
    }
    
    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-internal-tools/hla_prep:1.1.0"
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
        File kouramiFastaGem3Index
        String r1FilePath = "~{sampleId}.first_pair"
        String r2FilePath = "~{sampleId}.second_pair"
        Float maxMismatches = 0.04
        Float alignmentGlobalMinIdentity = 0.80
        String outputFormat = "MAP"
        String alignmentHistoPath = "~{sampleId}.alignment.pdf"
        String r1MappedFastqPath = "~{sampleId}.R1_mapped.fastq"
        String r2MappedFastqPath = "~{sampleId}.R2_mapped.fastq"
        
        File describeAlignments = "gs://nygc-comp-s-fd4e-input/describe_alignments.py"
        File gemToFastq = "gs://nygc-comp-s-fd4e-input/gem_to_fastq.py"
    }

    command {
        set -e -o pipefail
        
        chmod 755 ~{describeAlignments}
        
        chmod 755 ~{gemToFastq}
        
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
        --threads ~{gemThreads} \
        --verbose \
        --index ~{kouramiFastaGem3Index} \
        --alignment-global-min-identity ~{alignmentGlobalMinIdentity} \
        --alignment-max-error ~{maxMismatches} \
        --output-format ~{outputFormat}
        --mapping-mode "fast"
        | ~{describeAlignments} \
        ~{alignmentHistoPath} \
        | ~{gemToFastq} \
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
        docker : "gcr.io/nygc-internal-tools/hla_prep:1.1.0"
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
        set -e -o pipefail
    
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
        docker : "gcr.io/nygc-internal-tools/hla_prep:1.1.0"
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
        docker : "gcr.io/nygc-internal-tools/hla_prep:1.1.0"
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
        File sortedFastq = "~{sortedFastqPath}"
    }

    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-internal-tools/hla_prep:1.1.0"
    }
}

task AlignToPanel {
    input {
        Int threads = 96
        Int bwaThreads = 80
        Int samtoolsThreads = 16
        Int memoryGb = 24
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
        | samtools sort \
        --threads ~{samtoolsThreads} \
        -O BAM \
        - \
        > ~{kouramiBamPath}
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
        -jar kourami.jar \
        -d "./" \
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


