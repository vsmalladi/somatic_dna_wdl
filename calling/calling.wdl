version 1.0

import "../wdl_structs.wdl"

# General tasks

task Gatk4MergeSortVcf {
    input {
        Int diskSize
        Int memoryGb
        String sortedVcfPath
        Array[File] tempChromVcfs
        IndexedReference referenceFa
    }

    command {
        gatk \
        SortVcf \
        --java-options "-Xmx8196m -XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempChromVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
                vcf : "~{sortedVcfPath}",
                vcfIndex : "~{sortedVcfPath}.idx"
            }
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    }
}

task ReorderVcfColumns {
    input {
        Int diskSize
        Int memoryGb
        String normal
        String tumor
        File rawVcf
        String orderedVcfPath
    }

    command {
        python2.7 \
        reorder_vcf.py \
        ~{rawVcf} \
        ~{orderedVcfPath} \
        ~{normal} ~{tumor}
    }

    output {
        File orderedVcf = "~{orderedVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.2"
    }
}

task AddVcfCommand {
    input {
        Int diskSize
        Int memoryGb
        File inVcf
        String outVcfPath = sub(inVcf, ".vcf$", "_w_command.vcf")
        File jsonLog
    }

    command {
        python2.7 \
        add_command.py \
        ~{inVcf} \
        ~{outVcfPath} \
        ~{jsonLog}
    }

    output {
        File outVcf = "~{outVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.2"
    }
}

# caller specific tasks

task MantaWgs {
    input {
        Int threads = 8
        Int memoryGb
        Int diskSize
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        Bam normalFinalBam
        IndexedTable callRegions
        Bam tumorFinalBam
    }

    command {
        mkdir "output" \
        && \
        configManta.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --runDir "output" \
        && \
        "output/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf candidateSmallIndels = object {
                vcf : "output/results/variants/candidateSmallIndels.vcf.gz", 
                index : "output/results/variants/candidateSmallIndels.vcf.gz.tbi"
            }
        IndexedVcf diploidSV = object {
                vcf : "output/results/variants/diploidSV.vcf.gz", 
                index : "output/results/variants/diploidSV.vcf.gz.tbi"
            }
        IndexedVcf somaticSV = object {
                vcf : "output/results/variants/somaticSV.vcf.gz", 
                index : "output/results/variants/somaticSV.vcf.gz.tbi"
            }
        IndexedVcf candidateSV = object {
                vcf : "output/results/variants/candidateSV.vcf.gz", 
                index : "output/results/variants/candidateSV.vcf.gz.tbi"
            }
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " SSD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/manta:1.4.0"
    }
}

task FilterNonpass {
    input {
        Int threads = 4
        Int memoryGb = 8
        Int diskSize
        String pairName
        String outVcfPath = "~{pairName}.manta.v1.4.0.filtered.vcf"
        IndexedReference referenceFa
        File vcf
    }

    command {
        gatk \
        SelectVariants \
        --java-options "-Xmx8g -XX:ParallelGCThreads=4" \
        -R ~{referenceFa.fasta} \
        -V ~{vcf} \
        -O ~{outVcfPath} \
        --exclude-filtered
    }

    output {
        File outVcf = "~{outVcfPath}"
    }

    runtime {
        cpu : 4
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    }
}

task Strelka2 {
    input {
        Int threads
        Int memoryGb
        Int diskSize
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        Bam normalFinalBam
        IndexedTable callRegions
        IndexedVcf candidateSmallIndels
        Bam tumorFinalBam
    }

    command {
        mkdir "output" \
        && \
        configureStrelkaSomaticWorkflow.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --indelCandidates ~{candidateSmallIndels.vcf} \
        --config "configureStrelkaSomaticWorkflow.py.ini" \
        --runDir "output" \
        && \
        "output/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf strelka2Snvs = object {
                vcf : "output/results/variants/somatic.snvs.vcf.gz", 
                index : "output/results/variants/somatic.snvs.vcf.gz.tbi"
            }
        IndexedVcf strelka2Indels = object {
                vcf : "output/results/variants/somatic.indels.vcf.gz", 
                index : "output/results/variants/somatic.indels.vcf.gz.tbi"
            }
    }

    runtime {        
        cpu : threads
        disks: "local-disk " + diskSize + " SSD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/strelka:v2.9.3"
    }
}

task LancetExome {
    input {
        Int threads
        Int diskSize
        Int memoryGb
        String pairName
        String chrom
        String lancetChromVcfPath = "~{pairName}_~{chrom}.lancet.v1.0.7.vcf"
        IndexedReference referenceFa
        File chromBed
        Bam normalFinalBam
        Bam tumorFinalBam
    }

    command {
        lancet \
        --normal ~{normalFinalBam.bam} \
        --tumor ~{tumorFinalBam.bam} \
        --bed ~{chromBed} \
        --ref ~{referenceFa.fasta} \
        --min-k 11 \
        --low-cov 1 \
        --min-phred-fisher 5 \
        --min-strand-bias 1 \
        --min-alt-count-tumor 3 \
        --min-vaf-tumor 0.04 \
        --num-threads ~{threads} \
        > ~{lancetChromVcfPath}
    }

    output {
        File lancetChromVcf = "~{lancetChromVcfPath}"
    }

    runtime {
        cpu : threads        
        disks: "local-disk " + diskSize + " SSD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet:v1.0.7"
    }
}

task Mutect2Wgs {
    input {
        Int memoryGb
        Int diskSize
        String chrom
        String tumor
        String normal
        String pairName
        String mutect2ChromRawVcfPath = "~{pairName}_~{chrom}.mutect2.v4.0.5.1.raw.vcf"
        IndexedReference referenceFa
        Bam normalFinalBam
        Bam tumorFinalBam
    }

    command {
        gatk \
        Mutect2 \
        --java-options "-Xmx8196m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -L ~{chrom} \
        -I ~{tumorFinalBam.bam} \
        -I ~{normalFinalBam.bam} \
        -tumor ~{tumor} \
        -normal ~{normal} \
        -O ~{mutect2ChromRawVcfPath}
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
    }

    runtime {
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.0.5.1"
        disks: "local-disk " + diskSize + " SSD"
    }
}

task Mutect2Filter {
    input {
        Int memoryGb
        Int diskSize
        String pairName
        String chrom
        String mutect2ChromVcfPath = "~{pairName}_~{chrom}.mutect2.v4.0.5.1.vcf"
        IndexedReference referenceFa
        File mutect2ChromRawVcf
    }

    command {
        gatk \
        FilterMutectCalls \
        --java-options "-Xmx8196m -XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -V ~{mutect2ChromRawVcf} \
        -O ~{mutect2ChromVcfPath}
    }

    output {
        File mutect2ChromVcf = "~{mutect2ChromVcfPath}"
    }

    runtime {
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.0.5.1"
        disks: "local-disk " + diskSize + " HDD"
    }
}

task SvabaWgs {
    input {
        Int threads
        Int memoryGb
        String pairName
        BwaReference bwaReference
        Bam normalFinalBam
        IndexedVcf dbsnp
        Bam tumorFinalBam
        Int diskSize
    }

    command {
        svaba \
        run \
        -t ~{tumorFinalBam.bam} \
        -n ~{normalFinalBam.bam} \
        -p ~{threads} \
        -D ~{dbsnp.vcf} \
        -a ~{pairName} \
        -G ~{bwaReference.fasta} \
        -z on
    }

    output {
        Array[File] svabaInternalInput = ["~{pairName}.svaba.unfiltered.somatic.indel.vcf.gz",
                        "~{pairName}.svaba.unfiltered.germline.indel.vcf.gz",
                        "~{pairName}.bps.txt.gz",
                        "~{pairName}.log",
                        "~{pairName}.svaba.germline.indel.vcf.gz",
                        "~{pairName}.contigs.bam",
                        "~{pairName}.discordant.txt.gz",
                        "~{pairName}.svaba.unfiltered.germline.sv.vcf.gz",
                        "~{pairName}.alignments.txt.gz",
                        "~{pairName}.svaba.germline.sv.vcf.gz",
                        "~{pairName}.svaba.unfiltered.somatic.sv.vcf.gz"
                        ]
        File svabaIndelGz = "~{pairName}.svaba.somatic.indel.vcf.gz"
        File svabaGz = "~{pairName}.svaba.somatic.sv.vcf.gz"
    }

    runtime {
        cpu : threads        
        disks: "local-disk " + diskSize + " SSD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/svaba:1.1.0"
    }
}

task UniqReads {
    input {
        Int threads
        Int memoryGb
        String tumor
        String normal
        String dockerImage
        String tempTumorPrefix = "~{tumor}/~{tumor}_"
        String tempNormalPrefix = "~{normal}/~{normal}"
        Array[String] tempNormalSeqsPaths
        Array[String] tempTumorSeqsPaths
        Bam tumorFinalBam
        Bam normalFinalBam
    }

    command {
        mkdir ~{tumor} \
        && \
        mkdir ~{normal} \
        && \
        samtools \
        view \
        -U "BWA,~{tempTumorPrefix},N,N" \
        ~{tumorFinalBam.bam} \
        && \
        samtools \
        view \
        -U "BWA,~{tempNormalPrefix},N,N" \
        ~{normalFinalBam.bam}
    }
    
    output {
        Array[File] tempNormalSeqs = tempNormalSeqsPaths
        Array[File] tempTumorSeqs = tempTumorSeqsPaths
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


task Bicseq2Norm {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        Int readLength
        Int medianInsertSize
        String GCvsRDPath = "~{sampleId}/~{sampleId}.GCvsRD.pdf"
        String paramsPath
        String sampleId
        Array[File] tempSeqs
        Array[File] uniqCoords
        Array[String] tempNormPaths
        File configFile
        Bam FinalBam
        Array[File] chromFastas
    }

    command {
        mkdir ~{sampleId} \
        && \
        perl \
        BICseq2-norm.pl \
        -l=~{readLength} \
        -s=~{medianInsertSize} \
        -fig=~{GCvsRDPath} \
        -tmp=~{sampleId} \
        ~{configFile} \
        ~{paramsPath}
    }

    output {
        Array[File] tempNorm = tempNormPaths
        File GCvsRD = "~{GCvsRDPath}"
        File params = paramsPath
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


task Bicseq2Wgs {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String bicseq2PngPath = "~{pairName}.bicseq2.v0.2.6.png"
        String bicseq2Path = "~{pairName}.bicseq2.v0.2.6.txt"
        Array[File] tempTumorNorms
        Array[File] tempNormalNorms
        File segConfigFile
    }

    command {
        perl \
        NBICseq-seg.pl \
        --control \
        --fig=~{bicseq2PngPath} \
        --title=~{pairName} \
        --lambda=4 \
        ~{segConfigFile} \
        ~{bicseq2Path}
    }

    output {
        File bicseq2Png = "~{pairName}.bicseq2.v0.2.6.png"
        File bicseq2 = "~{pairName}.bicseq2.v0.2.6.txt"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}



