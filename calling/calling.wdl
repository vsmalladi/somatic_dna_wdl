version 1.0

import "../wdl_structs.wdl"

# General tasks

task Gatk4MergeSortVcf {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String pairName
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
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task ReorderVcfColumns{
    input {
        Int threads
        Int memory_gb
        String dockerImage
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
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task AddVcfCommand{
    input {
        Int threads
        Int memory_gb
        String dockerImage
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
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

# caller specific tasks

task MantaWgs{
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String intHVmem
        IndexedReference referenceFa
        Bam normalFinalBam
        File callRegions
        Bam tumorFinalBam
    }

    command {
        mkdir "output" \
        && \
        configManta.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions} \
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
                vcfIndex : "output/results/variants/candidateSmallIndels.vcf.gz.tbi"
            }
        IndexedVcf diploidSV = object {
                vcf : "output/results/variants/diploidSV.vcf.gz", 
                vcfIndex : "output/results/variants/diploidSV.vcf.gz.tbi"
            }
        IndexedVcf somaticSV = object {
                vcf : "output/results/variants/somaticSV.vcf.gz", 
                vcfIndex : "output/results/variants/somaticSV.vcf.gz.tbi"
            }
        IndexedVcf candidateSV = object {
                vcf : "output/results/variants/candidateSV.vcf.gz", 
                vcfIndex : "output/results/variants/candidateSV.vcf.gz.tbi"
            }
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task FilterNonpass {
    input {
        Int threads
        Int memory_gb
        String dockerImage
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
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task Strelka2 {
    input {
        Int threads
        Int memory_gb
        String dockerImage
        String intHVmem
        IndexedReference referenceFa
        Bam normalFinalBam
        File callRegions
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
        --callRegions ~{callRegions} \
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
                vcfIndex : "output/results/variants/somatic.snvs.vcf.gz.tbi"
            }
        IndexedVcf strelka2Indels = object {
                vcf : "output/results/variants/somatic.indels.vcf.gz", 
                vcfIndex : "output/results/variants/somatic.indels.vcf.gz.tbi"
            }
    }

    runtime {
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task LancetExome {
    input {
        Int threads
        Int memory_gb
        String dockerImage
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
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task Mutect2Wgs {
    input {
        Int threads
        Int memory_gb
        String dockerImage
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
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

task Mutect2Filter {
    input {
        Int threads
        Int memory_gb
        String dockerImage
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
        cpu : threads
        memory : memory_gb + "GB"
        docker : dockerImage
    }
}

