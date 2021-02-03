version 1.0

import "../wdl_structs.wdl"

task MultiMerge {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String flipAD = "'s/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./'"
        String flipBackAD = "'s/##FORMAT=<ID=AD,Number=./##FORMAT=<ID=AD,Number=R/'"
        String pairName
        String multiallelicVcfPath = "~{pairName}.haplotypecaller.v3.5.0.multiallelic.vcf"
        File finalGermlineVcf
    }

    command {
        sed ~{flipAD} ~{finalGermlineVcf} \
        | bcftools \
        norm \
        -O v \
        -m+both - \
        -o - \
        | sed ~{flipBackAD} \
        > ~{multiallelicVcfPath}
    }

    output {
        File multiallelicVcf = "~{multiallelicVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task FilterForHetSnps {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sellectionString = "'vc.getGenotype(\"NORMAL\").isHet()'"
        String pairName
        String hetVcfPath = "~{pairName}.haplotypecaller.v3.5.0.het.vcf"
        IndexedReference referenceFa
        File multiallelicVcf
    }

    command {
        gatk \
        SelectVariants \
        --java-options "-Xmx20g -XX:ParallelGCThreads=2" \
        -restrict-alleles-to BIALLELIC \
        -select-type SNP \
        -select ~{sellectionString} \
        -R ~{referenceFa.fasta} \
        --exclude-filtered \
        -V ~{multiallelicVcf} \
        -O ~{hetVcfPath}
    }

    output {
        File hetVcf = "~{hetVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task FilterBaf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String knownHetVcfPath = "~{pairName}.haplotypecaller.v3.5.0.known.het.vcf"
        File hetVcf
    }

    command {
        python2.7 \
        filter_baf.py \
        ~{hetVcf} \
        ~{knownHetVcfPath}
    }

    output {
        File knownHetVcf = "~{knownHetVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task AlleleCounts {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String alleleCountsTxtPath = "~{pairName}haplotypecaller.v3.5.0.alleles.txt"
        IndexedReference referenceFa
        Bam normalFinalBam
        File knownHetVcf
        Bam tumorFinalBam
    }

    command {
        python2.7 \
        parse_bam_generate_features_v3.py \
        --tumor_bam ~{tumorFinalBam.bam} \
        --normal_bam ~{normalFinalBam.bam} \
        --vcf ~{knownHetVcf} \
        --output ~{alleleCountsTxtPath} \
        --reference ~{referenceFa.fasta}
    }

    output {
        File alleleCountsTxt = "~{alleleCountsTxtPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task CalcBaf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String bafTxtPath = "~{pairName}.haplotypecaller.v3.5.0.baf.txt"
        File alleleCountsTxt
    }

    command {
        python2.7 \
        calc_baf.py \
        ~{alleleCountsTxt} \
        ~{bafTxtPath}
    }

    output {
        File bafTxt = "~{bafTxtPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


