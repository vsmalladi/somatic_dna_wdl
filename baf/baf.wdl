version 1.0

import "../wdl_structs.wdl"

task FilterForHetSnps {
    input {
        String sampleId
        String hetVcfPath = "~{sampleId}.haplotypecaller.gatk.het.vcf"
        String sellectionString = "'vc.getGenotype(\"~{sampleId}\").isHet()'"
        IndexedReference referenceFa
        # require file!
        # marked as optional so that pipeline can be dependent on input that may not
        # pass QC
        # do not run with out input file!
        File? germlineVcf

        Int memoryGb = 24
        Int diskSize = (ceil( size(germlineVcf, "GB") )  * 2 )
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        gatk \
        SelectVariants \
        --java-options "-Xmx~{jvmHeap}m" \
        -restrict-alleles-to BIALLELIC \
        -select-type SNP \
        -select ~{sellectionString} \
        -R ~{referenceFa.fasta} \
        --exclude-filtered \
        -V ~{germlineVcf} \
        -O ~{hetVcfPath}
    }

    output {
        File hetVcf = "~{hetVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

task FilterBaf {
    input {
        String sampleId
        String knownHetVcfPath = "~{sampleId}.haplotypecaller.gatk.known.het.vcf"
        File hetVcf
        Int memoryGb = 4
        Int diskSize = (ceil( size(hetVcf, "GB") )  * 2 )
    }

    command {
        python \
        /filter_baf.py \
        ~{hetVcf} \
        ~{knownHetVcfPath}
    }

    output {
        File knownHetVcf = "~{knownHetVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:20a48e2c422a43ce35e197243bda8dbf06c9a7b3175094524f74f8835cce85b6"
    }
}

task AlleleCounts {
    input {
        String pairName
        String alleleCountsChromTxtPath = "~{pairName}.~{chrom}.haplotypecaller.gatk.alleles.txt"
        IndexedReference referenceFa
        Bam normalFinalBam
        File knownHetVcf
        Bam tumorFinalBam
        String chrom
        Int memoryGb = 4

        Int diskSize = (ceil( size(knownHetVcf, "GB") )  * 2 ) + ceil( size(tumorFinalBam.bam, "GB")) + ceil( size(normalFinalBam.bam, "GB")) + 10
    }

    command {
        python \
        /parse_bam_generate_features.py \
        --tumor_bam ~{tumorFinalBam.bam} \
        --normal_bam ~{normalFinalBam.bam} \
        --vcf ~{knownHetVcf} \
        --output ~{alleleCountsChromTxtPath} \
        --reference ~{referenceFa.fasta} \
        --chrom ~{chrom}
    }

    output {
        File alleleCountsTxt = "~{alleleCountsChromTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:20a48e2c422a43ce35e197243bda8dbf06c9a7b3175094524f74f8835cce85b6"
    }
}

task CalcBaf {
    input {
        String pairName
        String bafTxtPath = "~{pairName}.haplotypecaller.gatk.baf.txt"
        File alleleCountsTxt
        Int memoryGb = 4
        Int diskSize = (ceil( size(alleleCountsTxt, "GB") )  * 2 ) + 1
    }

    command {
        python \
        /calc_baf.py \
        ~{alleleCountsTxt} \
        ~{bafTxtPath}
    }

    output {
        File bafTxt = "~{bafTxtPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:20a48e2c422a43ce35e197243bda8dbf06c9a7b3175094524f74f8835cce85b6"
    }
}
