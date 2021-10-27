version 1.0

import "../wdl_structs.wdl"

# General tasks

task GetInsertSize {
    input {
        File insertSizeMetrics
    }

    command {
        grep -A 1 \
        "MEDIAN_INSERT_SIZE" \
        ~{insertSizeMetrics} \
        | tail -n 1 \
        | cut -f 1
    }

    output {
        Int insertSize = read_int(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
    }
}

task Gatk4MergeSortCompressVcf {
    input {
        Int diskSize = 10
        Int memoryGb = 8
        String sortedVcfPath
        Array[File] tempChromVcfs
        IndexedReference referenceFa
    }

    command {
        gatk \
        SortVcf \
        --java-options "-XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempChromVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedIndexedVcf = object {
                vcf : "~{sortedVcfPath}",
                index : "~{sortedVcfPath}.tbi"
            }
    }

    parameter_meta {
        sortedVcfPath: {
            description: "Output VCF filename for file in compressed format. Must end in .gz. Corresponding index file will end in .tbi",
            category: "other"
        }

        tempChromVcfs: {
            description: "Input VCF filenames",
            category: "required"
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required"
        }
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    }
}

task Gatk4MergeSortVcf {
    input {
        Int diskSize = 10
        Int memoryGb = 8
        String sortedVcfPath
        Array[File] tempChromVcfs
        IndexedReference referenceFa
    }

    command {
        gatk \
        SortVcf \
        --java-options "-XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempChromVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
                vcf : "~{sortedVcfPath}",
                index : "~{sortedVcfPath}.idx"
            }
    }

    parameter_meta {
        sortedVcfPath: {
            description: "Output VCF filename for uncompressed output file. Must not end in .gz. Corresponding index file will end in .idx",
            category: "other"
        }

        tempChromVcfs: {
            description: "Input VCF filenames",
            category: "required"
        }

        referenceFa: {
            description: "Fasta object (the dictionary will be used to order the final VCF)",
            category: "required"
        }
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    }
}

task AddCommandReorderColumnsVcf {
    input {
        Int diskSize
        Int memoryGb
        String normal
        String tumor
        File inVcf
        String orderedVcfPath
        String outVcfPath = sub(sub(basename(inVcf), ".gz$", ""), ".vcf$", "_w_command.vcf")
        File jsonLog
    }

    command {
    
        python \
        /add_command.py \
        ~{inVcf} \
        ~{outVcfPath} \
        ~{jsonLog}
        
        python \
        /reorder_vcf.py \
        ~{outVcfPath} \
        ~{orderedVcfPath} \
        ~{normal} ~{tumor}
    }

    output {
        File orderedVcf = "~{orderedVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
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
        python \
        /reorder_vcf.py \
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
        docker : "gcr.io/nygc-public/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
}

task AddVcfCommand {
    input {
        Int diskSize
        Int memoryGb
        File inVcf
        String outVcfPath = sub(sub(basename(inVcf), ".gz$", ""), ".vcf$", "_w_command.vcf")
        File jsonLog
    }

    command {
        python \
        /add_command.py \
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
        docker : "gcr.io/nygc-public/somatic_tools@sha256:46ab81b8dc09d6f8cf90c81f7d5692f23d73c134df6dbcd5298abde7f414dce3"
    }
}

# caller specific tasks

task MantaWgs {
    input {
        Int threads = 8
        Int memoryGb = 4
        Int diskSize
        String pairName
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        Bam normalFinalBam
        IndexedTable callRegions
        Bam tumorFinalBam
    }

    command {
        set -e -o pipefail

        mkdir ~{pairName}.MantaRaw

        configManta.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --runDir ~{pairName}.MantaRaw

        "~{pairName}.MantaRaw/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf candidateSmallIndels = object {
                vcf : "~{pairName}.MantaRaw/results/variants/candidateSmallIndels.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/candidateSmallIndels.vcf.gz.tbi"
            }
        IndexedVcf diploidSV = object {
                vcf : "~{pairName}.MantaRaw/results/variants/diploidSV.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/diploidSV.vcf.gz.tbi"
            }
        IndexedVcf somaticSV = object {
                vcf : "~{pairName}.MantaRaw/results/variants/somaticSV.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/somaticSV.vcf.gz.tbi"
            }
        IndexedVcf candidateSV = object {
                vcf : "~{pairName}.MantaRaw/results/variants/candidateSV.vcf.gz",
                index : "~{pairName}.MantaRaw/results/variants/candidateSV.vcf.gz.tbi"
            }
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/manta@sha256:e171112cccf6758693b7a8aab80000fe6121d6969ce83e14a4e15fbc5f2f3662"
    }
}

task MantaWgsPon {
    input {
        Int threads = 8
        Int memoryGb = 4
        Int diskSize
        String sampleId
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        IndexedTable callRegions
        Bam tumorFinalBam
    }

    command {
        set -e -o pipefail

        mkdir ~{sampleId}.MantaRaw

        configManta.py \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --runDir ~{sampleId}.MantaRaw

        "~{sampleId}.MantaRaw/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf somaticSV = object {
                vcf : "~{sampleId}.MantaRaw/results/variants/somaticSV.vcf.gz",
                index : "~{sampleId}.MantaRaw/results/variants/somaticSV.vcf.gz.tbi"
            }

    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/manta@sha256:e171112cccf6758693b7a8aab80000fe6121d6969ce83e14a4e15fbc5f2f3662"
    }
}

task FilterNonpass {
    input {
        Int threads = 4
        Int memoryGb = 8
        Int diskSize
        String pairName
        String outVcfPath = "~{pairName}.manta.v1.4.0.filtered.unorder.vcf"
        IndexedReference referenceFa
        File vcf
    }

    command {
        gatk \
        SelectVariants \
        --java-options "-XX:ParallelGCThreads=4" \
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
        String pairName
        String intHVmem = "unlimited"
        IndexedReference referenceFa
        Bam normalFinalBam
        IndexedTable callRegions
        IndexedVcf candidateSmallIndels
        Bam tumorFinalBam
        File configureStrelkaSomaticWorkflow
    }

    command {
        set -e -o pipefail

        mkdir ~{pairName}.Strelka2Raw

        configureStrelkaSomaticWorkflow.py \
        --normalBam ~{normalFinalBam.bam} \
        --tumorBam ~{tumorFinalBam.bam} \
        --referenceFasta ~{referenceFa.fasta} \
        --callRegions ~{callRegions.table} \
        --indelCandidates ~{candidateSmallIndels.vcf} \
        --config ~{configureStrelkaSomaticWorkflow} \
        --runDir ~{pairName}.Strelka2Raw \

        "~{pairName}.Strelka2Raw/runWorkflow.py" \
        --mode local \
        --job ~{threads} \
        --memGb ~{intHVmem}
    }

    output {
        IndexedVcf strelka2Snvs = object {
                vcf : "~{pairName}.Strelka2Raw/results/variants/somatic.snvs.vcf.gz",
                index : "~{pairName}.Strelka2Raw/results/variants/somatic.snvs.vcf.gz.tbi"
            }
        IndexedVcf strelka2Indels = object {
                vcf : "~{pairName}.Strelka2Raw/results/variants/somatic.indels.vcf.gz",
                index : "~{pairName}.Strelka2Raw/results/variants/somatic.indels.vcf.gz.tbi"
            }
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/strelka@sha256:eb71db1fbe25d67c025251823ae5d5e9dbf5a6861a98298b47ecd722dfa5bd14"
    }
}

task LancetWGSRegional {
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
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet@sha256:25169d34b41de9564e03f02ebcbfb4655cf536449592b0bd58773195f9376e61"
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
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet@sha256:25169d34b41de9564e03f02ebcbfb4655cf536449592b0bd58773195f9376e61"
    }
}

task Mutect2Wgs {
    input {
        Int memoryGb = 4
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
        --java-options "-XX:ParallelGCThreads=4" \
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
        disks: "local-disk " + diskSize + " LOCAL"
    }
}

task Mutect2WgsPon {
    input {
        Int memoryGb = 4
        Int diskSize
        String chrom
        String tumor
        String mutect2ChromRawVcfPath = "~{tumor}_~{chrom}.mutect2.v4.0.5.1.raw.vcf"
        IndexedReference referenceFa
        Bam tumorFinalBam
    }

    command {
        gatk \
        Mutect2 \
        --java-options "-XX:ParallelGCThreads=4" \
        --reference ~{referenceFa.fasta} \
        -L ~{chrom} \
        -I ~{tumorFinalBam.bam} \
        -tumor ~{tumor} \
        -O ~{mutect2ChromRawVcfPath}
    }

    output {
        File mutect2ChromRawVcf = "~{mutect2ChromRawVcfPath}"
    }

    runtime {
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.0.5.1"
        disks: "local-disk " + diskSize + " LOCAL"
    }
}

task Mutect2Filter {
    input {
        Int memoryGb = 4
        Int diskSize = 5
        String pairName
        String chrom
        String mutect2ChromVcfPath = "~{pairName}_~{chrom}.mutect2.v4.0.5.1.vcf"
        IndexedReference referenceFa
        File mutect2ChromRawVcf
    }

    command {
        gatk \
        FilterMutectCalls \
        --java-options "-XX:ParallelGCThreads=4" \
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
        Int memoryGb = 16
        String pairName
        IndexedTable callRegions
        BwaReference bwaReference
        Bam normalFinalBam
        File dbsnpIndels
        Bam tumorFinalBam
        Int diskSize
    }

    command {
        svaba \
        run \
        -t ~{tumorFinalBam.bam} \
        -n ~{normalFinalBam.bam} \
        -p ~{threads} \
        --region ~{callRegions.table} \
        -D ~{dbsnpIndels} \
        -a ~{pairName} \
        -G ~{bwaReference.fasta} \
        -z on
    }

    output {
        File svabaRawGermlineIndel = "~{pairName}.svaba.germline.indel.vcf.gz"
        File svabaRawGermlineSv = "~{pairName}.svaba.germline.sv.vcf.gz"
        File svabaIndelGz = "~{pairName}.svaba.somatic.indel.vcf.gz"
        File svabaGz = "~{pairName}.svaba.somatic.sv.vcf.gz"
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/svaba@sha256:48f6bd86e933ca88fd74d8effc66e93eee5b40945ee37612b80d7edaadc567f3"
        preemptible: 0
    }
}

task SvabaWgsPon {
    input {
        Int threads
        Int memoryGb = 16
        String sampleId
        IndexedTable callRegions
        BwaReference bwaReference
        File dbsnpIndels
        Bam tumorFinalBam
        Int diskSize
    }

    command {
        svaba \
        run \
        -t ~{tumorFinalBam.bam} \
        -p ~{threads} \
        --region ~{callRegions.table} \
        -D ~{dbsnpIndels} \
        -a ~{sampleId} \
        -G ~{bwaReference.fasta} \
        -z on
    }

    output {
        File svabaIndelGz = "~{sampleId}.svaba.somatic.indel.vcf.gz"
        File svabaGz = "~{sampleId}.svaba.somatic.sv.vcf.gz"
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/svaba@sha256:48f6bd86e933ca88fd74d8effc66e93eee5b40945ee37612b80d7edaadc567f3"
        preemptible: 1
    }
}

task ReheaderVcf {
    input {
        Array[String] sampleIds
        String outVcfPath
        IndexedReference referenceFa
        File inVcf
        Int memoryGb = 16
        Int diskSize = (ceil( size(inVcf, "GB") )  * 2 ) + 1
    }

    command {
        bcftools \
        reheader \
        --samples ~{sep=" " sampleIds} \
        --output ~{outVcfPath} \
        ~{inVcf}
    }

    output {
        File reheaderedVcf = "~{outVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
    }
}

task UniqReads {
    input {
        Int memoryGb
        String sampleId
        String tempPrefix = "~{sampleId}_"
        Array[String] tempSeqsPaths
        Bam finalBam
        Int diskSize
    }

    command {
        /samtools-0.1.7a_getUnique-0.1.3/samtools \
        view \
        -U "BWA,~{tempPrefix},N,N" \
        ~{finalBam.bam}

    }

    output {
        Array[File] tempSeqs = tempSeqsPaths
    }

    runtime {
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bicseq2@sha256:3d110b672df0385f761fb64fcf63e98685e6d810c5560043efed5ab94961f3a9"
    }
}


task Bicseq2Norm {
    input {
        Int memoryGb
        Int readLength
        Int medianInsertSize
        String GCvsRDPath = "~{sampleId}.GCvsRD.pdf"
        String paramsPath = "~{sampleId}.params.out"
        String sampleId

        Array[File] tempSeqs
        Array[String] tempNormPaths

        File bicseq2ConfigFile
        String configFilePath = "~{sampleId}.bicseq2.config"

        Array[File] uniqCoordsFiles
        Array[File] chromFastasFiles
        Int diskSize = 70
    }

    command {
        set -e -o pipefail

        python3 \
        /bicseq2_config_writer.py \
        --fa-files ~{sep=" " chromFastasFiles} \
        --mappability-files ~{sep=" " uniqCoordsFiles} \
        --temp-seqs ~{sep=" " tempSeqs} \
        --temp-norm-paths ~{sep=" " tempNormPaths} \
        --norm-bicseq2-config ~{bicseq2ConfigFile} \
        --sample-id ~{sampleId} \
        --out-file ~{configFilePath}

        mkdir -p ~{sampleId}

        /NBICseq-norm_v0.2.4/NBICseq-norm.pl \
        -l=~{readLength} \
        -s=~{medianInsertSize} \
        -fig=~{GCvsRDPath} \
        -tmp=~{sampleId} \
        ~{configFilePath} \
        ~{paramsPath}
    }

    output {
        Array[File] tempNorm = tempNormPaths
        File params = paramsPath
    }

    runtime {
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bicseq2@sha256:3d110b672df0385f761fb64fcf63e98685e6d810c5560043efed5ab94961f3a9"
    }
}


task Bicseq2Wgs {
    input {
        Int memoryGb
        String pairName
        String bicseq2PngPath = "~{pairName}.bicseq2.v0.2.6.png"
        String bicseq2Path = "~{pairName}.bicseq2.v0.2.6.txt"
        Array[File] tempTumorNorms
        Array[File] tempNormalNorms
        File bicseq2SegConfigFile
        String segConfigFilePath = "~{pairName}.bicseq2.seg.config"
        Int lambda = 4
        Int diskSize = 10
    }

    command {
        set -e -o pipefail

        python3 \
        /bicseq2_seg_config_writer.py \
        --tumor-norms ~{sep=" " tempTumorNorms} \
        --normal-norms ~{sep=" " tempNormalNorms} \
        --seg-bicseq2-config ~{bicseq2SegConfigFile} \
        --out-file ~{segConfigFilePath} \
        --pair-id ~{pairName}

        perl /NBICseq-seg_v0.7.2/NBICseq-seg.pl \
        --control \
        --fig=~{bicseq2PngPath} \
        --title=~{pairName} \
        --lambda=4 \
        ~{segConfigFilePath} \
        ~{bicseq2Path}
    }

    output {
        File bicseq2Png = "~{pairName}.bicseq2.v0.2.6.png"
        File bicseq2 = "~{pairName}.bicseq2.v0.2.6.txt"
    }

    runtime {
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bicseq2@sha256:3d110b672df0385f761fb64fcf63e98685e6d810c5560043efed5ab94961f3a9"
    }
}


task GridssPreprocess {
    input {
        Int threads
        Int memoryGb = 16
        Int diskSize = 700
        Bam finalBam

        BwaReference bwaReference
        Array[File] gridssAdditionalReference
    }

    String bamBase = basename(finalBam.bam)
    String svBamPath = bamBase + ".gridss.working/" + bamBase + ".sv.bam"
    String cigarMetricsPath = bamBase + ".gridss.working/" + bamBase + ".cigar_metrics"
    String idsvMetricsPath = bamBase + ".gridss.working/" + bamBase + ".idsv_metrics"
    String tagMetricsPath = bamBase + ".gridss.working/" + bamBase + ".tag_metrics"
    String mapqMetricsPath = bamBase + ".gridss.working/" + bamBase + ".mapq_metrics"
    String insertSizeMetricsPath = bamBase + ".gridss.working/" + bamBase + ".insert_size_metrics"

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir

        mkdir -p /scratch/
        working=$( pwd )

        gridss.sh \
        --steps preprocess \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.11.1-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        ~{finalBam.bam} \
        && cat *.log
    }

    output {
        Bam svBam = object {
            bam : svBamPath,
            bamIndex : sub(svBamPath, ".bam$", ".bam.bai")
        }
        File cigarMetrics = "~{cigarMetricsPath}"
        File idsvMetrics = "~{idsvMetricsPath}"
        File tagMetrics = "~{tagMetricsPath}"
        File mapqMetrics = "~{mapqMetricsPath}"
        File insertSizeMetrics = "~{insertSizeMetricsPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " LOCAL"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:0dd813126af95f511f0576045165a8c14fca46ea99fe44110df2f8e6b4febf72"
    }
}

task GridssAssembleChunk {
    input {
        Int threads
        Int memoryGb = 48
        Int diskSize = 700
        String pairName

        String gridssassemblyBamPath = "~{pairName}.gridssassembly.bam"
        Int jobIndex
        Int assembleChunks
        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Bam tumorFinalBam
        Bam normalFinalBam

        Bam normalSvBam
        File normalCigarMetrics
        File normalIdsvMetrics
        File normalTagMetrics
        File normalMapqMetrics
        File normalInsertSizeMetrics

        Bam tumorSvBam
        File tumorCigarMetrics
        File tumorIdsvMetrics
        File tumorTagMetrics
        File tumorMapqMetrics
        File tumorInsertSizeMetrics
    }

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir


        # reposition unnamed input
        bash gridss_arrange.sh \
        --tumorFinalBam ~{tumorFinalBam.bam} \
        --normalFinalBam ~{normalFinalBam.bam} \
        --tumorSvBam ~{tumorSvBam.bam} \
        --normalSvBam ~{normalSvBam.bam}

        working=$( pwd )

        bash gridss.sh \
        --steps assemble \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.11.1-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --assembly ~{gridssassemblyBamPath} \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        --jobindex ~{jobIndex} \
        --jobnodes ~{assembleChunks} \
        ~{normalFinalBam.bam} ~{tumorFinalBam.bam} \
        && cat *.log
    }

    output {
        File downsampled = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.downsampled_~{jobIndex}.bed"
        File excluded = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.excluded_~{jobIndex}.bed"
        File subsetCalled = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.subsetCalled_~{jobIndex}.bed"
    }

    runtime {
        disks: "local-disk " + diskSize + " LOCAL"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:0dd813126af95f511f0576045165a8c14fca46ea99fe44110df2f8e6b4febf72"
    }
}


task GridssAssemble {
    input {
        Int threads
        Int memoryGb = 48
        Int diskSize = 700
        String pairName
        String gridssassemblyBamPath = "~{pairName}.gridssassembly.bam"
        String gridssassemblySvBamPath = "~{pairName}.gridssassembly.bam.gridss.working/~{pairName}.gridssassembly.bam.sv.bam"
        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Bam tumorFinalBam
        Bam normalFinalBam
        Array[File] downsampled
        Array[File] excluded
        Array[File] subsetCalled

        # from earlier
        Bam normalSvBam
        File normalCigarMetrics
        File normalIdsvMetrics
        File normalTagMetrics
        File normalMapqMetrics
        File normalInsertSizeMetrics
        Bam tumorSvBam
        File tumorCigarMetrics
        File tumorIdsvMetrics
        File tumorTagMetrics
        File tumorMapqMetrics
        File tumorInsertSizeMetrics
    }

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir

        # link preprocess results
        # reposition unnamed input
        bash gridss_arrange.sh \
        --tumorFinalBam ~{tumorFinalBam.bam} \
        --normalFinalBam ~{normalFinalBam.bam} \
        --tumorSvBam ~{tumorSvBam.bam} \
        --normalSvBam ~{normalSvBam.bam}

        # link chunk assembly results
        sub_dir=~{gridssassemblyBamPath}.gridss.working/
        mkdir -p $sub_dir
        copied_sub_dir=$( dirname ~{downsampled[0]} )
        ln -s $copied_sub_dir/* $sub_dir

        working=$( pwd )

        bash gridss.sh \
        --steps assemble \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.11.1-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --assembly ~{gridssassemblyBamPath} \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        ~{normalFinalBam.bam} ~{tumorFinalBam.bam} \
        && cat *.log
    }

    output {
        File gridssassemblyBam = gridssassemblyBamPath
        Bam gridssassemblySvBam = object {
            bam : gridssassemblySvBamPath,
            bamIndex : sub(gridssassemblySvBamPath, ".bam$", ".bam.bai")
        }
    }

    runtime {
        disks: "local-disk " + diskSize + " LOCAL"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:0dd813126af95f511f0576045165a8c14fca46ea99fe44110df2f8e6b4febf72"
    }
}

task GridssCalling {
    input {
        Int threads
        Int memoryGb = 48
        Int diskSize = 700
        String pairName
        String gridssassemblyBamPath = "~{pairName}.gridssassembly.bam"
        String gridssUnfilteredVcfPath = "~{pairName}.sv.gridss.v2.10.2.unfiltered.vcf"
        BwaReference bwaReference
        Array[File] gridssAdditionalReference
        Array[File] downsampled
        Array[File] excluded
        Array[File] subsetCalled
        File gridssassemblyBam
        Bam gridssassemblySvBam
        Bam tumorFinalBam
        Bam normalFinalBam

        # required from Preprocessing
        Bam normalSvBam
        File normalCigarMetrics
        File normalIdsvMetrics
        File normalTagMetrics
        File normalMapqMetrics
        File normalInsertSizeMetrics

        Bam tumorSvBam
        File tumorCigarMetrics
        File tumorIdsvMetrics
        File tumorTagMetrics
        File tumorMapqMetrics
        File tumorInsertSizeMetrics
    }

    command {
        set -e -o pipefail

        # mv gridss fasta refs into position in the fasta dir
        fasta_dir=$( dirname ~{bwaReference.fasta} )
        mv ~{sep=" " gridssAdditionalReference} $fasta_dir

        # link assembly results
        sub_dir=~{gridssassemblyBamPath}.gridss.working/
        mkdir -p $sub_dir
        copied_sub_dir=$( dirname ~{gridssassemblySvBam.bam} )
        ln -s $copied_sub_dir/* $sub_dir

        ln -s \
        ~{gridssassemblyBam} \
        ~{gridssassemblyBamPath}
        

        # link chunk assembly results
        sub_dir=~{gridssassemblyBamPath}.gridss.working/
        mkdir -p $sub_dir
        copied_sub_dir=$( dirname ~{downsampled[0]} )
        ln -s $copied_sub_dir/* $sub_dir

        # link preprocess results
        # reposition unnamed input
        bash gridss_arrange.sh \
        --tumorFinalBam ~{tumorFinalBam.bam} \
        --normalFinalBam ~{normalFinalBam.bam} \
        --tumorSvBam ~{tumorSvBam.bam} \
        --normalSvBam ~{normalSvBam.bam}

        working=$( pwd )

        bash gridss.sh \
        --steps call \
        --reference ~{bwaReference.fasta} \
        --jar /opt/gridss/gridss-2.11.1-gridss-jar-with-dependencies.jar \
        --threads ~{threads} \
        --workingdir $working \
        --assembly ~{gridssassemblyBamPath} \
        --output ~{gridssUnfilteredVcfPath} \
        --picardoptions VALIDATION_STRINGENCY=LENIENT \
        ~{normalFinalBam.bam} ~{tumorFinalBam.bam} \
        && cat *.log
    }

    output {
        File gridssUnfilteredVcf = "~{pairName}.sv.gridss.v2.10.2.unfiltered.vcf"
    }

    runtime {
        disks: "local-disk " + diskSize + " LOCAL"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:0dd813126af95f511f0576045165a8c14fca46ea99fe44110df2f8e6b4febf72"
    }
}

task GridssFilter {
    input {
        Int threads
        Int memoryGb = 16
        Int diskSize = 4
        String bsGenome
        String pairName
        File ponTarGz
        String gridssVcfPath = "~{pairName}.sv.gridss.v2.10.2.vcf"
        String gridssVcfPathOut = "~{pairName}.sv.gridss.v2.10.2.vcf.bgz"
        String tumourordinal = 2
        File gridssUnfilteredVcf
    }

    command {
        set -e -o pipefail

        tar -zxvf ~{ponTarGz}

        working=$( pwd )

        Rscript /opt/gridss/gridss_somatic_filter.R \
        --ref ~{bsGenome} \
        --input ~{gridssUnfilteredVcf} \
        --output ~{gridssVcfPath} \
        --tumourordinal ~{tumourordinal} \
        --plotdir $working \
        --scriptdir /opt/gridss/ \
        --configdir /opt/gridss/ \
        --pondir pon/
    }

    output {
        IndexedVcf gridssVcf = object {
                vcf : "~{gridssVcfPathOut}",
                index : "~{gridssVcfPathOut}.tbi"
            }
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        cpu : threads
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/gridss@sha256:0dd813126af95f511f0576045165a8c14fca46ea99fe44110df2f8e6b4febf72"
    }
}
