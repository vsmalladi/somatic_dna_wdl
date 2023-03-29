version 1.0

import "../wdl_structs.wdl"

# General tasks

task getIntersection {
    input {
        String pairName
        String intersectedVcfPath = "~{pairName}.wgs_in_exonic.vcf"
        File intervalListBed
        File wgsVcf
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bedtools \
        intersect \
        -header \
        -a ~{wgsVcf} \
        -b ~{intervalListBed} \
        -F 1 \
        > ~{intersectedVcfPath}
    }

    output {
        File intersectedVcf = "~{intersectedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bedtools@sha256:9e737f5c96c00cf3b813d419d7a7b474c4013c9aa9dfe704eb36417570c6474e"
    }
}

task RenameExomeWgsMetadata {
    input {
        String pairName
        File callerVcf
        String renameMetaVcfPath = sub(basename(callerVcf), "$", ".rename_metadata.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
        # public copy avilable in
        # https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/tools/rename_center_metadata.py?at=refs%2Fheads%2Fexome_wgs_merge
        File renameExomeWgsMetadata = "/gpfs/commons/projects/TCGA/gdc-awg/MILD/WGS/scripts/exome/somatic_dna_wdl/tools/rename_exome_wgs_metadata.py"
    }

    command {
        python \
        ~{renameExomeWgsMetadata} \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    }

    output {
        File renameMetaVcf = "~{renameMetaVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task MergePrepExomeWgs {
    input {
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath = sub(basename(renameMetaVcf, ".gz"), ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(renameMetaVcf, "GB") )  * 2 ) + 4
        # public copy avilable in
        # https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_wdl/browse/tools/merge_exome_wgs_prep.py?at=refs%2Fheads%2Fexome_wgs_merge
        File mergeExomeWgsPrep = "/gpfs/commons/projects/TCGA/gdc-awg/MILD/WGS/scripts/exome/somatic_dna_wdl/tools/merge_exome_wgs_prep.py"
    }

    command {
        python \
        ~{mergeExomeWgsPrep} \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task MergeExomeWgsCallers {
    input {
        String chrom
        String pairName
        String mergedChromVcfPath = "~{pairName}.merged_supported.v7.~{chrom}.vcf"
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bcftools \
        merge \
        -r ~{chrom} \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -f PASS \
        -F x \
        -m none \
        -o ~{mergedChromVcfPath} \
        -i called_by:join,num_callers:sum,MNV_ID:join \
        ~{sep=" " allVcfCompressedList}
    }

    output {
        File mergedChromVcf = "~{mergedChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
    }
}

task SnvstomnvsAnnotateExomeWgsCalled {
    input {
        String pairName
        String chrom
        String finalChromVcfPath = "~{pairName}.mnv.final.v7.filtered.~{chrom}.vcf"
        File filteredOutFile
        Int memoryGb = 16
        Int diskSize = 4
        
        File snvsToMnvs = "/gpfs/commons/projects/TCGA/gdc-awg/MILD/WGS/scripts/exome/somatic_dna_wdl/tools/SNVsToMNVs_AnnotateExomeWgsCalled.py"
        File classes = "/gpfs/commons/projects/TCGA/gdc-awg/MILD/WGS/scripts/exome/somatic_dna_tools/Classes.py"
    }

    command {
        set -e -o pipefail
        
        cp \
        ~{snvsToMnvs} \
        SNVsToMNVs_AnnotateExomeWgsCalled.py
        
        cp \
        ~{classes} \
        Classes.py
        
        python \
        SNVsToMNVs_AnnotateExomeWgsCalled.py \
        -i ~{filteredOutFile} \
        -o ~{finalChromVcfPath} \
    }

    output {
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}


task CompressVcf {
    input {
        File vcf
        String vcfCompressedPath = sub(basename(vcf), "$", ".gz")
        Int memoryGb = 4
        Int diskSize = (ceil( size(vcf, "GB") )  * 2 ) + 4
    }

    command {
        bgzip \
        -c \
        ~{vcf} \
        > ~{vcfCompressedPath} \
    }

    output {
        File vcfCompressed = "~{vcfCompressedPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/samtools@sha256:963b0b2f24908832efab8ddccb7a7f3ba5dca9803bc099be7cf3a455766610fd"
    }
}

task IndexVcf {
    input {
        File vcfCompressed
        String vcfCompressedPath = basename(vcfCompressed)
        String vcfIndexedPath = sub(basename(vcfCompressed), "$", ".tbi")
        Int threads = 4
        Int memoryGb = 8
        Int diskSize = (ceil( size(vcfCompressed, "GB") )  * 2 ) + 4
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
        set -e -o pipefail

        cp \
        ~{vcfCompressed} \
        ~{vcfCompressedPath}

        gatk \
        IndexFeatureFile \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        --input ~{vcfCompressedPath} \
        -O ~{vcfIndexedPath}
    }

    output {
        IndexedVcf vcfCompressedIndexed = object {
                vcf : "~{vcfCompressedPath}",
                index : "~{vcfIndexedPath}"
            }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}


# Prep tasks

task RenameMetadata {
    input {
        String pairName
        File callerVcf
        String renameMetaVcfPath = sub(basename(callerVcf), "$", ".rename_metadata.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(callerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /rename_metadata.py \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    }

    output {
        File renameMetaVcf = "~{renameMetaVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task MergePrepSupport {
    input {
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath =  sub(basename(renameMetaVcf), ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(renameMetaVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool} \
        --support
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task MergePrep {
    input {
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath = sub(basename(renameMetaVcf, ".gz"), ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(renameMetaVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task RenameVcf {
    input {
        File prepCallerVcf
        String pairName
        String renameVcfPath = sub(basename(prepCallerVcf, ".gz"), ".merge_prep.vcf$", ".rename.vcf")
        String normalId
        String tumorId
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(prepCallerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /rename_vcf.py \
        ~{prepCallerVcf} \
        ~{renameVcfPath} \
        ~{normalId} \
        ~{tumorId} \
        ~{tool}
    }

    output {
        File renameVcf = "~{renameVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task RenameVcfPon {
    input {
        File prepCallerVcf
        String renameVcfPath = sub(basename(prepCallerVcf, ".gz"), ".vcf$", ".rename.vcf")
        String tumor
        String tool
        File renameVcfPon = "gs://nygc-comp-s-fd4e-input/scripts/rename_vcf_pon.py"
        Int memoryGb = 16
        Int diskSize = (ceil( size(prepCallerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        ~{renameVcfPon} \
        ~{prepCallerVcf} \
        ~{renameVcfPath} \
        ~{tumor} \
        ~{tool}
    }

    output {
        File renameVcf = "~{renameVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task SplitMultiAllelic {
    input {
        String pairName
        String splitVcfPath
        IndexedReference referenceFa
        IndexedVcf vcfCompressedIndexed
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = (ceil( size(vcfCompressedIndexed.vcf, "GB") )  * 2 ) + 10
    }

    command {
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{vcfCompressedIndexed.vcf}
    }

    output {
        File splitVcf = "~{splitVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
    }
}

task SplitMultiAllelicRegions {
    input {
        String pairName
        String splitVcfPath
        IndexedReference referenceFa
        IndexedVcf vcfCompressedIndexed
        Array[String]+ listOfChroms
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = (ceil( size(vcfCompressedIndexed.vcf, "GB") )  * 3 ) + 10
    }

    command {
        bcftools \
        norm \
        -m \
        -any \
        --threads ~{threads} \
        --regions ~{sep="," listOfChroms} \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{vcfCompressedIndexed.vcf}
    }

    output {
        File sortedVcf = "~{splitVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
    }
}

task SplitMnv {
    input {
        File splitVcf
        String mnvVcfPath = sub(basename(splitVcf), ".split.vcf", ".split_mnvs.vcf")
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(splitVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /split_mnv.py \
        ~{splitVcf} \
        ~{mnvVcfPath} \
        ~{tool}
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task RemoveContig {
    input {
        String mnvVcfPath
        String removeChromVcfPath = "~{mnvVcfPath}"
        File removeChromVcf
        String dollarSign = "$"
        Int memoryGb = 16
        Int diskSize = (ceil( size(removeChromVcf, "GB") )  * 2 ) + 4
    }

    command <<<
        set -e -o pipefail

        filename=$(basename -- "~{dollarSign}~{removeChromVcf}")
        extension="~{dollarSign}{filename##*.}"
        filename="~{dollarSign}{filename%.*}"

        if [[ ~{dollarSign}extension == gz ]]; then
            input_path=~{dollarSign}filename

            gunzip -c \
            ~{removeChromVcf} \
            > ~{dollarSign}input_path
        else
            input_path=~{removeChromVcf}
        fi

        python \
        /remove_contig.py \
        ~{dollarSign}input_path \
        ~{removeChromVcfPath}
    >>>

    output {
        File removeContigVcf = "~{removeChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task Gatk4MergeSortVcf {
    input {
        String sortedVcfPath
        Array[File] tempVcfs
        IndexedReference referenceFa
        Boolean gzipped = false
        String suffix = if gzipped then ".tbi" else ".idx"
        Int threads = 4
        Int memoryGb = 8
        Int diskSize = 10
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        gatk \
        SortVcf \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
            vcf : "~{sortedVcfPath}",
            index : "~{sortedVcfPath}~{suffix}"
        }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

# Merge Callers Section

task MergeCallers {
    input {
        String chrom
        String pairName
        String mergedChromVcfPath = "~{pairName}.merged_supported.v7.~{chrom}.vcf"
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bcftools \
        merge \
        -r ~{chrom} \
        --force-samples \
        --no-version \
        --threads ~{threads} \
        -f PASS,SUPPORT \
        -F x \
        -m none \
        -o ~{mergedChromVcfPath} \
        -i called_by:join,num_callers:sum,MNV_ID:join,supported_by:join \
        ~{sep=" " allVcfCompressedList}
    }

    output {
        File mergedChromVcf = "~{mergedChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
    }
}


task StartCandidates {
    input {
        String pairName
        String chrom
        String startChromVcfPath = "~{pairName}.start.merged.v7.~{chrom}.vcf"
        File intervalListBed
        File mergedChromVcf
        Int memoryGb = 16
        Int diskSize = 20
    }

    command {
        bedtools \
        intersect \
        -header \
        -a ~{mergedChromVcf} \
        -b ~{intervalListBed} \
        -v \
        > ~{startChromVcfPath}
    }

    output {
        File startChromVcf = "~{startChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bedtools@sha256:9e737f5c96c00cf3b813d419d7a7b474c4013c9aa9dfe704eb36417570c6474e"
    }
}

task GetCandidates {
    input {
        String pairName
        String chrom
        String candidateChromVcfPath = "~{pairName}.candidate.merged.v7.~{chrom}.vcf"
        File startChromVcf
        Int memoryGb = 16
        Int diskSize = (ceil( size(startChromVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /get_candidates.py \
        ~{startChromVcf} \
        ~{candidateChromVcfPath}
    }

    output {
        File candidateChromVcf = "~{candidateChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task VcfToBed {
    input {
        String pairName
        String chrom
        String candidateChromBedPath = "~{pairName}.candidate.merged.v7.~{chrom}.bed"
        File candidateChromVcf
        Int memoryGb = 16
        Int diskSize = (ceil( size(candidateChromVcf, "GB") )  * 2 ) + 4
    }

    command {
        set -e -o pipefail

        python \
        /vcf_to_bed.py \
        ~{candidateChromVcf} \
        | bedtools \
        merge \
        > ~{candidateChromBedPath}
    }

    output {
        File candidateChromBed = "~{candidateChromBedPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task LancetConfirm {
    input {
        String pairName
        String chrom
        String lancetChromVcfPath = "~{pairName}.lancet.merged.v7.~{chrom}.vcf"
        IndexedReference referenceFa
        Bam normalFinalBam
        File candidateChromBed
        Bam tumorFinalBam
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = (ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) ) + 10
    }

    command {
        set -e -o pipefail

        mkdir ~{chrom}

        cd ~{chrom}

        lancet \
        --normal ~{normalFinalBam.bam} \
        --tumor ~{tumorFinalBam.bam} \
        --bed ~{candidateChromBed} \
        --ref ~{referenceFa.fasta} \
        --min-k 11 \
        --low-cov 1 \
        --min-phred-fisher 5 \
        --min-strand-bias 1 \
        --min-alt-count-tumor 3 \
        --min-vaf-tumor 0.04 \
        --padding 250 \
        --window-size 2000 \
        --num-threads ~{threads} \
        > ~{lancetChromVcfPath}
    }

    output {
        File lancetChromVcf = "~{chrom}/~{lancetChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet@sha256:25169d34b41de9564e03f02ebcbfb4655cf536449592b0bd58773195f9376e61"
    }
}

task IntersectVcfs {
    input {
        Int threads = 8
        Int memoryGb = 16
        Int diskSize = 4
        String pairName
        String chrom
        String vcfConfirmedCandidatePath = "~{pairName}.confirmed_lancet.merged.v7.~{chrom}.vcf"
        IndexedVcf vcfCompressedLancet
        IndexedVcf vcfCompressedCandidate
    }

    command {
        bcftools \
        isec \
        -w 1 \
        -c none \
        -n =2 \
        --threads ~{threads} \
        ~{vcfCompressedLancet.vcf} \
        ~{vcfCompressedCandidate.vcf} \
        > ~{vcfConfirmedCandidatePath}
    }

    output {
        File vcfConfirmedCandidate = "~{vcfConfirmedCandidatePath}"
    }

    runtime {
        mem: memoryGb + "G"
        cpus: threads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools@sha256:d9b84254b8cc29fcae76b728a5a9a9a0ef0662ee52893d8d82446142876fb400"
    }
}

task MergeColumns {
    input {
        String pairName
        String chrom
        String columnChromVcfPath = "~{pairName}.single_column.v7.~{chrom}.vcf"
        String tumorId
        String normalId
        File supportedChromVcf
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /merge_columns.py \
        ~{supportedChromVcf} \
        ~{columnChromVcfPath} \
        ~{tumorId} \
        ~{normalId}
    }

    output {
        File columnChromVcf = "~{columnChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task MergeColumnsPon {
    input {
        String chrom
        String tumor
        String columnChromVcfPath = "~{tumor}.single_column.v7.~{chrom}.vcf"
        File supportedChromVcf
        Int memoryGb = 16
        Int diskSize = 4
        File MergeColumnsPon = "gs://nygc-comp-s-fd4e-input/scripts/merge_columns_pon.py"
    }

    command {
        python \
        ~{MergeColumnsPon} \
        ~{supportedChromVcf} \
        ~{columnChromVcfPath} \
        ~{tumor}
    }

    output {
        File columnChromVcf = "~{columnChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task AddNygcAlleleCountsToVcf {
    input {
        String pairName
        String chrom
        String preCountsChromVcfPath = "~{pairName}.pre_count.v7.~{chrom}.vcf"
        Bam normalFinalBam
        Bam tumorFinalBam
        File columnChromVcf
        Int memoryGb = 40
        Int diskSize = ceil( size(tumorFinalBam.bam, "GB") + size(normalFinalBam.bam, "GB")) + 20
    }

    command {
        python \
        /add_nygc_allele_counts_to_vcf.py \
        -t ~{tumorFinalBam.bam} \
        -n ~{normalFinalBam.bam} \
        -v ~{columnChromVcf} \
        -b 10 \
        -m 10 \
        -o ~{preCountsChromVcfPath}
    }

    output {
        File preCountsChromVcf = "~{preCountsChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task AddFinalAlleleCountsToVcf {
    input {
        String pairName
        String chrom
        String countsChromVcfPath = "~{pairName}.final.v7.~{chrom}.vcf"
        File preCountsChromVcf
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /add_final_allele_counts_to_vcf.py \
        -v ~{preCountsChromVcf} \
        -o ~{countsChromVcfPath} \
    }

    output {
        File countsChromVcf = "~{countsChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task FilterPon {
    input {
        String chrom
        String pairName
        String ponOutFilePath = "~{pairName}.pon.final.v7.~{chrom}.vcf"
        File countsChromVcf
        File ponFile
        Int memoryGb = 60
        Int diskSize = 16
    }

    command {
        python \
        /filter_pon.py \
        --bed ~{ponFile} \
        --chrom ~{chrom} \
        --vcf ~{countsChromVcf} \
        --out ~{ponOutFilePath} \
    }

    output {
        File ponOutFile = "~{ponOutFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task FilterVcf {
    input {
        String pairName
        String chrom
        String filteredOutFilePath = "~{pairName}.final.v7.filtered.~{chrom}.vcf"
        IndexedVcf germFile
        File ponOutFile
        Int memoryGb = 60
        Int diskSize = 16
    }

    command {
        python \
        /filter_vcf.py \
        ~{germFile.vcf} \
        ~{ponOutFile} \
        ~{filteredOutFilePath}
    }

    output {
        File filteredOutFile = "~{filteredOutFilePath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}

task SnvstomnvsCountsbasedfilterAnnotatehighconf {
    input {
        String pairName
        String chrom
        String finalChromVcfPath = "~{pairName}.mnv.final.v7.filtered.~{chrom}.vcf"
        File filteredOutFile
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /SNVsToMNVs_CountsBasedFilter_AnnotateHighConf.py \
        -i ~{filteredOutFile} \
        -o ~{finalChromVcfPath} \
    }

    output {
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:ffc8cfe7dd42baf87b56f121e2f3ad72b6897a31c56328c989d56379c9fbc7ae"
    }
}
