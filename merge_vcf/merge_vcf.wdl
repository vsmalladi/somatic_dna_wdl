version 1.0

import "../wdl_structs.wdl"

# General tasks

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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/samtools:1.9.1"
    }
}

task IndexVcf {
    input {
        File vcfCompressed
        String vcfCompressedPath = basename(vcfCompressed)
        String vcfIndexedPath = sub(basename(vcfCompressed), "$", ".tbi")
        Int threads = 4
        Int memoryGb = 24
        Int diskSize = (ceil( size(vcfCompressed, "GB") )  * 2 ) + 20
    }

    command {
        set -e -o pipefail
        
        cp \
        ~{vcfCompressed} \
        ~{vcfCompressedPath}
        
        gatk \
        IndexFeatureFile \
        --java-options "-XX:ParallelGCThreads=4" \
        --feature-file ~{vcfCompressedPath} \
        -O ~{vcfIndexedPath}
    }

    output {
        IndexedVcf vcfCompressedIndexed = object {
                vcf : "~{vcfCompressedPath}", 
                index : "~{vcfIndexedPath}"
            }
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
    }
}

task MergePrep {
    input {
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath = sub(basename(renameMetaVcf), ".rename_metadata.vcf$", ".merge_prep.vcf")
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
    }
}

task RenameVcf {
    input {
        File prepCallerVcf
        String pairName
        String renameVcfPath = sub(basename(prepCallerVcf), ".merge_prep.vcf$", ".rename.vcf")
        String normal
        String tumor
        String tool
        Int memoryGb = 16
        Int diskSize = (ceil( size(prepCallerVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /rename_vcf.py \
        ~{prepCallerVcf} \
        ~{renameVcfPath} \
        ~{normal} \
        ~{tumor} \
        ~{tool}
    }

    output {
        File renameVcf = "~{renameVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
    }
}

task SplitMultiAllelic {
    input {
        String pairName
        String splitVcfPath
        IndexedReference referenceFa
        IndexedVcf vcfCompressedIndexed
        Int threads = 16
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
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools:1.5"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
    }
}

task RemoveContig {
    input {
        String mnvVcfPath
        String removeChromVcfPath = "~{mnvVcfPath}"
        File removeChromVcf
        Int memoryGb = 16
        Int diskSize = (ceil( size(removeChromVcf, "GB") )  * 2 ) + 4
    }

    command {
        python \
        /remove_contig.py \
        ~{removeChromVcf} \
        ~{removeChromVcfPath}
    }

    output {
        File removeContigVcf = "~{removeChromVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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

    command {        
        gatk \
        SortVcf \
        --java-options "-Xmx8g -XX:ParallelGCThreads=4" \
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
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/broad-gatk/gatk:4.1.1.0"
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
        Int threads = 16
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
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools:1.5"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bedtools:v2.26.0"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        Int threads = 16
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/lancet:v1.0.7"
    }
}

task IntersectVcfs {
    input {
        Int threads = 16
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
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools:1.5"
    }
}

task MergeColumns {
    input {
        String pairName
        String chrom
        String columnChromVcfPath = "~{pairName}.single_column.v7.~{chrom}.vcf"
        String tumor
        String normal
        File supportedChromVcf
        Int memoryGb = 16
        Int diskSize = 4
    }

    command {
        python \
        /merge_columns.py \
        ~{supportedChromVcf} \
        ~{columnChromVcfPath} \
        ~{tumor} \
        ~{normal}
    }

    output {
        File columnChromVcf = "~{columnChromVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:0.9.4"
    }
}









