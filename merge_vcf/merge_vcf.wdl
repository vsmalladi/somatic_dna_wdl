version 1.0

import "../wdl_structs.wdl"

# General tasks

task CompressVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        File vcf
        String vcfCompressedPath = sub(vcf, "$", ".gz")
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
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task IndexVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        File vcfCompressed
    }

    command {
        gatk \
        IndexFeatureFile \
        --java-options "-Xmx24576m -XX:ParallelGCThreads=4" \
        --feature-file ~{vcfCompressed}
    }

    output {
        IndexedVcf vcfCompressedIndexed = object {
                vcf : "${vcfCompressed}", 
                index : "${vcfCompressed}.tbi"
            }
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


# Prep tasks

task RenameMetadata {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        File callerVcf
        String renameMetaVcfPath = sub(callerVcf, "$", ".rename_metadata.vcf")
        String tool
        File callerVcf
    }

    command {
        python \
        rename_metadata.py \
        ~{callerVcf} \
        ~{renameMetaVcfPath} \
        ~{tool}
    }

    output {
        File renameMetaVcf = "~{renameMetaVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task MergePrepSupport {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath =  sub(renameMetaVcf, ".rename_metadata.vcf$", ".merge_prep.vcf") 
        String tool
    }

    command {
        python \
        merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool} \
        --support
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task MergePrep{
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        File renameMetaVcf
        String prepCallerVcfPath = sub(renameMetaVcf, ".rename_metadata.vcf$", ".merge_prep.vcf")
        String tool
    }

    command {
        python \
        merge_prep.py \
        --vcf ~{renameMetaVcf} \
        --out ~{prepCallerVcfPath} \
        --tool ~{tool}
    }

    output {
        File prepCallerVcf = "~{prepCallerVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task RenameVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        File prepCallerVcf
        String pairName
        String renameVcfPath = sub(prepCallerVcf, ".merge_prep.vcf$", ".rename.vcf")
        String normal
        String tumor
        String tool
    }

    command {
        python \
        rename_vcf.py \
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
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task SplitMultiAllelic {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String splitVcfPath
        IndexedReference referenceFa
        IndexedVcf vcfCompressedIndexed
    }

    command {
        bcftools \
        norm \
        -m \
        -any \
        --no-version \
        -f ~{referenceFa.fasta} \
        -o ~{splitVcfPath} \
        ~{vcfCompressedIndexed}
    }

    output {
        File splitVcf = "~{splitVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task SplitMnv {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        File splitVcf
        String mnvVcfPath = sub(splitVcf, ".split.vcf", ".split_mnvs.vcf")
        String tool
    }

    command {
        python \
        split_mnv.py \
        ~{splitVcf} \
        ~{mnvVcfPath} \
        ~{tool}
    }

    output {
        File mnvVcf = "~{mnvVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task RemoveContig {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String mnvVcfPath
        String removeChromVcfPath = "~{mnvVcfPath}"
        File removeChromVcf
    }

    command {
        python2.7 \
        remove_contig.py \
        ~{removeChromVcf} \
        ~{removeChromVcfPath}
    }

    output {
        File removeContigVcf = "~{removeChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task Gatk4MergeSortVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String sortedVcfPath
        Array[File] tempVcfs
        IndexedReference referenceFa
    }

    command {
        gatk \
        SortVcf \
        --java-options "-Xmx8196m -XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
                vcf : "~{sortedVcfPath}", 
                index : "~{sortedVcfPath}.idx"
            }
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

# Merge Callers Section

task MergeCallers {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String chrom
        String pairName
        String mergedChromVcfPath = "~{pairName}.merged_supported.v6.~{chrom}.vcf"
        Array[IndexedVcf] allVcfCompressed
        Array[File] allVcfCompressedList
    }

    command {
        bcftools \
        merge \
        -r ~{chrom} \
        --force-samples \
        --no-version \
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
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}


task StartCandidates {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String startChromVcfPath = "~{pairName}.start.merged.v6.~{chrom}.vcf"
        File knownGeneBed
        File mergedChromVcf
    }

    command {
        bedtools \
        intersect \
        -header \
        -a ~{mergedChromVcf} \
        -b ~{knownGeneBed} \
        -v \
        > ~{startChromVcfPath}
    }

    output {
        File startChromVcf = "~{startChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task GetCandidates {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String candidateChromVcfPath = "~{pairName}.candidate.merged.v6.~{chrom}.vcf"
        File startChromVcf
    }

    command {
        python \
        get_candidates.py \
        ~{startChromVcf} \
        ~{candidateChromVcfPath}
    }

    output {
        File candidateChromVcf = "~{candidateChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task VcfToBed {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String candidateChromBedPath = "~{pairName}.candidate.merged.v6.~{chrom}.bed"
        File candidateChromVcf
    }

    command {
        python \
        vcf_to_bed.py \
        ~{candidateChromVcf} \
        | bedtools \
        merge \
        > ~{candidateChromBedPath}
    }

    output {
        File candidateChromBed = "~{candidateChromBedPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task LancetConfirm {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String lancetChromVcfPath = "~{pairName}.lancet.merged.v6.~{chrom}.vcf"
        IndexedReference referenceFa
        Bam normalFinalBam
        File candidateChromBed
        Bam tumorFinalBam
    }

    command {
        mkdir ~{chrom} \
        && \
        cd ~{chrom} \
        && \
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
        File lancetChromVcf = "~{lancetChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task IntersectVcfs {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String vcfConfirmedCandidatePath = "~{pairName}.confirmed_lancet.merged.v6.~{chrom}.vcf"
        IndexedVcf vcfCompressedLancet
        IndexedVcf vcfCompressedCandidate
    }

    command {
        bcftools \
        isec \
        -w 1 \
        -c none \
        -n =2 \
        ~{vcfCompressedLancet.vcf} \
        ~{vcfCompressedCandidate.vcf} \
        > ~{vcfConfirmedCandidatePath}
    }

    output {
        File vcfConfirmedCandidate = "~{vcfConfirmedCandidatePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task MergeColumns {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String columnChromVcfPath = "~{pairName}.single_column.v6.~{chrom}.vcf"
        String tumor
        String normal
        File supportedChromVcf
    }

    command {
        python \
        merge_columns.py \
        ~{supportedChromVcf} \
        ~{columnChromVcfPath} \
        ~{tumor} \
        ~{normal}
    }

    output {
        File columnChromVcf = "~{columnChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task AddNygcAlleleCountsToVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String preCountsChromVcfPath = "~{pairName}.pre_count.v6.~{chrom}.vcf"
        Bam normalFinalBam
        Bam tumorFinalBam
        File columnChromVcf
    }

    command {
        python \
        add_nygc_allele_counts_to_vcf.py \
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
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task AddFinalAlleleCountsToVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String countsChromVcfPath = "~{pairName}.final.v6.~{chrom}.vcf"
        File preCountsChromVcf
    }

    command {
        python \
        add_final_allele_counts_to_vcf.py \
        -v ~{preCountsChromVcf} \
        -o ~{countsChromVcfPath} \
    }

    output {
        File countsChromVcf = "~{countsChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task FilterPon {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String chrom
        String pairName
        String ponOutFilePath = "~{pairName}.pon.final.v6.~{chrom}.vcf"
        File countsChromVcf
        File ponFile
    }

    command {
        python \
        filter_pon.py \
        --bed ~{ponFile} \
        --chrom ~{chrom} \
        --vcf ~{countsChromVcf} \
        --out ~{ponOutFilePath} \
    }

    output {
        File ponOutFile = "~{ponOutFilePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task FilterVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String filteredOutFilePath = "~{pairName}.final.v6.filtered.~{chrom}.vcf"
        IndexedVcf germFile
        File ponOutFile
    }

    command {
        python \
        filter_vcf.py \
        ~{germFile.vcf} \
        ~{ponOutFile} \
        ~{filteredOutFilePath}
    }

    output {
        File filteredOutFile = "~{filteredOutFilePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task SnvstomnvsCountsbasedfilterAnnotatehighconf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String chrom
        String finalChromVcfPath = "~{pairName}.mnv.final.v6.filtered.~{chrom}.vcf"
        File filteredOutFile
    }

    command {
        python2.7 \
        SNVsToMNVs_CountsBasedFilter_AnnotateHighConf.py \
        -i ~{filteredOutFile} \
        -o ~{finalChromVcfPath} \
    }

    output {
        File finalChromVcf = "~{finalChromVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}









