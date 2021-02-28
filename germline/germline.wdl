version 1.0

import "../wdl_structs.wdl"

task Haplotypecaller {
    input {
        Bam finalBam
        IndexedReference referenceFa
        String chrom
        String sampleId
        String haplotypecallerChromVcfPath = "~{sampleId}.~{chrom}.haplotypeCalls.er.raw.vcf"
        Int threads = 18
        Int memoryGb = 24
        Int diskSize
    }

    command {
        java \
        -Xmx 24g \
        -XX:ParallelGCThreads=2 \
        -jar /usr/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        --genotyping_mode DISCOVERY \
        -A AlleleBalanceBySample \
        -A DepthPerAlleleBySample \
        -A DepthPerSampleHC \
        -A InbreedingCoeff \
        -A MappingQualityZeroBySample \
        -A StrandBiasBySample \
        -A Coverage \
        -A FisherStrand \
        -A HaplotypeScore \
        -A MappingQualityRankSumTest \
        -A MappingQualityZero \
        -A QualByDepth \
        -A RMSMappingQuality \
        -A ReadPosRankSumTest \
        -A VariantType \
        -l INFO \
        --emitRefConfidence GVCF \
        -rf BadCigar \
        --variant_index_parameter 128000 \
        --variant_index_type LINEAR \
        -R ~{referenceFa.fasta} \
        -nct 16 \
        -I ~{finalBam.bam} \
        -L ~{chrom} \
        -o ~{haplotypecallerChromVcfPath}
    }

    output {
        File haplotypecallerChromVcf = "~{haplotypecallerChromVcfPath} "
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk3:3.5-0"
    }
}

task GentotypeGvcfs {
    input {
        IndexedReference referenceFa
        IndexedVcf sortedVcf
        String sampleId
        String haplotypecallerGenoVcfPath = "~{sampleId}.genotypeGVCFs.vcf"
        Int threads = 12
        Int memoryGb = 24
        Int diskSize = (ceil( size(sortedVcf.vcf, "GB") )  * 2 ) + 4
    }

    command {
        java \
        -Xmx 24g \
        -XX:ParallelGCThreads=4 \
        -jar /usr/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R ~{referenceFa.fasta} \
        -nt 8 \
        --disable_auto_index_creation_and_locking_when_reading_rods \
        --variant ~{sortedVcf.vcf} \
        -o ~{haplotypecallerGenoVcfPath}
    }

    output {
        File haplotypecallerGenoVcf = "~{haplotypecallerGenoVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk3:3.5-0"
    }
}

task RecalVcfsSnp {
    input {
        String sampleId
        String recalSnpPath = "~{sampleId}.recalibrate_SNP.recal"
        String tranchesSnpPath = "~{sampleId}.recalibrate_SNP.tranches"
        String recalibratePlotsPath = "~{sampleId}.recalibrate_SNP_plots.R"
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedReference referenceFa
        IndexedVcf onekG
        IndexedVcf dbsnp
        File haplotypecallerGenoVcf
        Int memoryGb = 32
        Int threads = 8
        Int diskSize = (ceil( size(haplotypecallerGenoVcf, "GB") )  * 2 ) + 20
    }

    command {
        java \
        -Xmx 32g \
        -XX:ParallelGCThreads=4 \
        -jar /usr/GenomeAnalysisTK.jar \
        -T VariantRecalibrator \
        --maxGaussians 4 \
        -tranche 100.0 \
        -tranche 99.8 \
        -tranche 99.6 \
        -tranche 99.4 \
        -tranche 99.2 \
        -tranche 99.0 \
        -tranche 95.0 \
        -tranche 90.0 \
        -an QD \
        -an MQ \
        -an FS \
        -an MQRankSum \
        -an ReadPosRankSum \
        -an SOR \
        -an DP \
        -mode SNP \
        -R ~{referenceFa.fasta} \
        -nt 4 \
        --input ~{haplotypecallerGenoVcf} \
        -recalFile ~{recalSnpPath} \
        -tranchesFile ~{tranchesSnpPath} \
        -rscriptFile ~{recalibratePlotsPath} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ~{dbsnp.vcf} \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ~{hapmap.vcf} \
        -resource:omni,known=false,training=true,truth=true,prior=12.0 ~{omni.vcf} \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 ~{onekG.vcf}
    }

    output {
        File tranchesSnp = "~{tranchesSnpPath}"
        File recalibratePlots = "~{recalibratePlotsPath}"
        File recalSnp = "~{recalSnpPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk3:3.5-0"
    }
}

task RecalVcfsIndel {
    input {
        String sampleId
        String recalIndelPath = "~{sampleId}.recalibrate_INDEL.recal"
        String tranchesIndelPath = "~{sampleId}.recalibrate_INDEL.tranches"
        String recalibratePlotsIndelPath = "~{sampleId}.recalibrate_INDEL_plots.R"
        IndexedReference referenceFa
        IndexedVcf dbsnp
        File haplotypecallerGenoVcf
        IndexedVcf MillsAnd1000G
        Int memoryGb = 32
        Int threads = 8
        Int diskSize = (ceil( size(haplotypecallerGenoVcf, "GB") )  * 2 ) + 20
    }

    command {
        java \
        -Xmx 32g \
        -XX:ParallelGCThreads=4 \
        -jar /usr/GenomeAnalysisTK.jar \
        -T VariantRecalibrator \
        --maxGaussians 4 \
        -tranche 100.0 \
        -tranche 99.0 \
        -tranche 95.0 \
        -tranche 92.0 \
        -tranche 90.0 \
        -an QD \
        -an FS \
        -an MQRankSum \
        -an ReadPosRankSum \
        -an SOR \
        -an DP \
        -mode INDEL \
        -R ~{referenceFa.fasta} \
        -nt 4 \
        --input ~{haplotypecallerGenoVcf} \
        -recalFile ~{recalIndelPath} \
        -tranchesFile ~{tranchesIndelPath} \
        -rscriptFile ~{recalibratePlotsIndelPath} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0  ~{dbsnp.vcf} \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 ~{MillsAnd1000G.vcf}
    }
    
    output {
        File recalibratePlotsIndel = "~{recalibratePlotsIndelPath}"
        File tranchesIndel = "~{tranchesIndelPath}"
        File recalIndel = "~{recalIndelPath}"
    }
    

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk3:3.5-0"
    }
}

task ApplyRecal {
    input {
        String mode
        String sampleId
        String haplotypecallerGenoVcfApplySnpIndelPath
        IndexedReference referenceFa
        File tranches
        File recal
        File haplotypecallerGenoVcfApply
        Int threads = 20
        Int memoryGb = 8
        Int diskSize = ceil( size(haplotypecallerGenoVcfApply, "GB") + size(recal, "GB") + size(tranches, "GB") ) + 20
    }

    command {
        java \
        -Xmx 8g \
        -XX:ParallelGCThreads=4 \
        -jar /usr/GenomeAnalysisTK.jar \
        -T ApplyRecalibration \
        --ts_filter_level 99.6 \
        -R ~{referenceFa.fasta} \
        -nt 16 \
        --input ~{haplotypecallerGenoVcfApply} \
        -mode ~{mode} \
        -recalFile ~{recal} \
        -tranchesFile ~{tranches} \
        -o ~{haplotypecallerGenoVcfApplySnpIndelPath} \
    }

    output {
        File haplotypecallerGenoVcfApplySnpIndel = "~{haplotypecallerGenoVcfApplySnpIndelPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk3:3.5-0"
    }
}

task VarFilter {
    input {
        String sampleId
        String haplotypecallerRecalVcfPath = "~{sampleId}.recalibrated.haplotypeCalls.vcf"
        IndexedReference referenceFa
        File haplotypecallerGenoVcfApplySnpIndel
        Int threads = 2
        Int memoryGb = 24
        Int diskSize = (ceil( size(haplotypecallerGenoVcfApplySnpIndel, "GB") ) * 2 ) + 10
    }

    command {
        set -o pipefail && \
        java \
        -Xmx 24g \
        -XX:ParallelGCThreads=2 \
        -jar /usr/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R ~{referenceFa.fasta} \
        -V ~{haplotypecallerGenoVcfApplySnpIndel} \
        --genotypeFilterName "GQbelow20" \
        --genotypeFilterExpression "GQ < 20.0" \
        --genotypeFilterName "DPbelow10" \
        --genotypeFilterExpression "DP < 10.0" \
        | varfilt_reheader \
        > ~{haplotypecallerRecalVcfPath} \
    }

    output {
        File haplotypecallerRecalVcf = "~{haplotypecallerRecalVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk3:3.5-0"
    }
}

task VcfNorm {
    input {
        String sampleId
        String haplotypecallerNormVcfPath = "~{sampleId}.recalibrated.haplotypeCalls.norm.vcf.gz"
        IndexedVcf haplotypecallerRecalVcf
        Int threads = 4
        Int memoryGb = 20
        Int diskSize = (ceil( size(haplotypecallerRecalVcf.vcf, "GB") ) * 2 ) + 5
    }

    command {
        bcftools \
        norm \
        --threads ~{threads} \
        -O z \
        -m-both ~{haplotypecallerRecalVcf.vcf} \
        -o ~{haplotypecallerNormVcfPath} \
    }

    output {
        File haplotypecallerNormVcf = "~{haplotypecallerNormVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/bcftools:1.5"
    }
}

task VarEval {
    input {
        String sampleId
        String varReportPath = "~{sampleId}.VariantEval.report.jg.txt"
        IndexedReference referenceFa
        IndexedVcf dbsnp
        File haplotypecallerGenoVcfApplySnpIndel
        Int threads = 4
        Int memoryGb = 24
        Int diskSize = (ceil( size(haplotypecallerGenoVcfApplySnpIndel, "GB") ) * 2 ) + 20    
    }

    command {
        java \
        -XX:ParallelGCThreads=4 \
        -Xmx 24g \
        -jar /usr/GenomeAnalysisTK.jar \
        -T VariantEval \
        -ST Sample \
        -ST Novelty \
        -EV CompOverlap \
        -EV CountVariants \
        -EV IndelSummary \
        -EV TiTvVariantEvaluator \
        -noEV \
        -noST \
        -R ~{referenceFa.fasta} \
        --eval ~{haplotypecallerGenoVcfApplySnpIndel} \
        -o ~{varReportPath} \
        --dbsnp ~{dbsnp.vcf} \
    }

    output {
        File varReport = "~{varReportPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk3:3.5-0"
    }
}

task VarSum {
    input {
        String sampleId
        String varSummaryPath = "~{sampleId}.VariantEval.summary.jg.txt"
        File varReport
        Int memoryGb = 4
        Int diskSize = ceil( size(varReport, "GB") ) + 4
    }

    command {
        python2.7 \
        /make_vareval_summary.py \
        ~{varReport} \
        ~{varSummaryPath} \
    }

    output {
        File varSummary = "~{sampleId}.VariantEval.summary.jg.txt"
    }

    runtime {
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
    }
}

