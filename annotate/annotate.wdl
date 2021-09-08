version 1.0

import "../wdl_structs.wdl"

# General tasks


task RemoveSpanning {
    input {
        String sampleId
        String noSpanningVcfPath = "~{sampleId}.v7.remove_spanning.vep.annotated.vcf"
        File vcfAnnotatedVep
        Int threads = 4
        Int memoryGb = 16
        Int diskSize = ceil( size(vcfAnnotatedVep, "GB") * 2) + 10
    }

    command {
        bcftools \
        view \
        --exclude 'ALT="*"' \
        --threads ~{threads} \
        -o ~{noSpanningVcfPath} \
        ~{vcfAnnotatedVep}
    }

    output {
        File noSpanningVcf = "~{noSpanningVcfPath}"
    }

    runtime {
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/bcftools:1.5"
    }
}

task AddCosmic {
    input {
        String pairName
        String vcfAnnotatedCancerGeneCensusPath = "~{pairName}.v7.cosmic_census.vep.annotated.vcf"
        File vcfAnnotatedVep
        File cosmicCensus
        Int memoryGb = 2
        Int diskSize = ceil( size(vcfAnnotatedVep, "GB") * 2) + 5
    }

    command {
        python \
        /add_cancer_gene_census.py \
        ~{cosmicCensus} \
        ~{vcfAnnotatedVep} \
        ~{vcfAnnotatedCancerGeneCensusPath}
    }

    output {
        File vcfAnnotatedCancerGeneCensus = "~{vcfAnnotatedCancerGeneCensusPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task AddCancerResistanceMutations {
    input {
        String pairName
        String vcfAnnotatedResistancePath = "~{pairName}.v7.resistance.vep.annotated.vcf"
        File vcfAnnotatedCancerGeneCensus
        File cancerResistanceMutations
        Int memoryGb = 2
        Int diskSize = ceil( size(vcfAnnotatedCancerGeneCensus, "GB") * 2) + ceil( size(cancerResistanceMutations, "GB")) + 5
    }

    command {
        python \
        /add_cancer_resistance_mutations.py \
        ~{cancerResistanceMutations} \
        ~{vcfAnnotatedCancerGeneCensus} \
        ~{vcfAnnotatedResistancePath}
    }

    output {
        File vcfAnnotatedResistance = "~{vcfAnnotatedResistancePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task AddCancerResistanceMutationsFinal {
    input {
        String pairName
        String vcfAnnotatedResistancePath = "~{pairName}.v7.resistance.vep.annotated.vcf"
        File vcfAnnotatedCancerGeneCensus
        File cancerResistanceMutations
        Int memoryGb = 2
        Int diskSize = ceil( size(vcfAnnotatedCancerGeneCensus, "GB") * 2) + ceil( size(cancerResistanceMutations, "GB")) + 5
    }

    command {
        python \
        /add_cancer_resistance_mutations.py \
        ~{cancerResistanceMutations} \
        ~{vcfAnnotatedCancerGeneCensus} \
        ~{vcfAnnotatedResistancePath}
    }

    output {
        File vcfAnnotatedResistance = "~{vcfAnnotatedResistancePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task AnnotateId {
    input {
        Int memoryGb = 2
        String pairName
        String vcfAnnotatedIdPath = "~{pairName}.v7.id.vep.annotated.vcf"
        File vcfAnnotatedResistance
        Int diskSize = ceil( size(vcfAnnotatedResistance, "GB") * 2) + 5
    }

    command {
        python \
        /annotate_id.py \
        ~{vcfAnnotatedResistance} \
        ~{vcfAnnotatedIdPath}
    }

    output {
        File vcfAnnotatedId = "~{vcfAnnotatedIdPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task RenameCsqVcf {
    input {
        Int memoryGb = 2
        String pairName
        String vcfCsqRenamedPath = "~{pairName}.snv.indel.supplemental.v7.annotated.vcf"
        File vcfAnnotatedId
        Int diskSize = ceil( size(vcfAnnotatedId, "GB") * 2) + 5
    }

    command {
        python \
        /rename_csq_vcf.py \
        ~{vcfAnnotatedId} \
        ~{vcfCsqRenamedPath}
    }

    output {
        File vcfCsqRenamed = "~{vcfCsqRenamedPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task MainVcf {
    input {
        Int memoryGb = 2
        String pairName
        String mainVcfPath = "~{pairName}.snv.indel.final.v7.annotated.vcf"
        File vcfAnnotated
        Int diskSize = ceil( size(vcfAnnotated, "GB") * 2) + 5
    }

    command {
        python \
        /make_main_vcf.py \
        ~{vcfAnnotated} \
        ~{mainVcfPath}
    }

    output {
        File mainVcf = "~{mainVcfPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task TableVcf {
    input {
        Int memoryGb = 2
        String pairName
        String vcfAnnotatedTxtPath = "~{pairName}.snv.indel.final.v7.annotated.txt"
        String tumor
        String normal
        File mainVcf
        Int diskSize = ceil( size(mainVcf, "GB") * 2) + 5
    }

    command {
        /make_txt.py \
        --vcf ~{mainVcf} \
        --txt ~{vcfAnnotatedTxtPath} \
        --tumor ~{tumor} \
        --normal ~{normal}
    }

    output {
        File vcfAnnotatedTxt = "~{vcfAnnotatedTxtPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task VcfToMaf {
    input {
        Int memoryGb = 2
        String pairName
        String mafPath = "~{pairName}.snv.indel.final.v7.annotated.maf"
        String library
        String vepGenomeBuild
        String tumor
        String normal
        File ensemblEntrez
        File mainVcf
        Int diskSize = ceil( size(mainVcf, "GB") * 2) + 5
    }

    command {
        python \
        /make_maf.py \
        --vcf ~{mainVcf} \
        --maf ~{mafPath} \
        --library ~{library} \
        --vep-version ~{vepGenomeBuild} \
        --tumor ~{tumor} \
        --normal ~{normal} \
        --ensembl-entrez ~{ensemblEntrez}
    }

    output {
        File maf = "~{mafPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/somatic_tools:v1.1.2"
    }
}

task annotateBicSeq2Cnv {
    input {
        Int memoryGb = 24
        String pairName
        Array[String] listOfChroms
        String tumor
        String normal
        File bicseq2
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed
        
        String cnvAnnotatedFinalBedPath  = "~{pairName}.cnv.annotated.v7.final.bed"
        String cnvAnnotatedSupplementalBedPath  = "~{pairName}.cnv.annotated.v7.supplemental.bed"
        Int diskSize = 4
    }

    command {
        Rscript \
        /annotate-cnv.r \
        --cnv=~{bicseq2} \
        --caller=bicseq2 \
        --tumor=~{tumor} \
        --normal=~{normal} \
        --cytoband=~{cytoBand} \
        --db_names=DGV,1000G,COSMIC \
        --db_files=~{dgv},~{thousandG},~{cosmicUniqueBed} \
        --cancer_census=~{cancerCensusBed} \
        --ensembl=~{ensemblUniqueBed} \
        --allowed_chr=~{sep="," listOfChroms} \
        --overlap_fraction=0.8 \
        --out_file_main=~{cnvAnnotatedFinalBedPath} \
        --out_file_supplemental=~{cnvAnnotatedSupplementalBedPath}
    }

    output {
        File cnvAnnotatedFinalBed = "~{cnvAnnotatedFinalBedPath}"
        File cnvAnnotatedSupplementalBed = "~{cnvAnnotatedSupplementalBedPath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}

task mergeSv {
    input {
        Int memoryGb = 4
        String pairName
        Array[String] listOfChroms
        String tumor
        String normal
        File filteredMantaSV
        File svabaSv
        IndexedVcf gridssVcf
        String vepGenomeBuild
        Int slop = 300
        String svMergedFinalBedPePath  = "~{pairName}.sv.merged.v7.final.bedpe"
        String svMergedSupplementalBedPePath  = "~{pairName}.sv.merged.v7.supplemental.bedpe"
        Int diskSize = 4
        Int minSvLength = 500

    }

    command {
        Rscript \
        /merge-caller-vcfs.r \
        --vcf=~{filteredMantaSV},~{svabaSv},~{gridssVcf.vcf} \
        --caller=manta,svaba,gridss \
        --tumor=~{tumor} \
        --normal=~{normal} \
        --build=~{vepGenomeBuild} \
        --slop=~{slop} \
        --allowed_chr=~{sep="," listOfChroms} \
        --min_sv_length=~{minSvLength} \
        --out_file=~{svMergedFinalBedPePath} \
        --out_file_supplemental=~{svMergedSupplementalBedPePath}
    }

    output {
        File svMergedFinalBedPe = "~{svMergedFinalBedPePath}"
        File svMergedSupplementalBedPe = "~{svMergedSupplementalBedPePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}

task annotateSv {
    input {
        Int memoryGb = 4
        String pairName
        String tumor
        String normal
        
        Int slop = 500
        
        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe
        
        File svMergedBedPe
        String svMergedAnnotatedBedPePath  = "~{pairName}.sv.merged.annotated.v7.bedpe"
        Int diskSize = 4

    }

    command {        
        Rscript \
        /annotate-bedpe-with-databases.r \
        --db_names=gap,DGV,1000G,PON,COSMIC \
        --db_files=~{gap},~{dgvBedpe},~{thousandGVcf},~{svPon},~{cosmicBedPe} \
        --slop=~{slop} \
        --db_ignore_strand=COSMIC \
        --bedpe=~{svMergedBedPe} \
        --out_file=~{svMergedAnnotatedBedPePath}
    }

    output {
        File svMergedAnnotatedBedPe = "~{svMergedAnnotatedBedPePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}

task annotateGenesSv {
    input {
        Int memoryGb = 4
        String pairName
        String tumor
        String normal
        
        File ensemblUniqueBed
        File cancerCensusBed
        
        File svMergedAnnotatedFinalBedPe
        String svGeneAnnotatedFinalBedPePath  = "~{pairName}.sv.merged.gene.annotated.v7.bedpe"
        Int diskSize = 4

    }

    command {        
        Rscript \
        /annotate-bedpe-with-genes.r \
        --ensembl=~{ensemblUniqueBed} \
        --cancer_census=~{cancerCensusBed} \
        --bedpe=~{svMergedAnnotatedFinalBedPe} \
        --out_file=~{svGeneAnnotatedFinalBedPePath}
    }

    output {
        File svGeneAnnotatedFinalBedPe = "~{svGeneAnnotatedFinalBedPePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}

task annotateGenesSvSupplemental {
    input {
        Int memoryGb = 4
        String pairName
        String tumor
        String normal
        
        File ensemblUniqueBed
        File cancerCensusBed
        
        File svMergedAnnotatedSupplementalBedPe
        String svGeneAnnotatedSupplementalBedPePath  = "~{pairName}.sv.merged.gene.annotated.v7.supplemental.bedpe"
        Int diskSize = 4
    }

    command {        
        Rscript \
        /annotate-bedpe-with-genes.r \
        --ensembl=~{ensemblUniqueBed} \
        --cancer_census=~{cancerCensusBed} \
        --bedpe=~{svMergedAnnotatedSupplementalBedPe} \
        --out_file=~{svGeneAnnotatedSupplementalBedPePath} \
        --supplemental
    }

    output {
        File svGeneAnnotatedSupplementalBedPe = "~{svGeneAnnotatedSupplementalBedPePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}

task annotateWithCnvSv {
    input {
        Int memoryGb = 4
        String pairName
        String tumor
        String normal
        
        File cnvAnnotatedFinalBed
        
        File svGeneAnnotatedBedPe
        String svCnvAnnotatedBedPePath = "~{pairName}.sv.merged.gene.cnv.annotated.v7.bedpe"
        Int diskSize = 4
    }

    command {
        Rscript \
        /annotate-bedpe-with-cnv.r \
        --cnv=~{cnvAnnotatedFinalBed} \
        --bedpe=~{svGeneAnnotatedBedPe} \
        --out_file=~{svCnvAnnotatedBedPePath}
    }

    output {
        File svCnvAnnotatedBedPe = "~{svCnvAnnotatedBedPePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}

task filterBedPe {
    input {
        Int memoryGb = 4
        String pairName
        
        File svCnvAnnotatedBedPe
        String svBedPePath
        String svHighConfidenceBedPePath
        Int diskSize = 4

    }

    command {
        Rscript \
        /filter-bedpe.r \
        --max_changepoint_distance=1000 \
        --filter_databases=DGV,1000G,PON \
        --bedpe=~{svCnvAnnotatedBedPe} \
        --out_file_somatic=~{svBedPePath} \
        --out_file_highconf=~{svHighConfidenceBedPePath}
    }

    output {
        File svBedPe = "~{svBedPePath}"
        File svHighConfidenceBedPe = "~{svHighConfidenceBedPePath}"
    }

    runtime {
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-internal-tools/sv_cnv:1.0.0"
    }
}







