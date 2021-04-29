version 1.0

import "../wdl_structs.wdl"

task vepSvnIndel {
    input {
        IndexedReference vepFastaReference
        String pairName
        IndexedVcf unannotatedVcf
        
        # Somatic
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        
        # Public
        File vepCache
        File annotations
        File plugins
        String vepGenomeBuild
        
        # NYGC-only
        IndexedVcf hgmdGene
        IndexedVcf hgmdUd10
        IndexedVcf hgmdPro
        IndexedVcf omimVcf
        
        # Public
        IndexedVcf chdGenesVcf
        IndexedVcf chdEvolvingGenesVcf
        IndexedVcf chdWhitelistVcf
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf masterMind

        String vcfAnnotatedVepPath = "~{pairName}.v7.vep.annotated.vcf"
        Int threads = 16
        Int disk_size
        
    }
    
    command {
        set -e -o pipefail
        
        # NOTE:task will not work with any other genome build as is
        # because of this section
        
        mkdir -p ensembl_vep/homo_sapiens_refseq
        tar -xzvf ~{vepCache}
        mv 97_GRCh38 ensembl_vep/homo_sapiens_refseq/
        tar -xzvf ~{annotations}
        mv annotations ensembl_vep/
        tar -xzvf ~{plugins}
        mv gpfs/commons/groups/clinical/ensembl_vep/Plugins ensembl_vep/

        /opt/vep/src/ensembl-vep/vep \
        --fork ~{threads} \
        --buffer_size 50000 \
        --format vcf \
        --no_stats \
        --no_escape \
        --offline \
        --assembly ~{vepGenomeBuild} \
        --cache \
        --dir_cache ensembl_vep \
        --refseq \
        --exclude_predicted \
        --fasta ~{vepFastaReference.fasta} \
        --symbol \
        --hgvs \
        --check_existing \
        --vcf \
        --pick_allele_gene \
        --dir_plugins ensembl_vep/Plugins \
        --plugin dbscSNV,ensembl_vep/Plugins/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
        --plugin MaxEntScan,ensembl_vep/Plugins/fordownload \
        --plugin dbNSFP,ensembl_vep/Plugins/dbNSFP4.0a/dbNSFP4.0a.gz,ensembl_vep/Plugins/dbNSFP_replacement_logic,REVEL_score,SIFT_pred,SIFT4G_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,MetaSVM_pred,PrimateAI_pred,fathmm-MKL_coding_pred,GERP++_RS,phyloP100way_vertebrate,CADD_phred,Polyphen2_HVAR_pred \
        --custom ensembl_vep/annotations/04142020_NYGC_samples.vcf.gz,NYGC,vcf,exact,0,AF,Samples,AC_Het,AC_Hom \
        --custom ensembl_vep/annotations/clinvar.vep.vcf.gz,CLN_Overlap,vcf,overlap,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
        --custom ensembl_vep/annotations/clinvar.vep.vcf.gz,CLN_Exact,vcf,exact,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
        --custom ensembl_vep/annotations/gnomad_exomes_subset_final.vcf.gz,GnomadExomes,vcf,exact,0,AF,nhomalt \
        --custom ensembl_vep/annotations/gnomad_genomes_subset_final.vcf.gz,GnomadGenomes,vcf,exact,0,AF,nhomalt \
        --custom ~{cosmicCoding.vcf},CosmicCoding,vcf,exact,0,GENOMIC_ID,LEGACY_ID,CNT,CDS,AA
        --custom ~{cosmicNoncoding.vcf},CosmicNonCoding,vcf,exact,0,GENOMIC_ID,LEGACY_ID,CNT,CDS,AA
        --custom ~{hgmdPro.vcf},HGMD,vcf,exact,0,CLASS,PHEN,DNA \
        --custom ~{hgmdUd10.vcf},HGMDUD10,vcf,overlap,0,CLASS,PHEN \
        --custom ~{hgmdGene.vcf},_,vcf,overlap,0,HGMD_PHENO \
        --custom ~{omimVcf.vcf},OMIM,vcf,overlap,0,PHENO,INHT \
        --custom ~{chdGenesVcf.vcf},CHD_GENES,vcf,overlap,0,GENE \
        --custom ~{chdEvolvingGenesVcf.vcf},CHD_EVOLVING,vcf,overlap,0,GENE \
        --custom ~{deepIntronicsVcf.vcf},INTRONIC,vcf,exact,0,INTRONIC \
        --custom ~{clinvarIntronicsVcf.vcf},CLINVAR_INTRONIC,vcf,exact,0,INTRONIC \
        --custom ensembl_vep/annotations/spliceai_scores.hg38.sorted.vcf.gz,SPLICEAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL \
        --custom ensembl_vep/annotations/pli_hg38.vcf.gz,PLI,vcf,overlap,0,pLI,mis_z \
        --custom ensembl_vep/annotations/domino_genes_38.vcf.gz,Domino,vcf,overlap,0,Domino_Score \
        --custom ensembl_vep/annotations/ar_extended.vcf.gz,AR,vcf,overlap,0,AR_GENE \
        --custom ensembl_vep/annotations/ACMG59_2017-09-28.vcf.gz,ACMG59,vcf,overlap,0,GENE,DISEASE \
        --custom ensembl_vep/annotations/dials_genes_b38.vcf.gz,DIALS,vcf,overlap,0,DIALS_GENE \
        --custom ensembl_vep/annotations/pgx_vep.vcf.gz,PGx,vcf,exact,0,pgx_rsid \
        --custom ~{chdWhitelistVcf.vcf},chd_whitelist,vcf,overlap,0,END \
        --custom ~{masterMind.vcf},mm,vcf,exact,0,GENE,HGVSG,MMCNT1,MMCNT2,MMCNT3,MMID3,MMURI3 \
        --custom ensembl_vep/annotations/sema4_immuno_genes_b38.vcf.gz,IMMUNO,vcf,overlap,0,IMMUNO_Gene \
        --custom ensembl_vep/annotations/sema4_neuro_genes_b38.vcf.gz,NEURO,vcf,overlap,0,NEURO_Gene \
        --custom ensembl_vep/annotations/sema4_cardio_genes_b38.vcf.gz,CARDIO,vcf,overlap,0,CARDIO_Gene \
        --custom ensembl_vep/annotations/nygc_curation_b38.vcf.gz,N19,vcf,overlap,0,NYGC_CUR \
        --custom ensembl_vep/annotations/nygc_reported_variants_b38.vcf.gz,R19,vcf,overlap,0,NYGC_REPORTED_SAMPLE,NYGC_CLASS,NYGC_DISEASE \
        --input_file ~{unannotatedVcf.vcf} \
        --output_file ~{vcfAnnotatedVepPath}

    }

    runtime {
        memory: "32 GB"
        cpu: "16"
        disks: "local-disk " + disk_size + " SSD"
        docker: "ensemblorg/ensembl-vep:release_97.4"
    }

    output {
        File vcfAnnotatedVep = "~{vcfAnnotatedVepPath}"
    }
}

task vepGermSvnIndel {
    input {
        IndexedReference vepFastaReference
        String sampleId
        File unannotatedVcf
        
        # Public
        File vepCache
        File annotations
        File plugins
        String vepGenomeBuild
        
        # Somatic
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        
        # NYGC-only
        IndexedVcf hgmdGene
        IndexedVcf hgmdUd10
        IndexedVcf hgmdPro
        IndexedVcf omimVcf
        
        # Public
        IndexedVcf chdGenesVcf
        IndexedVcf chdEvolvingGenesVcf
        IndexedVcf chdWhitelistVcf
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf masterMind

        String vcfAnnotatedVepPath = "~{sampleId}.v7.vep.annotated.vcf"
        Int threads = 16
        Int disk_size
        
    }
    
    command {
        set -e -o pipefail
        
        # NOTE:task will not work with any other genome build as is
        # because of this section
        
        mkdir -p ensembl_vep/homo_sapiens_refseq
        tar -xzvf ~{vepCache}
        mv 97_GRCh38 ensembl_vep/homo_sapiens_refseq/
        tar -xzvf ~{annotations}
        mv annotations ensembl_vep/
        tar -xzvf ~{plugins}
        mv gpfs/commons/groups/clinical/ensembl_vep/Plugins ensembl_vep/

        /opt/vep/src/ensembl-vep/vep \
        --fork ~{threads} \
        --buffer_size 50000 \
        --format vcf \
        --no_stats \
        --no_escape \
        --offline \
        --assembly ~{vepGenomeBuild} \
        --cache \
        --dir_cache ensembl_vep \
        --refseq \
        --exclude_predicted \
        --fasta ~{vepFastaReference.fasta} \
        --symbol \
        --hgvs \
        --check_existing \
        --vcf \
        --pick_allele_gene \
        --dir_plugins ensembl_vep/Plugins \
        --plugin dbscSNV,ensembl_vep/Plugins/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
        --plugin MaxEntScan,ensembl_vep/Plugins/fordownload \
        --plugin dbNSFP,ensembl_vep/Plugins/dbNSFP4.0a/dbNSFP4.0a.gz,ensembl_vep/Plugins/dbNSFP_replacement_logic,REVEL_score,SIFT_pred,SIFT4G_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,MetaSVM_pred,PrimateAI_pred,fathmm-MKL_coding_pred,GERP++_RS,phyloP100way_vertebrate,CADD_phred,Polyphen2_HVAR_pred \
        --custom ensembl_vep/annotations/04142020_NYGC_samples.vcf.gz,NYGC,vcf,exact,0,AF,Samples,AC_Het,AC_Hom \
        --custom ensembl_vep/annotations/clinvar.vep.vcf.gz,CLN_Overlap,vcf,overlap,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
        --custom ensembl_vep/annotations/clinvar.vep.vcf.gz,CLN_Exact,vcf,exact,0,CLIN_ID,CLNSIG,CLNREVSTAT,CLNDN \
        --custom ensembl_vep/annotations/gnomad_exomes_subset_final.vcf.gz,GnomadExomes,vcf,exact,0,AF,nhomalt \
        --custom ensembl_vep/annotations/gnomad_genomes_subset_final.vcf.gz,GnomadGenomes,vcf,exact,0,AF,nhomalt \
        --custom ~{cosmicCoding.vcf},CosmicCoding,vcf,exact,0,GENOMIC_ID,LEGACY_ID,CNT,CDS,AA
        --custom ~{cosmicNoncoding.vcf},CosmicNonCoding,vcf,exact,0,GENOMIC_ID,LEGACY_ID,CNT,CDS,AA
        --custom ~{hgmdPro.vcf},HGMD,vcf,exact,0,CLASS,PHEN,DNA \
        --custom ~{hgmdUd10.vcf},HGMDUD10,vcf,overlap,0,CLASS,PHEN \
        --custom ~{hgmdGene.vcf},_,vcf,overlap,0,HGMD_PHENO \
        --custom ~{omimVcf.vcf},OMIM,vcf,overlap,0,PHENO,INHT \
        --custom ~{chdGenesVcf.vcf},CHD_GENES,vcf,overlap,0,GENE \
        --custom ~{chdEvolvingGenesVcf.vcf},CHD_EVOLVING,vcf,overlap,0,GENE \
        --custom ~{deepIntronicsVcf.vcf},INTRONIC,vcf,exact,0,INTRONIC \
        --custom ~{clinvarIntronicsVcf.vcf},CLINVAR_INTRONIC,vcf,exact,0,INTRONIC \
        --custom ensembl_vep/annotations/spliceai_scores.hg38.sorted.vcf.gz,SPLICEAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL \
        --custom ensembl_vep/annotations/pli_hg38.vcf.gz,PLI,vcf,overlap,0,pLI,mis_z \
        --custom ensembl_vep/annotations/domino_genes_38.vcf.gz,Domino,vcf,overlap,0,Domino_Score \
        --custom ensembl_vep/annotations/ar_extended.vcf.gz,AR,vcf,overlap,0,AR_GENE \
        --custom ensembl_vep/annotations/ACMG59_2017-09-28.vcf.gz,ACMG59,vcf,overlap,0,GENE,DISEASE \
        --custom ensembl_vep/annotations/eye_disease_genes_b38.vcf.gz,EYE_DISEASE,vcf,overlap,0,EYE_DISEASE_GENE \
        --custom ensembl_vep/annotations/dials_genes_b38.vcf.gz,DIALS,vcf,overlap,0,DIALS_GENE \
        --custom ensembl_vep/annotations/renal_disease_genes_b38.vcf.gz,RENAL_DISEASE,vcf,overlap,0,RENAL_DISEASE_GENE \
        --custom ensembl_vep/annotations/pgx_vep.vcf.gz,PGx,vcf,exact,0,pgx_rsid \
        --custom ~{chdWhitelistVcf.vcf},chd_whitelist,vcf,overlap,0,END \
        --custom ~{masterMind.vcf},mm,vcf,exact,0,GENE,HGVSG,MMCNT1,MMCNT2,MMCNT3,MMID3,MMURI3 \
        --custom ensembl_vep/annotations/sema4_immuno_genes_b38.vcf.gz,IMMUNO,vcf,overlap,0,IMMUNO_Gene \
        --custom ensembl_vep/annotations/sema4_neuro_genes_b38.vcf.gz,NEURO,vcf,overlap,0,NEURO_Gene \
        --custom ensembl_vep/annotations/sema4_cardio_genes_b38.vcf.gz,CARDIO,vcf,overlap,0,CARDIO_Gene \
        --custom ensembl_vep/annotations/nygc_curation_b38.vcf.gz,N19,vcf,overlap,0,NYGC_CUR \
        --custom ensembl_vep/annotations/nygc_reported_variants_b38.vcf.gz,R19,vcf,overlap,0,NYGC_REPORTED_SAMPLE,NYGC_CLASS,NYGC_DISEASE \
        --input_file ~{unannotatedVcf} \
        --output_file ~{vcfAnnotatedVepPath}

    }

    runtime {
        memory: "32 GB"
        cpu: "16"
        disks: "local-disk " + disk_size + " SSD"
        docker: "ensemblorg/ensembl-vep:release_97.4"
    }

    output {
        File vcfAnnotatedVep = "~{vcfAnnotatedVepPath}"
    }
}

