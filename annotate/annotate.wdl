version 1.0

import "../wdl_structs.wdl"

# General tasks


task Vep{
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String vepVersion
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        IndexedVcf vepClinvarXMLTrait
        IndexedTable vepCaddSnp
        IndexedVcf vepGnomadExomes
        IndexedVcf vcfCompressed
        IndexedTable vepCaddIndel
        IndexedVcf vepGnomadGenomes
        IndexedReference vepFastaReference
        String pluginFathmm = "FATHMM,\" python fathmm.py \""
        String pluginFathmmSomatic = "FATHMMSOMATIC,\" python fathmm.py -w Cancer \""
        String pluginCadd = "CADD, ~{vepCaddSnp.table},~{vepCaddIndel.table}"
        String customGnomadExome = "~{vepGnomadExomes.vcf},GnomadExomes,vcf,exact,0,AF,AN,Hom,AFAFR,AFAMR,AFASJ,AFEAS,AFFIN,AFNFE,AFOTH"
        String customGnomadGenome = "~{vepGnomadGenomes.vcf},GnomadGenomes,vcf,exact,0,AF,AN,Hom,AFAFR,AFAMR,AFASJ,AFEAS,AFFIN,AFNFE,AFOTH"
        String customClinvar = "~{vepClinvarXMLTrait.vcf},CLN,vcf,exact,0,VARIATIONID,MOLECULARCONSEQUENCE,CLINICALSIGNIFICANCE,CONFLICTED,REVIEWSTATUS,TRAITS,PMIDS,XREFS,ORIGIN"
        String customCosmicCoding = "~{cosmicCoding.vcf},CosmicCoding,vcf,exact,0,CNT,CDS,AA"
        String customCosmicNoncoding = "~{cosmicNoncoding.vcf},CosmicNonCoding,vcf,exact,0"
        String pairName
        String vcfAnnotatedVepPath = "~{pairName}.v6.vep.annotated.vcf"
    }

    command {
        vep \
        --force_overwrite \
        --buffer_size 10000000 \
        --fork 8 \
        --no_stats \
        --use_transcript_ref \
        --offline \
        --assembly ~{vepVersion} \
        --cache \
        --dir_cache "/opt/vep/.vep" \
        --gencode_basic \
        --exclude_predicted \
        --fasta ~{vepFastaReference.fasta} \
        --sift p \
        --polyphen p \
        --hgvs \
        --symbol \
        --numbers \
        --domains \
        --regulatory \
        --nearest symbol \
        --biotype \
        --tsl \
        --af \
        --af_1kg \
        --appris \
        --gene_phenotype \
        --pubmed \
        --variant_class \
        --vcf \
        --pick_allele \
        --dir_plugins "/opt/vep/.vep/Plugins/" \
        --plugin fathmm.py " \
        --plugin fathmm.py -w Cancer " \
        --plugin ~{pluginCadd} \
        --custom ~{customGnomadExome} \
        --custom ~{customGnomadGenome} \
        --custom ~{customClinvar} \
        --custom ~{customCosmicCoding} \
        --custom ~{customCosmicNoncoding} \
        --input_file ~{vcfCompressed.vcf} \
        --output_file ~{vcfAnnotatedVepPath}
    }

    output {
        File vcfAnnotatedVep = "~{vcfAnnotatedVepPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task AddCosmic {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String vcfAnnotatedCancerGeneCensusPath = "~{pairName}.v6.cosmic_census.vep.annotated.vcf"
        File vcfAnnotatedVep
        File cosmicCensus
    }

    command {
        python \
        add_cancer_gene_census.py \
        ~{cosmicCensus} \
        ~{vcfAnnotatedVep} \
        ~{vcfAnnotatedCancerGeneCensusPath}
    }

    output {
        File vcfAnnotatedCancerGeneCensus = "~{vcfAnnotatedCancerGeneCensusPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task AddCancerResistanceMutations {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String genome
        String pairName
        String vcfAnnotatedResistancePath = "~{pairName}.v6.resistance.vep.annotated.vcf"
        File vcfAnnotatedCancerGeneCensus
        File cancerResistanceMutations
    }

    command {
        python \
        add_cancer_resistance_mutations.py \
        ~{cancerResistanceMutations} \
        ~{genome} \
        ~{vcfAnnotatedCancerGeneCensus} \
        ~{vcfAnnotatedResistancePath}
    }

    output {
        File vcfAnnotatedResistance = "~{vcfAnnotatedResistancePath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task AnnotateId {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String vcfAnnotatedIdPath = "~{pairName}.v6.id.vep.annotated.vcf"
        File vcfAnnotatedResistance
    }

    command {
        python2.7 \
        annotate_id.py \
        ~{vcfAnnotatedResistance} \
        ~{vcfAnnotatedIdPath}
    }

    output {
        File vcfAnnotatedId = "~{vcfAnnotatedIdPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task RenameCsqVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String vcfCsqRenamedPath = "~{pairName}.snv.indel.supplemental.v6.annotated.vcf"
        File vcfAnnotatedId
    }

    command {
        python \
        rename_csq_vcf.py \
        ~{vcfAnnotatedId} \
        ~{vcfCsqRenamedPath} \
    }

    output {
        File vcfCsqRenamed = "~{vcfCsqRenamedPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task MainVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String mainVcfPath = "~{pairName}.snv.indel.final.v6.annotated.vcf"
        File vcfAnnotated
    }

    command {
        python \
        make_main_vcf.py \
        ~{vcfAnnotated} \
        ~{mainVcfPath} \
    }

    output {
        File mainVcf = "~{mainVcfPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task TableVcf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String vcfAnnotatedTxtPath = "~{pairName}.snv.indel.final.v6.annotated.txt"
        String tumor
        String normal
        File mainVcf
    }

    command {
        make_txt.py \
        --vcf ~{mainVcf} \
        --txt ~{vcfAnnotatedTxtPath} \
        --tumor ~{tumor} \
        --normal ~{normal} \
    }

    output {
        File vcfAnnotatedTxt = "~{vcfAnnotatedTxtPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task VcfToMaf {
    input {
        Int threads
        Int memoryGb
        String dockerImage
        String pairName
        String mafPath = "~{pairName}.snv.indel.final.v6.annotated.maf"
        String library
        String vepVersion
        String tumor
        String normal
        File ensemblEntrez
        File mainVcf
    }

    command {
        python \
        make_maf.py \
        --vcf ~{mainVcf} \
        --maf ~{mafPath} \
        --library ~{library} \
        --vep-version ~{vepVersion} \
        --tumor ~{tumor} \
        --normal ~{normal} \
        --ensembl-entrez ~{ensemblEntrez} \
    }

    output {
        File maf = "~{mafPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}

task GermVep {
    input {    
        Int threads
        Int memoryGb
        String dockerImage
        String vepVersion
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        IndexedVcf vepClinvarXMLTrait
        IndexedTable vepCaddSnp
        IndexedVcf vepGnomadExomes
        IndexedVcf vcfCompressed
        IndexedTable vepCaddIndel
        IndexedVcf vepGnomadGenomes
        IndexedReference vepFastaReference
        
        File phylop100
        File phastcons100
        IndexedVcf clinvarXMLDisease
        IndexedVcf acmg
        IndexedVcf pgx
        IndexedVcf arExtended
        IndexedVcf revel
        IndexedVcf ccrs
        IndexedVcf mitimpact
        IndexedVcf mitomapDisease
        IndexedVcf mitomapPolymorphisms
        IndexedVcf mtFuncLoc
        
        String pluginFathmm = "FATHMM,\" python fathmm.py \""
        String pluginFathmmSomatic = "FATHMMSOMATIC,\" python fathmm.py -w Cancer \""
        String custmoCosmicCoding = "~{cosmicCoding.vcf},CosmicCoding,vcf,exact,0,CNT,CDS,AA"
        String customCosmicNoncoding = "~{cosmicNoncoding.vcf},CosmicNonCoding,vcf,exact,0"
        String customClinvarXMLTrait = "~{vepClinvarXMLTrait.vcf},CLN,vcf,exact,0,VARIATIONID,MOLECULARCONSEQUENCE,CLINICALSIGNIFICANCE,CONFLICTED,REVIEWSTATUS,TRAITS,PMIDS,XREFS,ORIGIN"
        String pluginCadd = "CADD, ~{vepCaddSnp.table},~{vepCaddIndel.table}"
        String customGnomadExome = "~{vepGnomadExomes.vcf},GnomadExomes,vcf,exact,0,AF,AN,Hom,AFAFR,AFAMR,AFASJ,AFEAS,AFFIN,AFNFE,AFOTH"
        String customGnomadGenomes = "~{vepGnomadGenomes.vcf},GnomadGenomes,vcf,exact,0,AF,AN,Hom,AFAFR,AFAMR,AFASJ,AFEAS,AFFIN,AFNFE,AFOTH"
        String customPhylop100 = "~{phylop100},phyloP100,bigwig"
        String customPhastcons100 = "~{phastcons100},phastcons100,bigwig"
        String customClinvarXMLDisease = "~{clinvarXMLDisease},CLN,vcf,overlap,0,DiseaseName"
        String customAcmg = "~{acmg},ACMG59,vcf,overlap,0,GENE,DISEASE"
        String customPgx = "~{pgx},PGx,vcf,exact,0,pgxRsid"
        String customArExtended = "~{arExtended},AR,vcf,overlap,0,ARGENE"
        String customRevel = "~{revel},REVEL,vcf,exact,0,REVELSCORE"
        String customCcrs = "~{ccrs},CCRS,vcf,overlap,0,ccrPct"
        String customMitimpact = "~{mitimpact},mitimpact,vcf,exact,0,OXPHOSComplex"
        String customMitomapPolymorphisms = "~{mitomapPolymorphisms},mitomap,vcf,exact,0,AC,AF"
        String customMitomapDisease = "~{mitomapDisease},mitomap,vcf,exact,0,AC,AF,homoplasmy,heteroplasmy,PubmedIDs,Disease,DiseaseStatus"
        String customMtFuncLoc = "~{mtFuncLoc},mtfl,vcf,overlap,0,FUNCLOC"
        String sampleId
        String vcfAnnotatedVepPath = "~{sampleId}.v6.vep.annotated.vcf"
    }

    command {
        vep \
        --force_overwrite \
        --buffer_size 10000000 \
        --fork 8 \
        --no_stats \
        --use_transcript_ref \
        --offline \
        --assembly ~{vepVersion} \
        --cache \
        --dir_cache "/opt/vep/.vep" \
        --gencode_basic \
        --exclude_predicted \
        --fasta ~{vepFastaReference.fasta} \
        --sift p \
        --polyphen p \
        --hgvs \
        --symbol \
        --numbers \
        --domains \
        --regulatory \
        --nearest symbol \
        --biotype \
        --tsl \
        --af \
        --af_1kg \
        --appris \
        --gene_phenotype \
        --pubmed \
        --variant_class \
        --vcf \
        --pick_allele \
        --dir_plugins "/opt/vep/.vep/Plugins/" \
        --plugin ~{pluginFathmm} \
        --plugin ~{pluginFathmmSomatic} \
        --plugin ~{pluginCadd} \
        --custom ~{customGnomadExome} \
        --custom ~{customGnomadGenomes} \
        --custom ~{customClinvarXMLTrait} \
        --custom ~{custmoCosmicCoding} \
        --custom ~{customCosmicNoncoding} \
        --custom ~{customPhylop100} \
        --custom ~{customPhastcons100} \
        --custom ~{customClinvarXMLDisease} \
        --custom ~{customAcmg} \
        --custom ~{customPgx} \
        --custom ~{customArExtended} \
        --custom ~{customRevel} \
        --custom ~{customCcrs} \
        --custom ~{customMitimpact} \
        --custom ~{customMitomapPolymorphisms} \
        --custom ~{customMitomapDisease} \
        --custom ~{customMtFuncLoc} \
        --input_file ~{vcfCompressed.vcf} \
        --output_file ~{vcfAnnotatedVepPath}
    }

    output {
        File vcfAnnotatedVep = "~{vcfAnnotatedVepPath}"
    }

    runtime {
        cpu : threads
        memory : memoryGb + "GB"
        docker : dockerImage
    }
}





