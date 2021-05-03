version 1.0

import "../merge_vcf/merge_vcf.wdl" as merge_vcf
import "annotate.wdl" as annotate
import "variantEffectPredictor.wdl" as variantEffectPredictor
import "../wdl_structs.wdl"

workflow Annotate {
    input {
        String tumor
        String normal
        String pairName
        File unannotatedVcf
        
        String vepGenomeBuild
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        
        # Public
        File vepCache
        File annotations
        File plugins
        IndexedReference vepFastaReference
        
        # NYGC-only
        IndexedVcf hgmdGene
        IndexedVcf hgmdUd10
        IndexedVcf hgmdPro
        IndexedVcf omimVcf
        
        # Public
        File cancerResistanceMutations
        IndexedVcf chdGenesVcf
        IndexedVcf chdEvolvingGenesVcf
        IndexedVcf chdWhitelistVcf
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf masterMind
        
        # post annotation
        File cosmicCensus
        
        File ensemblEntrez
        String library
    
        IndexedReference referenceFa
        Int vepDiskSize = ceil(size(vepCache, "GB") + size(plugins, "GB") + size(annotations, "GB") + size(hgmdGene.vcf, "GB") + size(hgmdUd10.vcf, "GB") + size(hgmdPro.vcf, "GB") + size(omimVcf.vcf, "GB") + size(chdGenesVcf.vcf, "GB") + size(chdEvolvingGenesVcf.vcf, "GB") + size(chdWhitelistVcf.vcf, "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + size(masterMind.vcf, "GB") + (size(unannotatedVcf, "GB") * 2)) + 500
    }
        
    call merge_vcf.CompressVcf as unannotatedCompressVcf {
        input:
            vcf = unannotatedVcf
    }

    call merge_vcf.IndexVcf as unannotatedIndexVcf {
        input:
            vcfCompressed = unannotatedCompressVcf.vcfCompressed
    }
    
    call variantEffectPredictor.vepSvnIndel {
        input:
            pairName = pairName,
            unannotatedVcf = unannotatedIndexVcf.vcfCompressedIndexed,
            vepCache = vepCache,
            annotations = annotations,
            plugins = plugins,
            vepGenomeBuild = vepGenomeBuild,
            vepFastaReference = vepFastaReference,
            # NYGC-only
            hgmdGene = hgmdGene,
            hgmdUd10 = hgmdUd10,
            hgmdPro = hgmdPro,
            omimVcf = omimVcf,
            # Public
            chdGenesVcf = chdGenesVcf,
            chdEvolvingGenesVcf = chdEvolvingGenesVcf,
            chdWhitelistVcf = chdWhitelistVcf,
            deepIntronicsVcf = deepIntronicsVcf,
            clinvarIntronicsVcf = clinvarIntronicsVcf,
            masterMind = masterMind,
            cosmicCoding = cosmicCoding,
            cosmicNoncoding = cosmicNoncoding,
            diskSize = vepDiskSize
    }
    
    call annotate.AddCosmic {
        input:
            pairName = pairName,
            cosmicCensus = cosmicCensus,
            vcfAnnotatedVep = vepSvnIndel.vcfAnnotatedVep
    }
    
    call annotate.AddCancerResistanceMutations {
        input:
            pairName = pairName,
            cancerResistanceMutations = cancerResistanceMutations,
            vcfAnnotatedCancerGeneCensus = AddCosmic.vcfAnnotatedCancerGeneCensus
    }
    
    call annotate.AnnotateId {
        input:
            pairName = pairName,
            vcfAnnotatedResistance = AddCancerResistanceMutations.vcfAnnotatedResistance
    }
    
    call annotate.RenameCsqVcf {
        input:
            pairName = pairName,
            vcfAnnotatedId = AnnotateId.vcfAnnotatedId
    }
    
    call annotate.MainVcf {
        input:
            pairName = pairName,
            vcfAnnotated = RenameCsqVcf.vcfCsqRenamed
    }
    
    call annotate.TableVcf {
        input:
            tumor = tumor,
            normal = normal,
            pairName = pairName,
            mainVcf = MainVcf.mainVcf
    }
    
    call annotate.VcfToMaf {
        input:
            tumor = tumor,
            normal = normal,
            pairName = pairName,
            mainVcf = MainVcf.mainVcf,
            library = library,
            vepGenomeBuild = vepGenomeBuild,
            ensemblEntrez = ensemblEntrez
    }
    
    output {
        PairVcfInfo pairVcfInfo  = object {
            pairId : "~{pairName}",
            tumor : "~{tumor}",
            normal : "~{normal}",
            mainVcf : "~{MainVcf.mainVcf}",
            supplementalVcf : "~{RenameCsqVcf.vcfCsqRenamed}",
            vcfAnnotatedTxt : "~{TableVcf.vcfAnnotatedTxt}",
            maf : "~{VcfToMaf.maf}",
        }
    }
}