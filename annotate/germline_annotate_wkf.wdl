version 1.0

import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "annotate.wdl" as annotate
import "variantEffectPredictor.wdl" as variantEffectPredictor
import "../wdl_structs.wdl"

workflow GermlineAnnotate {
    input {
        String normal
        String sampleId = "~{normal}"
        IndexedVcf unannotatedVcf
        Boolean production = true
        
        Array[String]+ listOfChroms
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
        Int vepDiskSize = ceil(size(vepCache, "GB") + size(plugins, "GB") + size(annotations, "GB") + size(hgmdGene.vcf, "GB") + size(hgmdUd10.vcf, "GB") + size(hgmdPro.vcf, "GB") + size(omimVcf.vcf, "GB") + size(chdGenesVcf.vcf, "GB") + size(chdEvolvingGenesVcf.vcf, "GB") + size(chdWhitelistVcf.vcf, "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + size(masterMind.vcf, "GB") + (size(unannotatedVcf.vcf, "GB") * 2)) + 500
    }
    
    # split multi alleleic sites and remove HLA contig calls
    call mergeVcf.SplitMultiAllelicRegions {
            input:
                pairName = sampleId,
                vcfCompressedIndexed = unannotatedVcf,
                listOfChroms = listOfChroms,
                splitVcfPath = sub(basename(unannotatedVcf.vcf), ".vcf.gz$", ".split.vcf"),
                referenceFa = referenceFa
        }
        
    call mergeVcf.CompressVcf {
        input:
            vcf = SplitMultiAllelicRegions.sortedVcf
    }
    
    call mergeVcf.IndexVcf {
        input:
            vcfCompressed = CompressVcf.vcfCompressed
    }
    
    if (production) {
        call variantEffectPredictor.vepPublicSvnIndel as productionVepSvnIndel {
            input:
                pairName = pairName,
                unannotatedVcf = unannotatedIndexVcf.vcfCompressedIndexed,
                vepCache = vepCache,
                annotations = annotations,
                plugins = plugins,
                vepGenomeBuild = vepGenomeBuild,
                vepFastaReference = vepFastaReference,
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
    }
    
    if (!production) {
        if ( size(omimVcf.vcf) > 0 ) {
            call variantEffectPredictor.vepSvnIndel as notProductionVepSvnIndel{
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
        }
    }
    
    File vcfAnnotatedVep = select_first([notProductionVepSvnIndel.vcfAnnotatedVep, productionVepSvnIndel.vcfAnnotatedVep])
    
    call annotate.RemoveSpanning {
        input:
            sampleId = sampleId,
            vcfAnnotatedVep = vcfAnnotatedVep
    }
    
    call annotate.AddCosmic {
        input:
            pairName = sampleId,
            cosmicCensus = cosmicCensus,
            vcfAnnotatedVep = RemoveSpanning.noSpanningVcf
    }
    
    call annotate.AddCancerResistanceMutations {
        input:
            pairName = sampleId,
            cancerResistanceMutations = cancerResistanceMutations,
            vcfAnnotatedCancerGeneCensus = AddCosmic.vcfAnnotatedCancerGeneCensus
    }
    
    call annotate.AnnotateId {
        input:
            pairName = sampleId,
            vcfAnnotatedResistance = AddCancerResistanceMutations.vcfAnnotatedResistance
    }
    
    call annotate.RenameCsqVcf {
        input:
            pairName = sampleId,
            vcfAnnotatedId = AnnotateId.vcfAnnotatedId,
            vcfCsqRenamedPath = "~{sampleId}.haplotypecaller.gatk.v4.1.8.0.annotated.vcf"
    }
    
    
    output {
    
        File haplotypecallerAnnotatedVcf = RenameCsqVcf.vcfCsqRenamed
    }
}