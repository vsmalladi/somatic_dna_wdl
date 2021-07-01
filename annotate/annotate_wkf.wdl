version 1.0

import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "annotate.wdl" as annotate
import "variantEffectPredictor.wdl" as variantEffectPredictor
import "../wdl_structs.wdl"

workflow Annotate {
    input {
        String tumor
        String normal
        String pairName
        File unannotatedVcf

        Boolean production = true

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
        Int vepDiskSize = ceil(size(vepCache, "GB") + size(plugins, "GB") + size(annotations, "GB") + size(hgmdGene.vcf, "GB") + size(hgmdUd10.vcf, "GB") + size(hgmdPro.vcf, "GB") + size(chdGenesVcf.vcf, "GB") + size(chdEvolvingGenesVcf.vcf, "GB") + size(chdWhitelistVcf.vcf, "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + size(masterMind.vcf, "GB") + (size(unannotatedVcf, "GB") * 2)) + 500

    }

    call mergeVcf.CompressVcf as unannotatedCompressVcf {
        input:
            vcf = unannotatedVcf
    }

    call mergeVcf.IndexVcf as unannotatedIndexVcf {
        input:
            vcfCompressed = unannotatedCompressVcf.vcfCompressed
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
                # NYGC-only
                hgmdGene = hgmdGene,
                hgmdUd10 = hgmdUd10,
                hgmdPro = hgmdPro,
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

    call annotate.AddCosmic {
        input:
            pairName = pairName,
            cosmicCensus = cosmicCensus,
            vcfAnnotatedVep = vcfAnnotatedVep
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
