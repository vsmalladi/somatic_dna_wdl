version 1.0


import "calling/mutect2_pon_wkf.wdl" as mutect2Pon
import "calling/manta_pon_wkf.wdl" as mantaPon
import "calling/svaba_pon_wkf.wdl" as svabaPon
import "merge_vcf/merge_vcf_pon_wkf.wdl" as MergeVcfPon
import "annotate/annotate.wdl" as annotate
import "annotate/variantEffectPredictor.wdl" as variantEffectPredictor
import "merge_vcf/merge_vcf.wdl" as mergeVcf

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Nico Robine (nrobine@nygenome.org)
#    Minita Shah (mshah@nygenome.org)
#    Timothy Chu (tchu@nygenome.org)
#    Will Hooper (whooper@nygenome.org)
#
# ================== /COPYRIGHT ===============================================

import "wdl_structs.wdl"

workflow CallingPon {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Boolean production = false
        Array[SampleBamInfo]+ tumorInfos
        # strelka2
        File strelkaJsonLog
        File configureStrelkaSomaticWorkflow
        #   mutect2
        File mutectJsonLog
        Array[String]+ listOfChroms
        IndexedReference referenceFa
        #   Manta
        IndexedTable callRegions
        File mantaJsonLog
        #   Svaba
        File dbsnpIndels
        File svabaJsonLog
        File mutectJsonLogFilter
        BwaReference bwaReference
        Boolean highMem = false
        
        # annotation:
        IndexedVcf cosmicCoding
        IndexedVcf cosmicNoncoding
        
        # Public
        File cancerResistanceMutations
        File vepCache
        File annotations
        File plugins
        String vepGenomeBuild
        IndexedReference vepFastaReference

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
        
        Int vepDiskSize = ceil(size(vepCache, "GB") + size(plugins, "GB") + size(annotations, "GB") + size(hgmdGene.vcf, "GB") + size(hgmdUd10.vcf, "GB") + size(hgmdPro.vcf, "GB") + size(chdGenesVcf.vcf, "GB") + size(chdEvolvingGenesVcf.vcf, "GB") + size(chdWhitelistVcf.vcf, "GB") + size(deepIntronicsVcf.vcf, "GB") + size(clinvarIntronicsVcf.vcf, "GB") + size(masterMind.vcf, "GB")) + 200
        
    }
    scatter(tumorInfo in tumorInfos) {
        
        call mutect2Pon.Mutect2Pon {
            input:
                mutectJsonLogFilter = mutectJsonLogFilter,
                tumor = tumorInfo.sampleId,
                listOfChroms = listOfChroms,
                referenceFa = referenceFa,
                tumorFinalBam = tumorInfo.finalBam,
                highMem = highMem
        }
        
        call mantaPon.MantaPon {
            input:
                mantaJsonLog = mantaJsonLog,
                tumor = tumorInfo.sampleId,
                callRegions = callRegions,
                referenceFa = referenceFa,
                tumorFinalBam = tumorInfo.finalBam,
                highMem = highMem
        }
        
        call svabaPon.SvabaPon {
            input:
                svabaJsonLog = svabaJsonLog,
                tumor = tumorInfo.sampleId,
                dbsnpIndels = dbsnpIndels,
                referenceFa = referenceFa,
                bwaReference = bwaReference,
                tumorFinalBam = tumorInfo.finalBam,
                highMem = highMem
        }
        call MergeVcfPon.MergeVcfPon as wgsMergeVcfPon {
                    input:
                        filteredMantaSV = MantaPon.filteredMantaSV,
                        mutect2 = Mutect2Pon.mutect2,
                        svabaIndel = SvabaPon.svabaIndel,
                        
                        referenceFa = referenceFa,
                        listOfChroms = listOfChroms
        }
        
        call mergeVcf.CompressVcf as unannotatedCompressVcf {
            input:
                vcf = wgsMergeVcfPon.mergedVcf
        }
    
        call mergeVcf.IndexVcf as unannotatedIndexVcf {
            input:
                vcfCompressed = unannotatedCompressVcf.vcfCompressed
        }
    
        if (production) {
            call variantEffectPredictor.vepPublicSvnIndel as productionVepSvnIndel {
                input:
                    pairName = tumorInfo.sampleId,
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
                    diskSize = vepDiskSize,
                    memoryGb = 8
            }
        }
        
        if (!production) {
            if ( size(omimVcf.vcf) > 0 ) {
                call variantEffectPredictor.vepSvnIndel as notProductionVepSvnIndel{
                    input:
                        pairName = tumorInfo.sampleId,
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
                        diskSize = vepDiskSize,
                        memoryGb = 8
                }
            }
        }
    
        File finalVcfAnnotatedVep = select_first([notProductionVepSvnIndel.vcfAnnotatedVep, productionVepSvnIndel.vcfAnnotatedVep])
        
    }
 
    output {
        # General 

        Array[File] vcfAnnotatedVep = finalVcfAnnotatedVep
        
        # Mutect2
        Array[File] mutect2 = Mutect2Pon.mutect2
        # Manta 
        Array[File] filteredMantaSV = MantaPon.filteredMantaSV
        # Strelka2

        # Svaba
        Array[File] svabaIndel = SvabaPon.svabaIndel
    }
}