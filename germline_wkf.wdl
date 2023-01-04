version 1.0

import "wdl_structs.wdl"
import "alignment_analysis/alignment_analysis.wdl" as alignmentAnalysis
import "germline/germline_wkf.wdl" as germline
import "annotate/germline_annotate_wkf.wdl" as germlineAnnotate

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

workflow GermlineAll {
    # command 
    input {
        Array[SampleBamInfo]+ normalSampleBamInfos

        IndexedReference referenceFa
        Array[String]+ listOfChroms

        Boolean production = true

        File excludeIntervalList
        Array[File] scatterIntervalsHcs

        IndexedVcf MillsAnd1000G
        IndexedVcf hapmap
        IndexedVcf omni
        IndexedVcf onekG
        IndexedVcf dbsnp

        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf

        # annotation:
        String vepGenomeBuild
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
        
        # post annotation
        File cosmicCensus
        
        File ensemblEntrez
        String library
        
        Boolean highMem = false
    }

    scatter (normalSampleBamInfo in normalSampleBamInfos) {
        # using small disk size because the file is not localized (on servers that support this)
        Int basicDiskSize = 4
        
        call alignmentAnalysis.GetSampleName {
            input:
                finalBam = normalSampleBamInfo.finalBam.bam,
                finalBai = normalSampleBamInfo.finalBam.bamIndex,
                diskSize = basicDiskSize
        }
        
        if (GetSampleName.bamSampleId != normalSampleBamInfo.sampleId ) {
            Int renameDiskSize = (ceil( size(normalSampleBamInfo.finalBam.bam, "GB") )  * 2 ) + 4
            
            call alignmentAnalysis.UpdateBamSampleName {
                input:
                    finalBam = normalSampleBamInfo.finalBam,
                    sampleId = normalSampleBamInfo.sampleId,
                    outputPrefix = normalSampleBamInfo.sampleId,
                    diskSize = renameDiskSize
            }
        }
        
        Bam normalSampleBam = select_first([UpdateBamSampleName.reheaderBam, normalSampleBamInfo.finalBam])
        
        call germline.Germline {
            input:
                finalBam = normalSampleBam,
                normal = normalSampleBamInfo.sampleId,
                outputPrefix = normalSampleBamInfo.sampleId,
                referenceFa = referenceFa,
                listOfChroms = listOfChroms,
                MillsAnd1000G = MillsAnd1000G,
                hapmap = hapmap,
                omni = omni,
                onekG = onekG,
                dbsnp = dbsnp,
                nygcAf = nygcAf,
                excludeIntervalList = excludeIntervalList,
                scatterIntervalsHcs = scatterIntervalsHcs,
                pgx = pgx,
                rwgsPgxBed = rwgsPgxBed,
                whitelist = whitelist,
                chdWhitelistVcf = chdWhitelistVcf,
                deepIntronicsVcf = deepIntronicsVcf,
                clinvarIntronicsVcf = clinvarIntronicsVcf,
                highMem = highMem
        }
        
        call germlineAnnotate.GermlineAnnotate as filteredGermlineAnnotate {
            input:
                unannotatedVcf = Germline.haplotypecallerFinalFiltered,
                production = production,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
                listOfChroms = listOfChroms,
                vepGenomeBuild = vepGenomeBuild,
                cosmicCoding = cosmicCoding,
                cosmicNoncoding = cosmicNoncoding,
                # Public
                cancerResistanceMutations = cancerResistanceMutations,
                vepCache = vepCache,
                annotations = annotations,
                plugins = plugins,
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
                # post annotation
                cosmicCensus = cosmicCensus,
                ensemblEntrez = ensemblEntrez,
                library = library  
        }

        call germlineAnnotate.GermlineAnnotate as unFilteredGermlineAnnotate {
            input:
                unannotatedVcf = Germline.haplotypecallerVcf,
                production = production,
                referenceFa = referenceFa,
                normal = normalSampleBamInfo.sampleId,
                listOfChroms = listOfChroms,
                vepGenomeBuild = vepGenomeBuild,
                cosmicCoding = cosmicCoding,
                cosmicNoncoding = cosmicNoncoding,
                # Public
                cancerResistanceMutations = cancerResistanceMutations,
                vepCache = vepCache,
                annotations = annotations,
                plugins = plugins,
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
                # post annotation
                cosmicCensus = cosmicCensus,
                ensemblEntrez = ensemblEntrez,
                library = library  
        }
    }

    output {
        Array[IndexedVcf] haplotypecallerVcf = Germline.haplotypecallerVcf 
        Array[IndexedVcf] haplotypecallerFinalFiltered = Germline.haplotypecallerFinalFiltered 
        Array[File] filteredHaplotypecallerAnnotatedVcf = filteredGermlineAnnotate.haplotypecallerAnnotatedVcf
        Array[File] haplotypecallerAnnotatedVcf = unFilteredGermlineAnnotate.haplotypecallerAnnotatedVcf
    }
    
}