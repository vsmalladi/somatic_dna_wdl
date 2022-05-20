version 1.0

import "alignment_analysis/kourami_wfk.wdl" as kourami
import "alignment_analysis/msi_wkf.wdl" as msi
import "alignment_analysis/fastngsadmix_wkf.wdl" as fastNgsAdmix
import "wdl_structs.wdl"

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

workflow AlignmentAnalysis {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[pairInfo]+ pairInfos
        
        # kourami
        BwaReference kouramiReference
        File kouramiFastaGem1Index
        
        # mantis
        File mantisBed
        File intervalListBed
        IndexedReference referenceFa

        #fastNgsAdmix
        File fastNgsAdmixChroms

        File fastNgsAdmixContinentalSites
        File fastNgsAdmixContinentalSitesBin
        File fastNgsAdmixContinentalSitesIdx
        File fastNgsAdmixContinentalRef
        File fastNgsAdmixContinentalNind

        File fastNgsAdmixPopulationSites
        File fastNgsAdmixPopulationSitesBin
        File fastNgsAdmixPopulationSitesIdx
        File fastNgsAdmixPopulationRef
        File fastNgsAdmixPopulationNind
    }
    
    scatter(pairInfo in pairInfos) {
        call kourami.Kourami as kouramiNormal {
            input:
                sampleId=pairInfo.normal,
                kouramiReference=kouramiReference,
                finalBam=pairInfo.normalFinalBam,
                kouramiFastaGem1Index=kouramiFastaGem1Index
        }
        
        call kourami.Kourami as kouramiTumor {
            input:
                sampleId=pairInfo.tumor,
                kouramiReference=kouramiReference,
                finalBam=pairInfo.tumorFinalBam,
                kouramiFastaGem1Index=kouramiFastaGem1Index
        }
        
        call msi.Msi {
            input:
                normal=pairInfo.normal,
                pairName=pairInfo.pairId,
                mantisBed=mantisBed,
                intervalListBed=intervalListBed,
                referenceFa=referenceFa,
                tumorFinalBam=pairInfo.tumorFinalBam,
                normalFinalBam=pairInfo.normalFinalBam
        }

        call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixContinental{
            input:
                normalFinalBam = pairInfo.normalFinalBam,
                fastNgsAdmixSites = fastNgsAdmixContinentalSites,
                fastNgsAdmixSitesBin = fastNgsAdmixContinentalSitesBin,
                fastNgsAdmixSitesIdx = fastNgsAdmixContinentalSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixRef = fastNgsAdmixContinentalRef,
                fastNgsAdmixNind = fastNgsAdmixContinentalNind,
                outprefix = pairInfo.normal + "_continental"
        }

        call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixPopulation{
            input:
                normalFinalBam = pairInfo.normalFinalBam,
                fastNgsAdmixSites = fastNgsAdmixPopulationSites,
                fastNgsAdmixSitesBin = fastNgsAdmixPopulationSitesBin,
                fastNgsAdmixSitesIdx = fastNgsAdmixPopulationSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixRef = fastNgsAdmixPopulationRef,
                fastNgsAdmixNind = fastNgsAdmixPopulationNind,
                outprefix = pairInfo.normal + "_population"
        }
    }
    
    output {
        Array[File] normalKouramiResult = kouramiNormal.result
        Array[File] tumorKouramiResult = kouramiTumor.result
        Array[File] mantisWxsKmerCountsFinal = Msi.mantisWxsKmerCountsFinal
        Array[File] mantisWxsKmerCountsFiltered = Msi.mantisWxsKmerCountsFiltered
        Array[File] mantisExomeTxt = Msi.mantisExomeTxt
        Array[File] mantisStatusFinal = Msi.mantisStatusFinal
        Array[File] beagleFileContinental = fastNgsAdmixContinental.beagleFile
        Array[File] fastNgsAdmixQoptContinental = fastNgsAdmixContinental.fastNgsAdmixQopt
        Array[File] beagleFilePopulation = fastNgsAdmixPopulation.beagleFile
        Array[File] fastNgsAdmixQoptPopulation = fastNgsAdmixPopulation.fastNgsAdmixQopt
    }
}