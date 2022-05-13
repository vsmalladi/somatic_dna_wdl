version 1.0

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

workflow ancestry {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        Array[SampleBamInfo]+ normalSampleBamInfos
        
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
    
    scatter (normalSampleBamInfo in normalSampleBamInfos) {
		String normalSampleId = normalSampleBamInfo.sampleId

        call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixContinental{
            input:
                normalFinalBam = normalSampleBamInfo.finalBam,
                fastNgsAdmixSites = fastNgsAdmixContinentalSites,
                fastNgsAdmixSitesBin = fastNgsAdmixContinentalSitesBin,
                fastNgsAdmixSitesIdx = fastNgsAdmixContinentalSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixRef = fastNgsAdmixContinentalRef,
                fastNgsAdmixNind = fastNgsAdmixContinentalNind,
                outprefix = normalSampleId
        }

        call fastNgsAdmix.FastNgsAdmix as fastNgsAdmixPopulation{
            input:
                normalFinalBam = normalSampleBamInfo.finalBam,
                fastNgsAdmixSites = fastNgsAdmixPopulationSites,
                fastNgsAdmixSitesBin = fastNgsAdmixPopulationSitesBin,
                fastNgsAdmixSitesIdx = fastNgsAdmixPopulationSitesIdx,
                fastNgsAdmixChroms = fastNgsAdmixChroms,
                fastNgsAdmixRef = fastNgsAdmixPopulationRef,
                fastNgsAdmixNind = fastNgsAdmixPopulationNind,
                outprefix = normalSampleId
        }
    }
    
    output {
		# ancestry
        Array[File] beagleFileContinental = fastNgsAdmixContinental.beagleFile
        Array[File] fastNgsAdmixQoptContinental = fastNgsAdmixContinental.fastNgsAdmixQopt
        Array[File] beagleFilePopulation = fastNgsAdmixPopulation.beagleFile
        Array[File] fastNgsAdmixQoptPopulation = fastNgsAdmixPopulation.fastNgsAdmixQopt
    }
}