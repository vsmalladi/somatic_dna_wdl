version 1.0

import "alignment_analysis/fastngsadmix_wkf.wdl" as fastngsadmix
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
        
        #fastngsadmix
        File fastngsadmixChroms

        File fastngsadmixContinentalSites
        File fastngsadmixContinentalSitesBin
        File fastngsadmixContinentalSitesIdx
        File fastngsadmixContinentalRef
        File fastngsadmixContinentalNind

        File fastngsadmixPopulationSites
        File fastngsadmixPopulationSitesBin
        File fastngsadmixPopulationSitesIdx
        File fastngsadmixPopulationRef
        File fastngsadmixPopulationNind

    }
    
    scatter (normalSampleBamInfo in normalSampleBamInfos) {
		String normalSampleIds = normalSampleBamInfo.sampleId
		
        call fastngsadmix.FastNGSadmix as fastngsadmixContinental{
            input:
                bam = normalSampleBamInfo.finalBam,
                fastngsadmixSites = fastngsadmixContinentalSites,
                fastngsadmixSitesBin = fastngsadmixContinentalSitesBin,
                fastngsadmixSitesIdx = fastngsadmixContinentalSitesIdx,
                fastngsadmixChroms = fastngsadmixChroms,
                fastngsadmixRef = fastngsadmixContinentalRef,
                fastngsadmixNind = fastngsadmixContinentalNind,
                outprefix = normalSampleIds
        }

        call fastngsadmix.FastNGSadmix as fastngsadmixPopulation{
            input:
                bam = normalSampleBamInfo.finalBam,
                fastngsadmixSites = fastngsadmixPopulationSites,
                fastngsadmixSitesBin = fastngsadmixPopulationSitesBin,
                fastngsadmixSitesIdx = fastngsadmixPopulationSitesIdx,
                fastngsadmixChroms = fastngsadmixChroms,
                fastngsadmixRef = fastngsadmixPopulationRef,
                fastngsadmixNind = fastngsadmixPopulationNind,
                outprefix = normalSampleIds
        }
    }
    
    output {
		# ancestry
        Array[File] beagleFileContinental = fastngsadmixContinental.beagleFile
        Array[File] fastngsadmixQoptContinental = fastngsadmixContinental.fastngsadmixQopt
        Array[File] beagleFilePopulation = fastngsadmixPopulation.beagleFile
        Array[File] fastngsadmixQoptPopulation = fastngsadmixPopulation.fastngsadmixQopt
    }
}