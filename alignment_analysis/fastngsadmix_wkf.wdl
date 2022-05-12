version 1.0

import "alignment_analysis.wdl" as alignmentAnalysis
import "../wdl_structs.wdl"

workflow FastNGSadmix {
    input {
        Bam normalFinalBam
        File fastngsadmixSites
        File fastngsadmixSitesBin
        File fastngsadmixSitesIdx
        File fastngsadmixChroms
        File fastngsadmixRef
        File fastngsadmixNind
        Int angsd_threads = 4 
        String outprefix
    }

    call alignmentAnalysis.Angsd {
        input:
            normalFinalBam = normalFinalBam,
            fastngsadmixSites = fastngsadmixSites,
            fastngsadmixSitesBin = fastngsadmixSitesBin,
            fastngsadmixSitesIdx = fastngsadmixSitesIdx,
            fastngsadmixChroms = fastngsadmixChroms,
            threads = angsd_threads,
            outprefix = outprefix
    }

    call alignmentAnalysis.FastNGSadmix {
        input:
            beagleFile = Angsd.beagleFile,
            fastngsadmixRef = fastngsadmixRef,
            fastngsadmixNind = fastngsadmixNind,
            outprefix = outprefix
    }

    output {
        File beagleFile = Angsd.beagleFile
        File fastngsadmixQopt = FastNGSadmix.fastngsadmixQopt
    }
}