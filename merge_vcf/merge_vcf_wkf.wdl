version 1.0

import "prep_merge_vcf_wkf.wdl" as prepMergeVcf
import "merge_callers_wkf.wdl" as mergeCallers
import "merge_chroms_wkf.wdl" as mergeChroms
import "../wdl_structs.wdl"

workflow MergeVcf {
    # command
    #   Call variants in BAMs
    #   merge and filter raw VCFs
    #   annotate
    input {
        PreMergedPairVcfInfo preMergedPairVcfInfo
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        
        # merge callers
        File intervalListBed
        File ponFile
        IndexedVcf germFile
    }
    
    call prepMergeVcf.PrepMergeVcf as filteredMantaSVPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.filteredMantaSV,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='manta',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as strelka2SnvPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.strelka2Snv,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='strelka2',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as strelka2IndelPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.strelka2Indel,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='strelka2',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as mutect2PrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.mutect2,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='mutect2',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as lancetPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.lancet,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='lancet',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call prepMergeVcf.PrepMergeVcf as svabaIndelPrepMergeVcf {
        input:
            callerVcf=preMergedPairVcfInfo.svabaIndel,
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            tool='svaba',
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa
            
    }
    
    call mergeCallers.MergeCallers {
        input:
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            pairName=preMergedPairVcfInfo.pairId,
            listOfChroms=listOfChroms,
            intervalListBed=intervalListBed,
            referenceFa=referenceFa,
            normalFinalBam=preMergedPairVcfInfo.normalFinalBam,
            tumorFinalBam=preMergedPairVcfInfo.tumorFinalBam,
            ponFile=ponFile,
            germFile=germFile,
            allVcfCompressed=[filteredMantaSVPrepMergeVcf.preppedVcf, 
                strelka2SnvPrepMergeVcf.preppedVcf,
                strelka2IndelPrepMergeVcf.preppedVcf,
                mutect2PrepMergeVcf.preppedVcf,
                lancetPrepMergeVcf.preppedVcf,
                svabaIndelPrepMergeVcf.preppedVcf]
        
    }
    
    call mergeChroms.MergeChroms {
        input:
            tumor=preMergedPairVcfInfo.tumor,
            normal=preMergedPairVcfInfo.normal,
            pairName=preMergedPairVcfInfo.pairId,
            referenceFa=referenceFa,
            finalChromVcf=MergeCallers.finalChromVcf,
            
    }
    
    output {
        File mergedVcf = MergeChroms.unannotatedVcf
    }
}