version 1.0

import "../wdl_structs.wdl"
import "germline.wdl" as germline
import "../merge_vcf/merge_vcf.wdl"
import "../calling/calling.wdl"

workflow Germline {
    # command 
    input {
        Bam finalBam
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        String normal
        
        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedReference referenceFa
        IndexedVcf onekG
        IndexedVcf dbsnp
        
        Int hcDiskSize = ceil( size(finalBam.bam, "GB") ) + 20
    }
    scatter(chrom in listOfChroms) {
        call germline.Haplotypecaller {
            input:
                referenceFa = referenceFa,
                finalBam = finalBam,
                chrom = chrom,
                sampleId = normal,
                diskSize = hcDiskSize
        }
    }
    
    call calling.Gatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{normal}.haplotypeCalls.er.raw.vcf",
            tempChromVcfs = Haplotypecaller.haplotypecallerChromVcf,
            referenceFa = referenceFa,
            memoryGb = 8,
            diskSize = 10
    }
    
    call merge_vcf.CompressVcf as haplotypecallerCompressVcf {
        input:
            vcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            memoryGb = 4
    }
    
    call merge_vcf.IndexVcf as haplotypecallerIndexVcf {
        input:
            vcfCompressed = haplotypecallerCompressVcf.vcfCompressed
    }
    
    call germline.GentotypeGvcfs {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            sortedVcf = haplotypecallerIndexVcf.vcfCompressedIndexed
    }
    
    call germline.RecalVcfsSnp {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            omni = omni,
            hapmap = hapmap,
            onekG = onekG,
            dbsnp = dbsnp,
            haplotypecallerGenoVcf = GentotypeGvcfs.haplotypecallerGenoVcf
    }
    
    call germline.RecalVcfsIndel {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            dbsnp = dbsnp,
            MillsAnd1000G = MillsAnd1000G,
            haplotypecallerGenoVcf = GentotypeGvcfs.haplotypecallerGenoVcf  
    }
    
    call germline.ApplyRecal as snpApplyRecal {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndelPath = "~{normal}.SNP.apply_vqsr.vcf",
            haplotypecallerGenoVcfApply = GentotypeGvcfs.haplotypecallerGenoVcf,
            mode = "SNP",
            tranches = RecalVcfsSnp.tranchesSnp,
            recal = RecalVcfsSnp.recalSnp    
    }
    
    call germline.ApplyRecal as indelApplyRecal {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndelPath = "~{normal}.SNP.INDEL.apply_vqsr.vcf",
            haplotypecallerGenoVcfApply = snpApplyRecal.haplotypecallerGenoVcfApplySnpIndel,
            mode = "INDEL",
            tranches = RecalVcfsIndel.tranchesIndel,
            recal = RecalVcfsIndel.recalIndel   
    }
    
    call germline.VarFilter {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndel = indelApplyRecal.haplotypecallerGenoVcfApplySnpIndel
    }
    
    call merge_vcf.CompressVcf as recalCompressVcf {
        input:
            vcf = VarFilter.haplotypecallerRecalVcf,
            memoryGb = 4
    }
    
    call merge_vcf.IndexVcf as recalIndexVcf {
        input:
            vcfCompressed = recalCompressVcf.vcfCompressed
    }
    
    call germline.VcfNorm {
        input:
            sampleId = normal,
            haplotypecallerRecalVcf = recalIndexVcf.vcfCompressedIndexed    
    }
    
    call merge_vcf.IndexVcf as normIndexVcf {
        input:
            vcfCompressed = VcfNorm.haplotypecallerNormVcf
    }
    
    call germline.VarEval {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndel = indelApplyRecal.haplotypecallerGenoVcfApplySnpIndel,
            dbsnp = dbsnp   
    }
    
#    call germline.VarSum {
#        input:
#            sampleId = normal,
#            varReport = VarEval.varReport  
#    }
    
    output {
#        File varSummary = VarSum.varSummary
#        File varReport = VarEval.varReport
        IndexedVcf haplotypecallerNormVcf = normIndexVcf.vcfCompressedIndexed
    }
    
}
