version 1.0

import "../wdl_structs.wdl"
import "../germline/germline.wdl"
import "../calling/calling.wdl"
import "../merge_vcf/merge_vcf.wdl"

workflow Germline {
    # command 
    input {
        Bam FinalBam
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        String normal
        
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedReference referenceFa
        IndexedVcf onekG
        IndexedVcf dbsnp

        Int threads
        Int memoryGb
        String gatkDockerImage
        String gatk3_5DockerImage
        String bgzipDockerImage
        String bcftoolsDockerImage
    }
    scatter(chrom in listOfChroms) {
        call germline.Haplotypecaller {
            input:
                referenceFa = referenceFa,
                finalBam = FinalBam,
                chrom = chrom,
                sampleId = normal,
                memoryGb = memoryGb,
                threads = threads,
                dockerImage = gatk3_5DockerImage
        }
    }
    
    call calling.Gatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{normal}.haplotypeCalls.er.raw.vcf",
            tempChromVcfs = Haplotypecaller.haplotypecallerChromVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call merge_vcf.CompressVcf as haplotypecallerCompressVcf {
        input:
            vcf = Gatk4MergeSortVcf.sortedVcf.vcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bgzipDockerImage
    }
    
    call merge_vcf.IndexVcf as haplotypecallerIndexVcf {
        input:
            vcfCompressed = haplotypecallerCompressVcf.vcfCompressed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call germline.GentotypeGvcfs {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            sortedVcf = haplotypecallerIndexVcf.vcfCompressedIndexed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    call germline.RecalVcfsSnp {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            omni = omni,
            hapmap = hapmap,
            onekG = onekG,
            dbsnp = dbsnp,
            haplotypecallerGenoVcf = GentotypeGvcfs.haplotypecallerGenoVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    call germline.RecalVcfsIndel {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcf = GentotypeGvcfs.haplotypecallerGenoVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    call germline.ApplyRecal as snpApplyRecal {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndelPath = "~{normal}.SNP.apply_vqsr.vcf",
            haplotypecallerGenoVcfApply = GentotypeGvcfs.haplotypecallerGenoVcf,
            mode = "SNP",
            tranches = RecalVcfsSnp.tranchesSnp,
            recal = RecalVcfsSnp.recalSnp,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    call germline.ApplyRecal as indelApplyRecal {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndelPath = "~{normal}.SNP.INDEL.apply_vqsr.vcf",
            haplotypecallerGenoVcfApply = snpApplyRecal.haplotypecallerGenoVcfApplySnpIndel,
            mode = "INDEL",
            tranches = RecalVcfsIndel.tranchesIndel,
            recal = RecalVcfsIndel.recalIndel,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    call germline.VarFilter {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndel = indelApplyRecal.haplotypecallerGenoVcfApplySnpIndel,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    call merge_vcf.CompressVcf as recalCompressVcf {
        input:
            vcf = VarFilter.haplotypecallerRecalVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bgzipDockerImage
    }
    
    call merge_vcf.IndexVcf as recalIndexVcf {
        input:
            vcfCompressed = recalCompressVcf.vcfCompressed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call germline.VcfNorm {
        input:
            sampleId = normal,
            haplotypecallerRecalVcf = recalIndexVcf.vcfCompressedIndexed,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = bcftoolsDockerImage    
    }
    
    call merge_vcf.IndexVcf as normIndexVcf {
        input:
            vcfCompressed = VcfNorm.haplotypecallerNormVcf,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatkDockerImage
    }
    
    call germline.VarEval {
        input:
            referenceFa = referenceFa,
            sampleId = normal,
            haplotypecallerGenoVcfApplySnpIndel = indelApplyRecal.haplotypecallerGenoVcfApplySnpIndel,
            dbsnp = dbsnp,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    call germline.VarSum {
        input:
            sampleId = normal,
            varReport = VarEval.varReport,
            memoryGb = memoryGb,
            threads = threads,
            dockerImage = gatk3_5DockerImage    
    }
    
    output {
        File varSummary = VarSum.varSummary
        File varReport = VarEval.varReport
        IndexedVcf haplotypecallerNormVcf = normIndexVcf.vcfCompressedIndexed
    }
    
}
