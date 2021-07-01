version 1.0

import "../wdl_structs.wdl"
import "germline.wdl" as germline
import "../merge_vcf/merge_vcf.wdl" as mergeVcf
import "../calling/calling.wdl" as calling

workflow Germline {
    # command

    input {
        Bam finalBam
        IndexedReference referenceFa
        Array[String]+ listOfChroms
        String normal

        File excludeIntervalList
        Array[File] scatterIntervalsHcs

        IndexedVcf MillsAnd1000G
        IndexedVcf hapmap
        IndexedReference referenceFa
        IndexedVcf onekG
        IndexedVcf dbsnp

        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf

        Int hcDiskSize = ceil( size(finalBam.bam, "GB") ) + 20
    }
    scatter (i in range(length(scatterIntervalsHcs))) {
        call germline.haplotypeCallerGatk4 {
            input:
                referenceFa = referenceFa,
                finalBam = finalBam,
                index = i,
                sampleId=normal,
                hapmap=hapmap,
                onekG=onekG,
                dbsnp=dbsnp,
                excludeIntervalList=excludeIntervalList,
                scatterIntervalsHc=scatterIntervalsHcs[i],
                diskSize = hcDiskSize
        }
    }

    scatter (vcf in haplotypeCallerGatk4.haplotypecallerIntervalVcf) {
        File haplotypecallerIntervalVcf = vcf.vcf
    }
    Array[File] haplotypecallerIntervalVcfs = haplotypecallerIntervalVcf

    call calling.Gatk4MergeSortVcf as haplotypecallerGatk4MergeSortVcf {
        input:
            sortedVcfPath = "~{normal}.haplotypecaller.g.vcf",
            tempChromVcfs = haplotypecallerIntervalVcfs,
            referenceFa = referenceFa,
            memoryGb = 16,
            diskSize = 200
    }

    call mergeVcf.CompressVcf as haplotypecallerCompressVcf {
        input:
            vcf = haplotypecallerGatk4MergeSortVcf.sortedVcf.vcf
    }
    
    call mergeVcf.IndexVcf as haplotypecallerIndexVcf {
        input:
            vcfCompressed = haplotypecallerCompressVcf.vcfCompressed
    }

    scatter (i in range(length(scatterIntervalsHcs))) {
        call germline.GentotypeGvcfsGatk4 {
            input:
                referenceFa = referenceFa,
                sampleId=normal,
                index=i,
                sortedVcf = haplotypecallerIndexVcf.vcfCompressedIndexed,
                scatterIntervalsHc=scatterIntervalsHcs[i],
                hapmap = hapmap,
                onekG = onekG,
                dbsnp = dbsnp
        }
    }

    scatter (vcf in GentotypeGvcfsGatk4.haplotypecallerFilteredGenoVcf) {
        File haplotypecallerFilteredGenoVcf = vcf.vcf
    }
    Array[File] haplotypecallerFilteredGenoVcfs = haplotypecallerFilteredGenoVcf

    call calling.Gatk4MergeSortVcf as genotypedFilteredMergeSortVcf {
        input:
            sortedVcfPath = "~{normal}.haplotypecaller.gatk.v4.1.8.0.filtered.genotypedGVCFs.vcf",
            tempChromVcfs = haplotypecallerFilteredGenoVcfs,
            referenceFa = referenceFa,
            memoryGb = 16,
            diskSize = 200
    }

    call mergeVcf.CompressVcf as genotypedFilteredCompressVcf {
        input:
            vcf = genotypedFilteredMergeSortVcf.sortedVcf.vcf
    }
    
    call mergeVcf.IndexVcf as genotypedFilteredIndexVcf {
        input:
            vcfCompressed = genotypedFilteredCompressVcf.vcfCompressed
    }

    call germline.genotypeRefinementWorkflow {
        input:
            genotypedGatk4=genotypedFilteredIndexVcf.vcfCompressedIndexed,
            sampleId=normal,
            referenceFa=referenceFa
    }

    call germline.filterHO as filterHO {
        input:
            sampleId=normal,
            nygcAf = nygcAf,
            haplotypecallerAfVcf = genotypeRefinementWorkflow.haplotypecallerAfVcf,
            pgx=pgx,
            rwgsPgxBed=rwgsPgxBed,
            whitelist=whitelist,
            chdWhitelistVcf=chdWhitelistVcf,
            deepIntronicsVcf=deepIntronicsVcf,
            clinvarIntronicsVcf=clinvarIntronicsVcf
    }

    output {
        # equal to mergeHO.output_vcf
        IndexedVcf haplotypecallerVcf = genotypeRefinementWorkflow.haplotypecallerAfVcf
        # pass on to annotation
        IndexedVcf haplotypecallerFinalFiltered = filterHO.haplotypecallerFinalFiltered
    }
}
