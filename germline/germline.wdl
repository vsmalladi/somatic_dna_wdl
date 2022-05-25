version 1.0

import "../wdl_structs.wdl"


task haplotypeCallerGatk4 {
    input {
        Bam finalBam
        IndexedReference referenceFa
        String sampleId
        Int index
        String haplotypecallerIntervalVcfPath = "~{sampleId}.~{index}.haplotypecaller.g.vcf.gz"
        String bamOutPath = "~{sampleId}.bamout.bam"
        Int memoryGb = 7
        Int diskSize

        ## Inputs for haplotypecaller
        File excludeIntervalList
        File scatterIntervalsHc

    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)
    command {
            /gatk/gatk \
            --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R ~{referenceFa.fasta} \
            -I ~{finalBam.bam} \
            -L ~{scatterIntervalsHc} \
            -XL ~{excludeIntervalList} \
            -O ~{haplotypecallerIntervalVcfPath} \
            -G StandardAnnotation \
            -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -ERC GVCF \
            -bamout ~{bamOutPath}
    }

    output {
        IndexedVcf haplotypecallerIntervalVcf = object {
                vcf : "~{haplotypecallerIntervalVcfPath}",
                index : "~{haplotypecallerIntervalVcfPath}.tbi"
            }
        File bamOut = "~{bamOutPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

task GentotypeGvcfsGatk4 {
    input {
        IndexedReference referenceFa
        IndexedVcf sortedVcf
        String sampleId
        String index
        String haplotypecallerGenoVcfPath = "~{sampleId}.~{index}.haplotypecaller.gatk.genotypedGVCFs.vcf.gz"
        String haplotypecallerFilteredGenoVcfPath = "~{sampleId}.~{index}.haplotypecaller.gatk.filtered.genotypedGVCFs.vcf.gz"
        File scatterIntervalsHc

        ## Inputs for haplotypecaller
        IndexedVcf hapmap
        IndexedVcf omni
        IndexedVcf onekG
        IndexedVcf dbsnp

        Int memoryGb = 16
        Int diskSize = (ceil( size(sortedVcf.vcf, "GB") )  * 2 ) + 20
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        ## GenotypeGvcfs
        /gatk/gatk \
        --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        GenotypeGVCFs \
        -V ~{sortedVcf.vcf} \
        -R ~{referenceFa.fasta} \
        -L ~{scatterIntervalsHc} \
        -O ~{haplotypecallerGenoVcfPath}

        ## Score variants using CNN
        /gatk/gatk \
        --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        CNNScoreVariants \
        -V ~{haplotypecallerGenoVcfPath} \
        -R ~{referenceFa.fasta} \
        -L ~{scatterIntervalsHc} \
        -O ~{sampleId}.haplotypecaller.annotated.vcf

        ## Filter variant tranches
        /gatk/gatk \
        --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        FilterVariantTranches \
        -L ~{scatterIntervalsHc} \
        -V ~{sampleId}.haplotypecaller.annotated.vcf \
        -O ~{haplotypecallerFilteredGenoVcfPath} \
        --snp-tranche 99.9 --snp-tranche 99.95 \
        --indel-tranche 99.0 --indel-tranche 99.4 \
        --resource ~{hapmap.vcf} \
        --resource ~{omni.vcf} \
        --resource ~{onekG.vcf} \
        --resource ~{dbsnp.vcf} \
        --info-key CNN_1D \
        --create-output-variant-index true

    }

    output {
        IndexedVcf haplotypecallerGenoVcf = object {
                vcf : "~{haplotypecallerGenoVcfPath}",
                index : "~{haplotypecallerGenoVcfPath}.tbi"
            }
        IndexedVcf haplotypecallerFilteredGenoVcf = object {
                vcf : "~{haplotypecallerFilteredGenoVcfPath}",
                index : "~{haplotypecallerFilteredGenoVcfPath}.tbi"
            }
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}


task genotypeRefinementWorkflow {
    input {
        String sampleId
        IndexedVcf genotypedGatk4
        String haplotypecallerAfVcfPath = "~{sampleId}.haplotypecaller.gatk.af.vcf.gz"
        String haplotypecallerAfGqFilteredVcfPath = "~{sampleId}.haplotypecaller.gatk.af-gq-filtered.vcf.gz"
        IndexedReference referenceFa
        Int memoryGb = 16
        Int diskSize = (ceil( size(genotypedGatk4.vcf, "GB") ) * 4 ) + 20
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        set -e -o pipefail

        ## Annotate FORMAT/AF (deliver)
        /gatk/gatk \
        --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        VariantAnnotator \
        -R ~{referenceFa.fasta} \
        -V ~{genotypedGatk4.vcf} \
        -O ~{haplotypecallerAfVcfPath} \
        -A AlleleFraction


        ## remove biallellic sites
        zcat ~{haplotypecallerAfVcfPath} \
        | awk '($5 !~ ",")' \
        > ~{sampleId}.biallellic.vcf

        ## Variant filtration
        /gatk/gatk \
        --java-options "-Xmx~{jvmHeap}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        VariantFiltration \
        -R ~{referenceFa.fasta} \
        -V ~{sampleId}.biallellic.vcf \
        -O ~{sampleId}.haplotypecaller.af-gq-filtered.vcf.gz \
        --genotype-filter-name "AlleleFraction" \
        --genotype-filter-expression "(AF < 0.25 && AF > 0.0) || AF > 0.75" \
        --genotype-filter-name "GQ20" \
        --genotype-filter-expression "GQ < 20"

        # filter with AF (deliver)
        zcat ~{sampleId}.haplotypecaller.af-gq-filtered.vcf.gz \
        | grep -v "AlleleFraction" \
        > ~{haplotypecallerAfGqFilteredVcfPath}
    }

    output {
        IndexedVcf haplotypecallerAfVcf = object {
                vcf : "~{haplotypecallerAfVcfPath}",
                index : "~{haplotypecallerAfVcfPath}.tbi"
            }
        File haplotypecallerAfGqFilteredVcf = "~{haplotypecallerAfGqFilteredVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
    }
}

task filterHO {
        input {
            String sampleId
            IndexedVcf haplotypecallerAfVcf
            IndexedVcf nygcAf
            IndexedVcf pgx
            IndexedTable rwgsPgxBed
            IndexedVcf whitelist
            IndexedVcf chdWhitelistVcf
            IndexedVcf deepIntronicsVcf
            IndexedVcf clinvarIntronicsVcf
            String haplotypecallerFinalFilteredPath = "~{sampleId}.haplotypecaller.gatk.final.filtered.vcf.gz"
            Int diskSize = 100
            Int threads = 2
            Int memoryGb = 8
    }

    command <<<
          ## Remove existing AF annotations from merged VCF
          bcftools annotate \
          -x INFO/AF \
          -Oz \
          ~{haplotypecallerAfVcf.vcf} \
          > noaf.vcf.gz

          tabix -p vcf noaf.vcf.gz

          ## Annotate with NYGC AF for filtering
          bcftools annotate \
          --annotations ~{nygcAf.vcf} \
          --columns 'INFO/AF,INFO/AC_Hom' \
          -Oz \
          noaf.vcf.gz \
          > ~{sampleId}.final.annotated.vcf.gz

          tabix -p vcf ~{sampleId}.final.annotated.vcf.gz

          ## filter variants >3% AF and >10 Homozygotes in NYGC vars
          bcftools filter \
          --exclude 'INFO/AF[*] > 0.03 || INFO/AC_Hom[*] > 10' \
          ~{sampleId}.final.annotated.vcf.gz \
          > ~{sampleId}.pop.filtered.vcf

          bgzip ~{sampleId}.pop.filtered.vcf
          tabix -p vcf ~{sampleId}.pop.filtered.vcf.gz

          ## select whitelist variants
          bcftools view \
          -Oz \
          -R ~{whitelist.vcf} \
          ~{haplotypecallerAfVcf.vcf} \
          > ~{sampleId}.whitelist.filtered.vcf.gz

          tabix -p vcf ~{sampleId}.whitelist.filtered.vcf.gz

          ## select pgx variants
          bcftools view \
          -Oz \
          -R ~{pgx.vcf} \
          ~{haplotypecallerAfVcf.vcf} \
          > ~{sampleId}.pgx.filtered.vcf.gz

          tabix -p vcf ~{sampleId}.pgx.filtered.vcf.gz

          ## select chd whitelist variants
          bcftools view \
          -Oz \
          -R ~{chdWhitelistVcf.vcf} \
          ~{haplotypecallerAfVcf.vcf} \
          > ~{sampleId}.chdwhitelist.filtered.vcf.gz

          tabix -p vcf ~{sampleId}.chdwhitelist.filtered.vcf.gz

          ## select rwgs pgx variants
          bcftools view \
          -Oz \
          -R ~{rwgsPgxBed.table} \
          ~{haplotypecallerAfVcf.vcf} \
          > ~{sampleId}.rwgspgx.filtered.vcf.gz

          tabix -p vcf ~{sampleId}.rwgspgx.filtered.vcf.gz

          ## Select deep intronics
          bcftools view \
          -Oz \
          -R ~{deepIntronicsVcf.vcf} \
          ~{haplotypecallerAfVcf.vcf} \
          > ~{sampleId}.deep_intronics.filtered.vcf.gz

          tabix -p vcf ~{sampleId}.deep_intronics.filtered.vcf.gz

          ## Select clinvar intronics
          bcftools view \
          -Oz \
          -R ~{clinvarIntronicsVcf.vcf} \
          ~{haplotypecallerAfVcf.vcf} \
          > ~{sampleId}.clinvar_intronics.filtered.vcf.gz

          tabix -p vcf ~{sampleId}.clinvar_intronics.filtered.vcf.gz

          echo ~{sampleId} > samples.txt

          ## merge all filtered files for further processing
          bcftools concat \
          -a \
          -d all \
          ~{sampleId}.pop.filtered.vcf.gz \
          ~{sampleId}.whitelist.filtered.vcf.gz \
          ~{sampleId}.pgx.filtered.vcf.gz \
          ~{sampleId}.chdwhitelist.filtered.vcf.gz \
          ~{sampleId}.rwgspgx.filtered.vcf.gz \
          ~{sampleId}.deep_intronics.filtered.vcf.gz \
          ~{sampleId}.clinvar_intronics.filtered.vcf.gz \
          | \
          bcftools view \
          -i 'GT[@samples.txt]="alt"' \
          | \
          bcftools sort \
          -Oz \
          > ~{haplotypecallerFinalFilteredPath}

          tabix -p vcf ~{haplotypecallerFinalFilteredPath}
    >>>

    runtime {
      mem: memoryGb + "G"
      cpus: threads
      memory: memoryGb + " GB"
      cpu: threads
      disks: "local-disk " + diskSize + " HDD"
      docker: "gcr.io/nygc-public/genome-utils@sha256:59603ab0aeda571c38811d6d1820d1b546a69fc342120bef75597bfd7905ea1f"
    }

    output {
        IndexedVcf haplotypecallerFinalFiltered = object {
                vcf : "~{haplotypecallerFinalFilteredPath}",
                index : "~{haplotypecallerFinalFilteredPath}.tbi"
        }
    }
}
