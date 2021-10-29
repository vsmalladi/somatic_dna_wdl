version 1.0

import "./wdl_structs.wdl"
import "calling/calling_wkf.wdl" as calling
import "calling/calling.wdl" as callingTasks
import "merge_vcf/merge_vcf_wkf.wdl" as mergeVcf
import "pre_process/qc.wdl" as qc
import "annotate/annotate_wkf.wdl" as annotate
import "annotate/annotate_cnv_sv_wkf.wdl" as annotate_cnv_sv

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


# for wdl version 1.0

task GetIndex {
    input {
        String sampleId
        Array[String] sampleIds
    }

    command {
        python /get_index.py \
        --sample-id ~{sampleId} \
        --sample-ids ~{sep=' ' sampleIds}
    }

    output {
        Int index = read_int(stdout())
    }

    runtime {
        docker: "gcr.io/nygc-internal-tools/workflow_utils:2.0"
    }
}


workflow SomaticBamWorkflow {
    input {
        BwaReference bwaReference
        IndexedReference referenceFa
        
        Boolean production = true
        Boolean external = false

        Array[pairInfo]+ pairInfos

        # calling
        Array[String]+ listOfChromsFull
        Array[String]+ listOfChroms
        IndexedTable callRegions
        File dbsnpIndels
        Map[String, File] chromBedsWgs
        File lancetJsonLog
        File mantaJsonLog
        File strelkaJsonLog
        File svabaJsonLog
        File mutectJsonLog
        File mutectJsonLogFilter
        File configureStrelkaSomaticWorkflow
        
        #   BicSeq2
        Int readLength
        Int coordReadLength
        Map[Int, Map[String, File]] uniqCoords
        File bicseq2ConfigFile
        File bicseq2SegConfigFile
        Map[String, File] chromFastas
        Int tumorMedianInsertSize = 400
        Int normalMedianInsertSize = 400
        Int lambda = 4
        
        # Gridss
        String bsGenome
        File ponTarGz
        Array[File] gridssAdditionalReference
        
        # merge callers
        File intervalListBed

        String library
        File ponWGSFile
        File ponExomeFile
        IndexedVcf gnomadBiallelic

        IndexedVcf germFile
        
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
        
        # annotate cnv
        File cytoBand
        File dgv
        File thousandG
        File cosmicUniqueBed
        File cancerCensusBed
        File ensemblUniqueBed
        
        # annotate sv
        String vepGenomeBuild
        # gap,DGV,1000G,PON,COSMIC
        File gap
        File dgvBedpe
        File thousandGVcf
        File svPon
        File cosmicBedPe
        
        # post annotation
        File cosmicCensus
        
        File ensemblEntrez
        
        IndexedVcf MillsAnd1000G
        IndexedVcf omni
        IndexedVcf hapmap
        IndexedVcf onekG
        IndexedVcf dbsnp
        
        IndexedVcf whitelist
        IndexedVcf nygcAf
        IndexedVcf pgx
        IndexedTable rwgsPgxBed
        IndexedVcf deepIntronicsVcf
        IndexedVcf clinvarIntronicsVcf
        IndexedVcf chdWhitelistVcf
        
        Boolean highMem = false
    }
    
    scatter(pairInfo in pairInfos) {
        
        # tumor insert size
        Int tumorDiskSize = ceil(size(pairInfo.tumorFinalBam.bam, "GB")) + 30
                      
        call qc.MultipleMetrics as tumorMultipleMetrics {
            input:
                referenceFa = referenceFa,
                finalBam = pairInfo.tumorFinalBam,
                sampleId = pairInfo.tumor,
                diskSize = tumorDiskSize
        }
        
        call callingTasks.GetInsertSize as tumorGetInsertSize {
            input:
                insertSizeMetrics = tumorMultipleMetrics.insertSizeMetrics
        }
        
        # normal insert size
        Int normalDiskSize = ceil(size(pairInfo.normalFinalBam.bam, "GB")) + 30
                      
        call qc.MultipleMetrics as normalMultipleMetrics {
            input:
                referenceFa = referenceFa,
                finalBam = pairInfo.normalFinalBam,
                sampleId = pairInfo.normal,
                diskSize = normalDiskSize
        }
        
        call callingTasks.GetInsertSize as normalGetInsertSize {
            input:
                insertSizeMetrics = normalMultipleMetrics.insertSizeMetrics
        }

        call calling.Calling {
            input:
                mantaJsonLog = mantaJsonLog,
                lancetJsonLog = lancetJsonLog,
                mutectJsonLog = mutectJsonLog,
                mutectJsonLogFilter = mutectJsonLogFilter,
                svabaJsonLog = svabaJsonLog,
                strelkaJsonLog = strelkaJsonLog,
                configureStrelkaSomaticWorkflow = configureStrelkaSomaticWorkflow,
                pairInfo = pairInfo,
                listOfChroms = listOfChroms,
                listOfChromsFull = listOfChromsFull,
                referenceFa = referenceFa,
                callRegions = callRegions,
                bwaReference = bwaReference,
                dbsnpIndels = dbsnpIndels,
                chromBedsWgs = chromBedsWgs,
                readLength = readLength,
                coordReadLength = coordReadLength,
                uniqCoords = uniqCoords,
                bicseq2ConfigFile = bicseq2ConfigFile,
                bicseq2SegConfigFile = bicseq2SegConfigFile,
                tumorMedianInsertSize = tumorGetInsertSize.insertSize,
                normalMedianInsertSize = normalGetInsertSize.insertSize,
                chromFastas = chromFastas,
                bsGenome = bsGenome,
                ponTarGz = ponTarGz,
                gridssAdditionalReference = gridssAdditionalReference,
                highMem = highMem
        }
        
        PreMergedPairVcfInfo preMergedPairVcfInfo = object {
            pairId : pairInfo.pairId,
            filteredMantaSV : Calling.filteredMantaSV,
            strelka2Snv : Calling.strelka2Snv,
            strelka2Indel : Calling.strelka2Indel,
            mutect2 : Calling.mutect2,
            lancet : Calling.lancet,
            svabaSv : Calling.svabaSv,
            svabaIndel : Calling.svabaIndel,
            tumor : pairInfo.tumor,
            normal : pairInfo.normal,
            tumorFinalBam : pairInfo.tumorFinalBam,
            normalFinalBam : pairInfo.normalFinalBam

        }

        PairRawVcfInfo pairRawVcfInfo = object {
            pairId : pairInfo.pairId,
            filteredMantaSV : Calling.filteredMantaSV,
            strelka2Snv : Calling.strelka2Snv,
            strelka2Indel : Calling.strelka2Indel,
            mutect2 : Calling.mutect2,
            lancet : Calling.lancet,
            svabaSv : Calling.svabaSv,
            svabaIndel : Calling.svabaIndel,
            gridssVcf : Calling.gridssVcf,
            bicseq2Png : Calling.bicseq2Png,
            bicseq2 : Calling.bicseq2,
            tumor : pairInfo.tumor,
            normal : pairInfo.normal,
            tumorFinalBam : pairInfo.tumorFinalBam,
            normalFinalBam : pairInfo.normalFinalBam

        }

        if (library == 'WGS') {
            call mergeVcf.MergeVcf as wgsMergeVcf {
                input:
                    external = external,
                    preMergedPairVcfInfo = preMergedPairVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponWGSFile,
                    germFile = germFile

            }
        }

        if (library == 'Exome') {
            call mergeVcf.MergeVcf as exomeMergeVcf {
                input:
                    external = external,
                    preMergedPairVcfInfo = preMergedPairVcfInfo,
                    referenceFa = referenceFa,
                    listOfChroms = listOfChroms,
                    intervalListBed = intervalListBed,
                    ponFile = ponExomeFile,
                    germFile = germFile

            }
        }

        File mergedVcf = select_first([wgsMergeVcf.mergedVcf, exomeMergeVcf.mergedVcf])
        
        call annotate.Annotate {
            input:
                unannotatedVcf = mergedVcf,
                referenceFa = referenceFa,
                production = production,
                tumor = pairInfo.tumor,
                normal = pairInfo.normal,
                pairName = pairInfo.pairId,
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
        
        call annotate_cnv_sv.AnnotateCnvSv {
            input:
                tumor = pairRawVcfInfo.tumor,
                normal = pairRawVcfInfo.normal,
                pairName = pairRawVcfInfo.pairId,
                listOfChroms =listOfChroms,
                bicseq2 = pairRawVcfInfo.bicseq2,
                cytoBand = cytoBand,
                dgv = dgv,
                thousandG = thousandG,
                cosmicUniqueBed = cosmicUniqueBed,
                cancerCensusBed = cancerCensusBed, 
                ensemblUniqueBed = ensemblUniqueBed,
                
                filteredMantaSV = pairRawVcfInfo.filteredMantaSV,
                svabaSv = pairRawVcfInfo.svabaSv,
                gridssVcf = pairRawVcfInfo.gridssVcf,
                vepGenomeBuild = vepGenomeBuild,
                gap = gap,
                dgvBedpe = dgvBedpe,
                thousandGVcf = thousandGVcf,
                svPon = svPon,
                cosmicBedPe = cosmicBedPe
        }
        
        FinalVcfPairInfo finalVcfPairInfo = object {
            pairId : pairInfo.pairId,
            tumor : pairInfo.tumor,
            normal : pairInfo.normal,
            mainVcf : Annotate.pairVcfInfo.mainVcf,
            supplementalVcf : Annotate.pairVcfInfo.supplementalVcf,
            vcfAnnotatedTxt : Annotate.pairVcfInfo.vcfAnnotatedTxt,
            maf : Annotate.pairVcfInfo.maf,
            filteredMantaSV : Calling.filteredMantaSV,
            strelka2Snv : Calling.strelka2Snv,
            strelka2Indel : Calling.strelka2Indel,
            mutect2 : Calling.mutect2,
            lancet : Calling.lancet,
            svabaSv : Calling.svabaSv,
            svabaIndel : Calling.svabaIndel,
            gridssVcf : Calling.gridssVcf,
            bicseq2Png : Calling.bicseq2Png,
            bicseq2 : Calling.bicseq2,
            cnvAnnotatedFinalBed : AnnotateCnvSv.cnvAnnotatedFinalBed,
            cnvAnnotatedSupplementalBed : AnnotateCnvSv.cnvAnnotatedSupplementalBed,
            svFinalBedPe : AnnotateCnvSv.svFinalBedPe,
            svHighConfidenceFinalBedPe : AnnotateCnvSv.svHighConfidenceFinalBedPe,
            svSupplementalBedPe : AnnotateCnvSv.svSupplementalBedPe,
            svHighConfidenceSupplementalBedPe : AnnotateCnvSv.svHighConfidenceSupplementalBedPe
        }
   }

    output {
        # CNV SV output and SNV INDELs
        Array[FinalVcfPairInfo] finalVcfPairInfos = finalVcfPairInfo

    }
}

