# NYGC Somatic Pipeline v7

## Docs for NYGC Somatic Pipeline v7
- [Deliverables](#deliverables)

![NYGC Somatic Pipeline overview](diagrams/pipeline_summary.png)

### Deliverables
<a name="deliverables"></a>

BAM alignments

- TUMOR.final.bam
- TUMOR.final.bai
- NORMAL.final.bam
- NORMAL.final.bai
 
SNV/INDEL/SV

- Raw SNV/INDEL/SV calls
  - TUMOR--NORMAL.indel.strelka2.v2.9.3.vcf
  - TUMOR--NORMAL.indel.svaba.v0.2.1.vcf
  - TUMOR--NORMAL.lancet.v1.0.7.vcf
  - TUMOR--NORMAL.manta.v1.4.0.filtered.vcf
  - TUMOR--NORMAL.mutect2.v4.0.5.1.vcf
  - TUMOR--NORMAL.snv.strelka2.v2.9.3.vcf
  - TUMOR--NORMAL.sv.gridss.v2.10.2.vcf.bgz
  - TUMOR--NORMAL.sv.gridss.v2.10.2.vcf.bgz.tbi
  - TUMOR--NORMAL.sv.svaba.v0.2.1.vcf
 
Merged SNV/INDEL

- SNV/INDEL calls annotated and filtered by PON and common germline variants
- HighConfidence calls
  - TUMOR--NORMAL.snv.indel.high_confidence.v7.annotated.vcf
- HighConfidence and AllSomatic (passing) calls
  - TUMOR--NORMAL.snv.indel.final.v7.annotated.vcf
- HighConfidence and AllSomatic (passing) calls (with all INFO fields from original callers)
  - TUMOR--NORMAL.snv.indel.supplemental.v7.annotated.vcf
- MAF formatted calls
  - TUMOR--NORMAL.snv.indel.final.v7.annotated.maf
- TABLE formatted calls
  - TUMOR--NORMAL.snv.indel.final.v7.annotated.txt

Merged SV

- SV calls annotated and filtered by PON and common germline variants
- HighConfidence calls
  - TUMOR--NORMAL.sv.annotated.v7.somatic.high_confidence.final.bedpe
- HighConfidence calls (with all INFO fields from original callers)
  - TUMOR--NORMAL.sv.annotated.v7.somatic.high_confidence.supplemental.bedpe
- HighConfidence and AllSomatic (passing) calls
  - TUMOR--NORMAL.sv.annotated.v7.somatic.final.bedpe
- HighConfidence and AllSomatic (passing) calls (with all INFO fields from original callers)
  - TUMOR--NORMAL.sv.annotated.v7.somatic.supplemental.bedpe

CNV

- Raw Bicseq2
  - TUMOR--NORMAL.bicseq2.v0.2.6.png
  - TUMOR--NORMAL.bicseq2.v0.2.6.txt
- CNV calls annotated and filtered by PON and common germline variants
  - TUMOR--NORMAL.cnv.annotated.v7.final.bed
  - TUMOR--NORMAL.cnv.annotated.v7.supplemental.bed

Germline

- All calls (including common germline variants)
  - NORMAL.haplotypecaller.gatk.v4.1.8.0.annotated.unfiltered.vcf
- Rare germline variants
  - NORMAL.haplotypecaller.gatk.v4.1.8.0.final.filtered.vcf.gz
  - NORMAL.haplotypecaller.gatk.v4.1.8.0.final.filtered.vcf.gz.tbi
- B-allele frequency
  - TUMOR--NORMAL.haplotypecaller.gatk.v4.1.8.0.baf.txt

MSI status

  - TUMOR--NORMAL.mantis.v1.0.4.WGS-targeted.status.final.tsv
 
HLA

- Kourami results
  - NORMAL.result
 
Mutational Signatures

- Raw deconstructSigs from HighConfidence variables
  - TUMOR--NORMAL.cosmic.v3.2.deconstructSigs.v1.9.0.signatures.highconfidence.counts.txt
  - TUMOR--NORMAL.cosmic.v3.2.deconstructSigs.v1.9.0.signatures.highconfidence.diff.txt
  - TUMOR--NORMAL.cosmic.v3.2.deconstructSigs.v1.9.0.signatures.highconfidence.input.txt
  - TUMOR--NORMAL.cosmic.v3.2.deconstructSigs.v1.9.0.signatures.highconfidence.reconstructed.txt
- Estimated proportion of each signature in the sample
  - TUMOR--NORMAL.cosmic.v3.2.deconstructSigs.v1.9.0.signatures.highconfidence.txt

