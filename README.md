# NYGC Somatic Pipeline v7



## WDL docs for NYGC Somatic Pipeline v7

- [Dependencies](#dependencies)
- [On prem environment setup](#environment)
- [Environment setup](#external_environment)
- [Available workflows](#workflows)
- [Write input and submit](#write_input)
- [Post run reports](#post_run)
- [Create new workflow](#create_new_workflow) 
- [Appendix](#appendix)
  - [Credentials](#credentials)
  - [Cram incompatible steps](#cram-incompatible)
  - [Cram compatible steps](#cram-compatible)
  - [Contact us](#contact_us)
  - [Release notes](#release_notes)

![NYGC Somatic Pipeline overview](diagrams/WDL_Pipeline.png)

### Dependencies
<a name="dependencies"></a>

- gcloud 
- gsutil
- cromwell-tools
- jq
- womtools

Python

- pandas
- jsonschema

Docker images
- [dockerfiles and images](https://bitbucket.nygenome.org/projects/DOC/repos/docker-images/browse)


### On prem environment setup
<a name="environment"></a>


Set up on prem:

```
module load wdl
```

### Environment setup
<a name="external_environment"></a>


Set up:

```
# Install gcloud cromwell-tools and jq and have them in you path

# Install conda (only need to do this install once)
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash \
Miniconda3-latest-Linux-x86_64.sh
# Then follow install instructions and then instructions to add conda to your path
# create wdl environment
conda env create -f wdl_port/tools/environment.yml

# activate the environment
conda activate wdl
```
gcloud cromwell-tools jq

### Available workflows
<a name="workflows"></a>

The pipeline is designed to be modular because there are times when we only run a segment 
and not the entire v7 pipeline. Below are available workflows:

```
alignment_analysis_wkf.wdl             Run Kourami (HLA typing) and Mantis (MSI status) on BAMs
annotate_cnv_sv_wkf.wdl                Merge SV calls, filter and annotation both SV and CNV calls
calling_wkf.wdl                        Run SNV, INDEL, SV and CNV callers on BAMs                                  
filter_intervals.wdl                   Prep reference files by filtering using a BED file of intervals
gdc_wkf.wdl                            Compare VCF files
germline_wkf.wdl                       Call germline SNVs and INDELs and output a filtered and unfiltered 
                                        annotated VCF
kourami_wkf.wdl                        Run Kourami (HLA typing)
merge_vcf_wkf.wdl                      Merge, filter and annotate v7 pipeline calls
report_mini_wkf.wdl                    In dev report writing workflow
report_wkf.wdl                         In dev full report writing workflow
somatic_bam_wkf.wdl                    Full v7 pipeline starting from BAMs: SNV, INDEL, SV and CNV calling 
                                        filtering and annotating, germline calls and BAF (HaplotypeCaller),
                                        DeconstructSigs (Mutational Signatures), contamination and concordance 
                                        (Conpair), Kourami (HLA typing) and 
                                        Mantis (MSI status)
somatic_wkf.wdl                        Full v7 pipeline starting from FASTQs: Alignment, QC, SNV, INDEL, SV 
                                        and CNV calling filtering and annotating, germline calls and BAF 
                                        (HaplotypeCaller), DeconstructSigs (Mutational Signatures), 
                                        contamination and concordance (Conpair), Kourami (HLA typing) and 
                                        Mantis (MSI status)
tests_wkf.wdl                          Automated comparison of prior pipeline run to current pipeline run output.
variant_analysis_wkf.wdl               Run DeconstructSigs (Mutational Signatures)
wdl_structs.wdl                        Custom Struct objects to reference primary and secondary files together 
                                        or group ids and related files
```

### Write input and submit
<a name="write_input"></a>


Script `wdl_port/run.sh` will first quickly validate the WDL workflow. Next it will use the WDL 
to determine which variables are required. All required variables will be defined from the 
Reference JSON in config, the pairing/sample info or the custom inputs JSON.

Most URIs for files will be validated to ensure you can read the file (unless the `--skip-validate` flag is used). If a pipeline has the variable `production` and/or `external` and these are set to true then the pipeline will skip tasks that require private NYGC files or licenses. Because of this if either `production` or `external` are true for a workflow private NYGC files will not have their uris validated.

Next, the input JSON will be compared to the original WDL and the program will exit if any required variable could not be found and will need to be provided in the custom input JSON.

Then the workflow will be submitted to the cromwell server.

In addition to submitting the command, this will create an output file that you should save. It contains information about the project, pipeline version, cromwell options, inputs. It will also contain the workflow UUID. This will be used after the run to agregate information about the run.
```
wdl_port/run.sh -h
run.sh [-h] --options OPTIONS --wdl-file WDL_FILE
               --url URL --log-dir LOG_DIR
               --project-name PROJECT_NAME
               [--library {WGS,Exome}]
               [--genome {Human_GRCh38_full_analysis_set_plus_decoy_hla, Human_GRCh38_tcga}]
               [--pairs-file PAIRS_FILE]
               [--samples-file SAMPLES_FILE]
               [--interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}]
               [--custom-inputs [CUSTOM_INPUTS [CUSTOM_INPUTS ...]]]
               [--skip-validate]
               [--dry-run]

DESCRIPTION: validate workflow, create input json and submit workflow to cromwell.
    Script requires jq, cromwell-tools, gcloud to be in the path.
    Script shows submission command and command to check status
    in the STDERR stream.
  Creation of input JSON:
    The WDL is used to determine which variables are required.
    Required or optional variables are defined from custom inputs JSON.
    Any required variable not defined in the custom inputs JSON will be defined from the
    reference JSONs in the config directory (as long as the variable names are identical).
    The pairing/sample info CSVs (--pairs-file/--samples-file) are used to create pairRelationships and
    (if columns named tumorBam and normalBam exist) map BAMs to pairs.
    If "production" or "external" inputs are true then validation of NYGC internal-only files is skipped
    The pipelines are written to skip tasks that localize these files if "production" or "external" are true
    so any inability to read these files will not negatively affect the run.
  -h, --help            show this help message and exit
  --url URL             Cromwell server URL (required)
  --log-dir LOG_DIR     Output directory for all logs and reports
                        related to this workflow UUID (required)
  --options OPTIONS     Options json file for cromwell (required)
  --wdl-file WDL_FILE   WDL workflow. An input JSON that matches this
                        WDL workflow will be created (required)
  --library {WGS,Exome}
                        Sequence library type.
  --genome {Human_GRCh38_full_analysis_set_plus_decoy_hla, Human_GRCh38_tcga}
                        Genome key to use for pipeline.
  --project-name PROJECT Project name associated with account.
  --pairs-file PAIRS_FILE
                        CSV file with items that are required to have
                        "tumor", "normal" and "pairId" in the columns.
                        Optionally, include "tumorBam", "normalBam" columns to create
                        "pairInfos" and "normalSampleBamInfos" automatically.
  --samples-file [SAMPLES_FILE]
                        Not generally required. If tasks only require
                        sampleId and do not use pairing information sample
                        info can be populated with a CSV file. The CSV file
                        requires a columns named ["sampleId"].
  --interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}
                        File basename for interval list.If not supplied the
                        default (the SureSelect interval list for your genome)
                        will be used (only needed for future Exome workflows)
  --custom-inputs [CUSTOM_INPUTS]
                        Optional JSON file with custom input variables. The
                        name of the variable in the input file must match the
                        name of the variable in the WDL workflow. It is not
                        required that the input specify the workflow. By
                        default the input will be added to the top-level
                        workflow. Any variable defined in in this JSON will
                        overwrite any reference variable in the the config
                        directory or workflow default.
  --skip-validate       Skip the step where input files are validated.
                        Otherwise all gs//: URIs will be checked to see that a
                        file exists. Disable with caution. Cromwell will launch
                        instances and run without checking. Test a small pairs
                        file to ensure all references exist and at least some
                        sample input files can be read by the current user.
  --dry-run             Skip the step where the job is submitted to cromwell-tools.
```

Command

```
# Create input json
cd ${working-dir}

../wdl_port/run.sh \
--log-dir ${working-dir} \
--url ${url} \
--project-name ${lab_quote_number} \
--pairs-file ${tumor_normal_pairs_csv} \
--library WGS \
--genome Human_GRCh38_full_analysis_set_plus_decoy_hla \
--wdl-file wdl_port/somatic_wkf.wdl \
--options options.json
```


#### Output:

  1. `somatic_wkfInput.json` - inputs for cromwell
  2. `${lab_quote_number}_projectInfo.json` - contains project info like the current list of samples/pairs and the library type as well as the pipeline version (tag and commit). 
  3. `${lab_quote_number}_project.<DATE>.RunInfo.json` - contains projectInfo, workflow Input, workflow UUID.


### Post run
<a name="post_run"></a>

Use the `cromwell-tools status` command printed to the screen when you submitted the workflow. Alternately, lookup the workflow UUID in the `lab-number_project.<DATE>.RunInfo.json` and run:

```
cromwell-tools status \
--url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
--uuid ${uuid}
```

After workflow finishes and the status is `SUCCEEDED` run:

```
cd ${working-dir}

bash ../wdl_port/run_summary.sh \
-u ${url} \
-d ${log_dir} \
-p ${gcp_project} \
-n ${project_name}

# optionally add -r ${run_info} or the script will use the most recent *RunInfo.json file in the ${log_dir}
```

If you do not have a `*RunInfo.json` you can start with the workflow uuid and optionally a sample or pair csv file

```
bash \
../wdl_port/run_summary.sh \
-u ${url} \
-d ${log_dir} \
-b ${billing_export} \
-n ${project_name} \
-g ${gcp_project} \
--uuid ${uuid} \
--samples-file ${sample_id_list}
```

NOTE: Optionally add (where `${billing_export}` is the name of the table for the SQL query) to add cost per instance_id to results:
```
-b ${billing_export}
```

#### Output:
If the workflow finishes and the status is `SUCCEEDED` and you have run `run_summary.sh`. It will output several log files:

`${lab_quote_number}<WORKFLOW_UUID>_outputInfo.json` :

In addition to the content of `${lab_quote_number}_project.<DATE>.RunInfo.json` the file includes:

  - `outputs`: map between workflow output object-name and object (including the file URIs)
  - `pair_association` : map of pair_ids to the `outputs` map for just that pair (searches for the pair_id followed by . _ or / )
  - `sample_association` : map of sample_ids to the `outputs` map for just that sample (searches for the sample_id followed by . _ or / )

Note: Pair association only works if the pair_id is used in the filename followed by a `.`, `/` or a `_`. Sample association only works if the sample_id is used in the filename followed by a `.` or a `_`. 


`${lab_quote_number}<WORKFLOW_UUID>_outputMetrics.csv`
the file includes the following metrics calculated from the BigQuery cromwell monitor:

    'id', 'project_id', 'zone', 'instance_name', 
    'preemptible', 'workflow_name', 'workflow_id', 
    'task_call_name', 'shard', 'attempt', 
    'start_time', 'end_time', 'execution_status', 
    'cpu_count', 'mem_total_gb', 'disk_mounts', 
    'disk_total_gb', 'disk_types', 'docker_image',
    'inputs', 'run_time', 'run_time_m', 
    'cpu_time_m', 'mean_task_core_h', 'mean_task_run_time_h', 
    'sample_task_run_time_h', 'max_mem_g', 'sample_task_core_h',
    'sample_subworkflow_core_h', 'sample_subworkflow_run_time_h', 
    'subworkflow_max_mem_g', 'sample_workflow_core_h', 
    'sample_workflow_run_time_h', 'workflow_max_mem_g'
    
`${lab_quote_number}<WORKFLOW_UUID>_outputMetrics.html`
This file includes plots of runtime metrics.
  
  

### Create new workflow
<a name="create_new_workflow"></a>

Use a [style guide](https://biowdl.github.io/styleGuidelines.html) to write you WDL files.

1. Write a new workflow (e.g. `wdl_port/new_wkf.wdl`) using structs from `wdl_structs.wdl` were needed. Use the variable from `config/fasta_references.json` and `config/interval_references.json` in your workflow (e.g. `referenceFa`)
  - keep your tasks and sub workflow in a separate WDL file in a subdirectory. That workflow should run a section of the pipeline on one sample/pair. 
  - In the main directory make a WDL that runs the sub workflow(s) for a list of `sampleInfo` or `pairInfo` objects.
  - Alternately add your subworkflow to an existing main workflow (e.g. `calling_wkf.wdl`)
2. Add any new required resource files to `config/fasta_references.json` and `config/interval_references.json`.
3. Upload any new resource files:

```
# modified or in-house files
gsutil cp \
${file} \
gs://${resources_project}/GRCh38_full_analysis_set_plus_decoy_hla/internal/
# external reference files
gsutil cp \
${file} \
gs://${resources_project}/GRCh38_full_analysis_set_plus_decoy_hla/external/
```
4. Confirm all workflows are still valid and commit changes

Also confirm that all workflows continue to be valid before commiting your changes

```
for file in */*wdl *wdl; do 
  echo $file; 
  womtool validate $file; 
done
# commit your changes and create new inputs file (with new branch and commit info)
git commit -m "feat: my new workflow"
```

5. [Run](#write_input) the new workflow as before

### Appendix
<a name="appendix"></a>

##### Credentials
<a name="credentials"></a>
Run the following once to generate a default credentials file

```
    $ gcloud auth application-default login
```

See for more details: https://google-auth.readthedocs.io/en/latest/reference/google.auth.html#google.auth.default.

##### Cram incompatible steps
<a name="cram-incompatible"></a>

- biqsec2 norm
- svaba?
- lancet
- mantis kmer counter ?


##### Cram compatible steps
<a name="cram-compatible"></a>

- strelka2
- mutect2
- gridss
- manta
- all samtools/pysam steps
- kourami w/ custom prep

### Contact us 
<a name="contact_us"></a>

We are in the process of setting up a public issue tracker. In the mean time  please send suggestions, questions or issues to jshelton@nygenome.org


# Release Notes
<a name="release_notes"></a>

[7.3.3] Refactor:

    - switch from tags to sha has for docker images
    - make bicseq2 config files reference files
    - begin pipeline output docs
    - switch to merge specific preMergedPairVcfInfo object for merge
    - add environment.yml
    - get rc file for more complete uuid list (for any not listed in subworkflow final output)
    - added SNV/INDEL only pipeline
    - added SNV/INDEL, SV and CNV only pipeline
    - added 'external' boolean to skip internal files
    - switch to public docker images
    - switch to public reference files
    - add input examples and example commands in docs
    - remove --read-length flag and replace with reading from input json
    - remove all private files from input JSON

[7.3.2](https://bitbucket.nygenome.org/rest/api/latest/projects/COMPBIO/repos/wdl_port/archive?at=refs%2Ftags%2F7.3.2&format=zip) Refactor:

    - adjust mem and disk size
    - finalize DeconstructSigs workflow
    - add highMem flag
    - calculate jvmHeap from mem
    - decrease mem and disk sizes
    - switch from local-disk SSD/HDD to local-disk LOCAL for steps that run on bams
    - add gridss arrange steps that works with cache
    - update resource usage scripts
    
[7.3.1](https://bitbucket.nygenome.org/rest/api/latest/projects/COMPBIO/repos/wdl_port/archive?at=refs%2Ftags%2F7.3.1&format=zip) Refactor:

    - add deconstructsigs
    - get chr6 coordinates from smaller file
    - adjust threads
    - speed up allele counts (chrom splits)
    

[7.2.0](https://bitbucket.nygenome.org/rest/api/latest/projects/COMPBIO/repos/wdl_port/archive?at=refs%2Ftags%2F7.2.0&format=zip) GDC references:

    - add Human_GRCh38_tcga
    - populate BAMs from a table
    - restrict to .bai index extensions
    - add workflow to reheader interval lists for external reference files

### Contact us 
<a name="contact_us"></a>

We are in the process of setting up a public issue tracker. In the mean time  please send suggestions, questions or issues to jshelton@nygenome.org
