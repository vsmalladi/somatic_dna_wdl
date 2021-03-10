# wdl_port

## WDL docs for v6 pipeline

- [Dependencies](#dependencies)
- [On prem environment setup](#environment)
- [Write input](#write_input)
- [Run](#run)
- [Post run](#post_run)
- [Create new workflow](#create_new_workflow) 

![NYGC Somatic Pipeline overview](diagrams/WDL_Pipeline.png)

### Dependencies
<a name="dependencies"/>

- gcloud 
- gsutil
- cromwell-tools
- jq
- womtools

Python

- pandas
- jsonschema

### On prem environment setup
<a name="environment"/>

Set up on prem:

```
# run in this order
module load gcloud cromwell-tools jq
. /gpfs/commons/groups/nygcfaculty/kancero/envs/miniconda3/etc/profile.d/conda.sh
conda activate wdl
```

### Write input
<a name="write_input"/>

```
usage: meta.py [-h] --options OPTIONS --wdl-file WDL_FILE
               [--library {WGS,Exome}]
               [--genome {Human_GRCh38_full_analysis_set_plus_decoy_hla}]
               [--project PROJECT] [--pairs-file PAIRS_FILE]
               [--samples-file SAMPLES_FILE]
               [--interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}]
               [--project-data PROJECT_DATA]
               [--custom-inputs [CUSTOM_INPUTS [CUSTOM_INPUTS ...]]]
               [--skip-validate]

optional arguments:
  -h, --help            show this help message and exit
  --options OPTIONS     Options json file (required)
  --wdl-file WDL_FILE   WDL workflow. To output an input JSON that matches a
                        WDL workflow parse the workflow file in as a flag.
  --library {WGS,Exome}
                        Sequence library type. If not supplied define library
                        using --project-data
  --genome {Human_GRCh38_full_analysis_set_plus_decoy_hla}
                        Genome key to use for pipeline. If not supplied define
                        genome using --project-data
  --project PROJECT     Project name associated with account. If not supplied
                        define genome using --project-data
  --pairs-file PAIRS_FILE
                        JSON file with items that are required to have
                        "tumor", "normal" sample_ids defined. If not supplied
                        define pairing using --project-data
  --samples-file SAMPLES_FILE
                        Not generally required. If steps run only require
                        sample_id and do not use pairing information sample
                        info can be populated with a CSV file. The CSV file
                        requires a columns named ["sampleId"]. If not supplied
                        define samples using --project-data
  --interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}
                        File basename for interval list.If not supplied the
                        default (the SureSelect interval list for your genome)
                        will be used
  --project-data PROJECT_DATA
                        Optional JSON file with project pairing, sample,
                        genome build, library and interval list information
  --custom-inputs [CUSTOM_INPUTS [CUSTOM_INPUTS ...]]
                        Optional JSON file with custom input variables. The
                        name of the variable in the input file must match the
                        name of the variable in the WDL workflow. It is not
                        required that the input specify the workflow. By
                        default the input will be added to the top-level
                        workflow.
  --skip-validate       Skip the step where input files are validated.
                        Otherwise all gs//: URIs will be checked to see that a
                        file exists. Disable with caution.Cromwell will launch
                        instances and run without checking. Test a small pairs
                        file to ensure all references exist and at least some
                        sample input files can be read by the current user.

```

Command

```
# Create input json
python wdl_port/tools/meta.py \
--project ${lab_quote_number} \
--pairs-file ${tumor_normal_pairs_csv} \
--library WGS \
--genome Human_GRCh38_full_analysis_set_plus_decoy_hla \
--wdl-file wdl_port/calling_wkf.wdl \
--options options.json
```

#### Output:

  1. `calling_wkfInput.json` - inputs for cromwell
  2. `lab-number_projectInfo.json` - contains project info like the current list of samples/pairs and the library type as well as the pipeline version (tag and commit). Subsequent runs can use `--project-data  lab-number_projectInfo.json` and skip defining pair, library, interval list, genome, etc. 

# zip dependencies

This must be done because cromwell says so :). You must use zip.

```
cd wdl_port
zip dependencies.zip wdl_structs.wdl */*.wdl
cd -
```

### Run
<a name="run"/>

Running `run.sh` will return the workflow UUID to STDOUT. 
It will also print the submitted command and the status command to the screen.

```
uuid=$( bash ../wdl_port/run.sh \
    -u https://cromwell-compbio01.nygenome.org \
    -w ../wdl_port/calling_wkf.wdl \
    -o options.json \
    -d ../wdl_port/dependencies.zip \
    -p lab-number_projectInfo.json \
    -i calling_wkfInput.json )
```

In addition to submitting the command, this will create an output file that you should save with information about the project, pipeline version, cromwell options, inputs. It will also contain the workflow UUID.

#### Output:

  1. `lab-number_project.<DATE>.RunInfo.json` - contains projectInfo, workflow Input, workflow UUID.

### Post run
<a name="post_run"/>

Use the `cromwell-tools status` command printed to the screen when you submitted the workflow. Alternately, lookup the workflow UUID in the `lab-number_project.<DATE>.RunInfo.json` and run:

```
cromwell-tools status \
--url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
--uuid ${uuid}
```

If the workflow finishes and the status is `SUCCEEDED` then run `collect.py`. It will output `lab-number<WORKFLOW_UUID>_outputInfo.json`. This file contains all the information you need about a run. 
In addition to the content of `lab-number_project.<DATE>.RunInfo.json` the file includes:

  - `outputs`: map between workflow output object-name and object (including the file URIs)
  - `unnamed`: list of files in the output bucket that are not from this workflow
  - `pair_association` : map of pair_ids to the `outputs` map for just that pair 
  - `sample_association` : map of sample_ids to the `outputs` map for just that sample 
  
Note: Pair association only works if the pair_id is used in the filename followed by a `.` or a `_`. Sample association only works if the sample_id is used in the filename followed by a `.` or a `_`. 

```
python \
../wdl_port/tools/collect.py \
--run-data lab-number_project.<DATE>.RunInfo.json \
--url ${url}
```

#### Output:

  1. `lab-number<WORKFLOW_UUID>_outputInfo.json` - contains all runInfo, outputs map, pair_association map, sample_association map and unnamed files.


### Create new workflow
<a name="create_new_workflow"/>

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
gs://nygc-comp-s-fd4e-resources/GRCh38_full_analysis_set_plus_decoy_hla/internal/
# external reference files
gsutil cp \
${file} \
gs://nygc-comp-s-fd4e-resources/GRCh38_full_analysis_set_plus_decoy_hla/external/
```
4. Create an input JSON. This will create `new_wkfInput.json`

```
python wdl_port/tools/meta.py \
--project ${lab_quote_number} \
--pairs-file tumor_normal_pairs.csv \
--library WGS \
--genome Human_GRCh38_full_analysis_set_plus_decoy_hla \
--wdl-file wdl_port/new_wkf.wdl  \
--options options.json
```

5. Validate you workflow and inputs

```
# validate new files
womtool validate \
--inputs new_wkfInput.json \
new_wkf.wdl

# Also confirm that all workflows continue to be valid before commiting your changes
for file in */*wdl *wdl; do 
  echo $file; 
  womtool validate $file; 
done

# commit your changes and create new inputs file (with new branch and commit info)
```

6. Update the dependencies zip file

```
cd wdl_port
rm dependencies.zip
zip dependencies.zip wdl_structs.wdl */*.wdl
cd -
```

7. [Run](#run) the new workflow as before

