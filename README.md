# wdl_port

## WDL docs for v6 pipeline

- [Dependencies](#dependencies)
- [On prem environment setup](#environment)
- [Write input](#write_input)
- [Run](#run)
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

Set up prem:
```
. /gpfs/commons/groups/nygcfaculty/kancero/envs/miniconda3/etc/profile.d/conda.sh
conda activate wdl
module load gcloud cromwell-tools jq
```

### Write input
<a name="write_input"/>

```
usage: meta.py [-h] [--library {WGS,Exome}] [--interval-list {WGS,Exome}]
               [--genome {Human_GRCh38_full_analysis_set_plus_decoy_hla}] 
               [--project PROJECT] [--pairs-file PAIRS_FILE]
               [--samples-file SAMPLES_FILE] [--wdl-file WDL_FILE] 
               [--project-data PROJECT_DATA]

optional arguments:
  -h, --help            show this help message and exit
  --library {WGS,Exome}
                        Sequence library type. If not supplied define library using 
                        --project-data
  --interval-list {SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla}
                        File basename for interval list.If not supplied the default
                        (the SureSelect interval list for your genome)
                        will be used
  --genome {Human_GRCh38_full_analysis_set_plus_decoy_hla}
                        Genome key to use for pipeline. If not supplied define 
                        genome using --project-data
  --project PROJECT     Project name associated with account. 
                        If not supplied define genome using --project-data
  --pairs-file PAIRS_FILE
                        JSON file with items that are required to have 
                        "tumor", "normal" sample_ids defined. If not supplied
                        define pairing using --project-data
  --samples-file SAMPLES_FILE
                        Not generally required. If steps run only require sample_id 
                        and do not use pairing information sample info
                        can be populated with a CSV file. The CSV file requires 
                        a columns named ["sampleId"]. If not supplied
                        define samples using --project-data
  --wdl-file WDL_FILE   WDL workflow. To output an input JSON that matches a 
                        WDL workflow parse the workflow file in as a flag.
  --project-data PROJECT_DATA
                        Optional JSON file with project pairing, sample, genome 
                        build, library and interval list information
```
Command
```
# Create input json
python wdl_port/tools/meta.py \
--project ${lab_quote_number} \
--pairs-file ${tumor_normal_pairs_csv} \
--library WGS \
--genome Human_GRCh38_full_analysis_set_plus_decoy_hla \
--wdl-file wdl_port/calling_wkf.wdl
```
Output:

  1. `calling_wkfInput.json` - inputs for cromwell
  2. `lab-number_projectInfo.json` - contains project info like the current list of samples/pairs and the library type as well as the pipeline version (tag and commit). Subsequent runs can use `--project-data  lab-number_projectInfo.json` and skip defining pair, library, interval list, genome, etc. 

### Run
<a name="run"/>

Running `run.sh` will return the UUID to STDOUT. 
It will also print the submitted command and the status command to the screen.
```
uuid=$( bash wdl_port/run.sh \
${url} \
wdl_port/calling_wkf.wdl \
options.json \
calling_wkfInput.json \
wdl_port/dependencies.zip )
```


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
--wdl-file wdl_port/new_wkf.wdl
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
```

6. Update the dependencies zip file

```
cd wdl_port
rm dependencies.zip
zip dependencies.zip wdl_structs.wdl */*.wdl
```

8. 

7. [Run](#run) the new workflow as before

