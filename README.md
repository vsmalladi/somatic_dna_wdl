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
Output is `calling_wkfInput.json`

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

