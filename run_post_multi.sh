#!/bin/bash
# USAGE: run_post.sh -u URL -p PROJECT_ID -d LOG_DIR -r RUNINFO_JSON
# DESCRIPTION: review workflow from cromwell.
# Script requires jq, cromwell-tools, gcloud to be in the path
# Script generates summary of any workflow.
# If no *.RunInfo.json file is provided then
# the most recent file in the log dir will be used

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

init_wdl_dir=$(dirname "$0")
cd ${init_wdl_dir}
wdl_dir=$( pwd )

print_help() {
  echo "USAGE: run_post_multi.sh -u URL -p PROJECT_ID -d LOG_DIR -r RUNINFO_CSV" >&2
  echo "DESCRIPTION: submit workflow to cromwell." >&2
  echo "Script requires jq, cromwell-tools, gcloud to be in the path." >&2
  echo "Script also requires a default_credentials_JSON to be created by the user" >&2
  echo "Run the following to generate a default credentials file " >&2
  echo "$ gcloud auth application-default login" >&2
  echo "# Script generates summary of any workflow."
  exit 1
}

print_usage() {
  echo "USAGE: run_post.sh -u URL -p PROJECT_ID -d log_dir -r RUNINFO_JSON " >&2
  exit 1
}

while getopts 'u:p:d:r:h' flag; do
  case "${flag}" in
    u) url="${OPTARG}" ;;
    p) project_id="${OPTARG}" ;;
    r) run_info_json="${OPTARG}" ;;
    d) log_dir="${OPTARG}" ;;
    h) print_help ;;
    \?) print_usage; echo "Unknown option: $OPTARG" >&2 ;;
    *) print_usage; echo "Unimplemented option: $OPTARG" >&2 ;;
  esac
done



if [ -z "$log_dir" ]; then
    echo "Error: Missing required value for -l log dir (for run logs)" >&2
    print_usage
    exit 1
fi

if [ -z "$url" ]; then
    echo "Error: Missing required value for -u cromwell server URL" >&2
    print_usage
    exit 1
fi

if [ -z "$project_id" ]; then
    echo "Error: Missing required value for -p project_id" >&2
    print_usage
    exit 1
fi



set -e 
set -o pipefail

cd ${log_dir}

output_info_manifest="${project_id}_outputInfoManifest.txt"
metrics_info_manifest="${project_id}_MetricsInfoManifest.csv"

plot_file="${log_dir}/${project_id}_outputMetrics.html"
nav="${wdl_dir}/pandoc/nav_wgs_v7"
pandoc_dir="${wdl_dir}/pandoc/"
md="${log_dir}/${project_id}_outputMetrics.md"
html="${log_dir}/${project_id}_outputMetrics.html"
header="${log_dir}/${project_id}_outputMetrics.header.txt"


echo "Collect subworkflow uuids..."
time python ${wdl_dir}/tools/collect.py \
--multi \
--name ${project_id} \
--run-data ${run_info_json} \
--url ${url}


echo "Gather usage metrics..."
time python ${wdl_dir}/tools/runtime.py \
    --multi \
    --name ${project_id} \
    --manifest ${output_info_manifest}

echo "Plot usage metrics..." 
python ${wdl_dir}/tools/plot_runtime.py \
    --name ${project_id} \
    --manifest ${metrics_info_manifest} \
    --multi \
    --plot ${plot_file} 

    
time bash ${wdl_dir}/tools/html_printer.sh \
    ${md} \
    ${html} \
    ${header} \
    ${nav} \
    ${pandoc_dir}

