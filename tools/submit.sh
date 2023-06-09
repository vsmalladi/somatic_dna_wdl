#!/bin/bash
# USAGE: uuid=$( submit.sh -u CROMWELL_URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS -d DEPENDENCIES -p PROJECT_DATA [-l LABLES] )
# DESCRIPTION: submit workflow to cromwell.
# Script requires jq, cromwell-tools, gcloud to be in the path
# Script returns the workflow uuid.
# Script shows submission command and command to check status .
# in the STDERR stream

print_help() {
  echo "USAGE: uuid=\$( submit.sh -u URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS -d DEPENDENCIES -p PROJECT_DATA [-l LABLES])"
  echo "DESCRIPTION: submit workflow to cromwell."
  echo "Script requires jq, cromwell-tools, gcloud to be in the path."
  echo "Script returns the workflow uuid."
  echo "Script shows submission command and command to check status"
  echo "in the STDERR stream."
  exit 1
}

print_usage() {
  echo "USAGE: uuid=\$( submit.sh -u URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS -d DEPENDENCIES -p PROJECT_DATA -f FILE_OUT_PREFIX [-l LABLES] )" >&2
  exit 1
}

while getopts 'u:w:o:d:p:f:i:l:h' flag; do
  case "${flag}" in
    u) url="${OPTARG}" ;;
    w) workflow="${OPTARG}" ;;
    o) options="${OPTARG}" ;;
    d) dependencies="${OPTARG}" ;;
    p) project_data="${OPTARG}" ;;
    f) file_out_prefix="${OPTARG}" ;;
    i) inputs="${OPTARG}" ;;
    l) labels="${OPTARG}" ;;
    h) print_help ;;
    \?) print_usage; echo "Unknown option: $OPTARG" >&2 ;;
    :) print_usage; echo "Missing option argument for option: $OPTARG" >&2 ;;
    *) print_usage; echo "Unimplemented option: $OPTARG" >&2 ;;
  esac
done

if [ -z "$project_data" ]; then
    echo "Missing required value for -p project_data json file" >&2
    print_usage
    exit 1
fi

set -e 
set -o pipefail

script_dir=$(dirname "$0")


# compose submission...
cmd="cromwell-tools submit --url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
-w ${workflow} -i ${inputs} -o ${options} -d ${dependencies}"

if [ ! -z "$labels" ]; then
    cmd="${cmd} \
    -l ${labels}"
fi

# submit submission...
# run workflow...
# echo $cmd >&2
response=$( eval $cmd )

# report status command and return UUID
uuid=$( echo ${response}  | jq -r ".id")
run_info_json="${file_out_prefix}.${uuid}.RunInfo.json"

# display submission status check command...
echo "# To print the current status of the run: " >&2
echo "cromwell-tools status --url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
--uuid ${uuid}" >&2

echo "# To print a detailed description of the run: " >&2
echo "cromwell-tools metadata --url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
--uuid ${uuid} \
| jq '.'" >&2

# return uuid
echo ${uuid}

# log uuid, project info, and inputs...
python ${script_dir}/log.py \
--project-data ${project_data} \
--uuid ${uuid} \
--inputs ${inputs} \
--file-out ${run_info_json}
