#!/bin/bash
# USAGE: uuid=$( run.sh -u CROMWELL_URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS... )
# DESCRIPTION: submit workflow to cromwell.
# Script requires jq, cromwell-tools, gcloud to be in the path
# Script returns the workflow uuid.
# Script shows submission command and command to check status .
# in the STDERR stream

print_help() {
  echo "USAGE: uuid=\$( run.sh -u URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS )"
  echo "DESCRIPTION: submit workflow to cromwell."
  echo "Script requires jq, cromwell-tools, gcloud to be in the path."
  echo "Script returns the workflow uuid."
  echo "Script shows submission command and command to check status"
  echo "in the STDERR stream."
  exit 1
}

print_usage() {
  echo "USAGE: uuid=\$( run.sh -u URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS )"
  exit 1
}

while getopts 'u:w:o:d:p:i:h' flag; do
  case "${flag}" in
    u) url="${OPTARG}" ;;
    w) workflow="${OPTARG}" ;;
    o) options="${OPTARG}" ;;
    d) dependencies="${OPTARG}" ;;
    p) project_data="${OPTARG}" ;;
    i) inputs="${OPTARG}" ;;
    h) print_help ;;
    \?) print_usage; echo "Unknown option: $OPTARG" >&2 ;;
    :) print_usage; echo "Missing option argument for option: $OPTARG" >&2 ;;
    *) print_usage; echo "Unimplemented option: $OPTARG" >&2 ;;
  esac
done

if [ -z "$project_data" ]; then
    echo "Missing required value for -p project_data json file"
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

# submit submission...
# run workflow...
echo $cmd >&2
response=$( eval $cmd )

# report status command and return UUID
uuid=$( echo ${response}  | jq -r ".id")

# display submission status check command...
echo "cromwell-tools status --url ${url} \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
--uuid ${uuid}" >&2

# return uuid
echo ${uuid}

# log uuid, project info, and inputs...
python ${script_dir}/tools/log.py \
--project-data ${project_data} \
--uuid ${uuid} \
--inputs ${inputs}
