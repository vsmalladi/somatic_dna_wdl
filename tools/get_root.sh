#!/bin/bash
# USAGE: uuid=$( get_outputs.sh -u CROMWELL_URL -w WORKFLOW [-h])
# DESCRIPTION: get outputs objects from a cromwell workflow.
# Script requires jq, cromwell-tools to be in the path
# Script returns the workflow outputs.


print_help() {
  echo "USAGE: uuid=\$( get_outputs.sh -u URL -w WORKFLOW [-h])"
  echo "DESCRIPTION: get outputs objects from a cromwell workflow."
  echo "Script requires jq, cromwell-tools to be in the path."
  echo "Script returns the workflow outputs."
  exit 0
}

print_usage() {
  echo "USAGE: uuid=\$( run.sh -u URL -w WORKFLOW -o OPTIONS_JSON -i INPUTS )" >&2
  exit 1
}


while getopts 'u:w:g:h' flag; do
  case "${flag}" in
    u) url="${OPTARG}" ;;
    w) uuid="${OPTARG}" ;;
    g) gcp_project="${OPTARG}" ;;
    h) print_help ;;
    \?) print_usage; echo "Unknown option: $OPTARG" >&2 ;;
    :) print_usage; echo "Missing option argument for option: $OPTARG" >&2 ;;
    *) print_usage; echo "Unimplemented option: $OPTARG" >&2 ;;
  esac
done

if [ -z "$url" ]; then
    echo "Missing required value for -u URL for cromwell" >&2
    print_usage
    exit 1
fi

if [ -z "$uuid" ]; then
    echo "Missing required value for -w WORKFLOW UUID for cromwell" >&2
    print_usage
    exit 1
fi

set -e
set -o pipefail

if [ -z "$gcp_project" ]; then
    cromwell-tools metadata --url ${url} \
    --username $(gcloud secrets versions access latest --secret="cromwell_username") \
    --password $(gcloud secrets versions access latest --secret="cromwell_password") \
    --uuid ${uuid} \
    | jq ".workflowRoot"
else
    cromwell-tools metadata --url ${url} \
    --username $(gcloud secrets versions access latest --project=${gcp_project} --secret="readonly_cromwell_username") \
    --password $(gcloud secrets versions access latest --project=${gcp_project} --secret="readonly_cromwell_password") \
    --uuid ${uuid} \
    | jq ".workflowRoot"
fi
