#!/bin/bash
# USAGE: get_uuids.sh -u CROMWELL_URL -w WORKFLOW [-h]
# DESCRIPTION: get outputs objects from a cromwell workflow.
# Script requires jq, cromwell-tools to be in the path
# Script returns the workflow outputs.


print_help() {
  echo "USAGE: get_uuids.sh -u URL -w WORKFLOW [-h])"
  echo "DESCRIPTION: get outputs objects from a cromwell workflow."
  echo "Script requires jq, cromwell-tools to be in the path."
  echo "Script returns the workflow outputs."
  exit 0
}

print_usage() {
  echo "USAGE: get_uuids.sh -u URL -w WORKFLOW"
  exit 1
}


while getopts 'u:w:h' flag; do
  case "${flag}" in
    u) url="${OPTARG}" ;;
    w) uuid="${OPTARG}" ;;
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


username=$(gcloud secrets versions access latest --secret="cromwell_username")
password=$(gcloud secrets versions access latest --secret="cromwell_password")
curl -X GET "${url}/api/workflows/v1/${uuid}/metadata?expandSubWorkflows=true&includeKey=outputs" \
-u ${username}:${password} \
| jq .

