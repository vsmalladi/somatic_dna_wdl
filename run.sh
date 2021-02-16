#!/bin/bash
workflow=$1
options=$2
inputs=$3
shift
shift
shift
dependencies=$@

module load gcloud cromwell-tools
# run this before...
#gcloud config set project project-id
cmd="cromwell-tools submit --url https://cromwell-internal.nygenome.org \
--username $(gcloud secrets versions access latest --secret="cromwell_username") \
--password $(gcloud secrets versions access latest --secret="cromwell_password") \
-w ${workflow} -i ${inputs} -o ${options} -d ${dependencies}"


echo $cmd
eval $cmd
