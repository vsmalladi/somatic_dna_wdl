import sys
from google.cloud import bigquery
from google.cloud import bigquery_storage
from google_auth_oauthlib import flow
import google.auth
import json
import datetime
import logging as log
import pandas as pd
import argparse
import os 
import numpy as np
import subprocess
import pprint
import ast


try:
    import Colorer
except ImportError:
    pass

pd.set_option('display.max_columns', 500)
log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class Cost():
    def __init__(self,
                 output_metrics,
                 url,
                 out_file_prefix,
                 billing_export,
                 limit=1000000,
                 gcp_project=False,
                 test=False):
        self.test = test
        self.billing_export = billing_export
        self.parent_dir = os.path.abspath(os.path.dirname(__file__))
        self.output_metrics = output_metrics
        self.out_file_prefix = out_file_prefix
        self.url = url
        self.limit = str(limit)
        # login
        self.credentials, self.gcp_project = self.login()
        if gcp_project:
            self.default_project = False
            self.gcp_query_project = gcp_project
        else:
            self.default_project = True
            self.gcp_query_project = self.gcp_project
        # gather big cloud tables
        self.bqclient = bigquery.Client(project=self.gcp_project , 
                                        credentials=self.credentials)
        self.bqstorageclient = bigquery_storage.BigQueryReadClient(credentials=self.credentials)
        for file in self.output_metrics:
            self.process_metadata(file)
    
    def process_metadata(self, file):
        '''Get costs for one workflow uuid'''
        self.metadata = pd.read_csv(file)
        # hack fix later by adding main uuid to metrics file
        self.workflow_uuid = file.split('.')[1].replace('_outputMetrics', '')
        # seg run_date and end_date
        self.get_times()
        log.info('lookup costs')
        self.load_cost()
        print(self.cost.head())
        if not self.test:
            log.info('write metrics')
            cost_file = self.out_file_prefix + '_outputCosts.csv'
            self.cost.to_csv(cost_file, index=False, float_format='{:f}'.format, encoding='utf-8')
        
    def login(self):
        '''Run the following to generate a default credentials file
    
        $ gcloud auth application-default login
    
        https://google-auth.readthedocs.io/en/latest/reference/google.auth.html#google.auth.default.'''
        credentials, gcp_project = google.auth.default()
        return credentials, gcp_project

    def load_end_date(self, date_record):
        try:
            return date_record.strftime(format='%Y-%m-%d')
        except ValueError:
            return False     

    def modify_date(self, run_date):
        return '-'.join(run_date.split('-')[0:3])
 
    
    def get_times(self):
        '''find start and endtime for workflow '''
        self.metadata['start_time'] = pd.to_datetime(self.metadata['start_time'])
        self.metadata['end_time'] = pd.to_datetime(self.metadata['end_time'])
        self.metadata['padded_start_time'] = self.metadata['start_time'] - datetime.timedelta(days=1)
        self.metadata['padded_end_time'] = self.metadata['end_time'] + datetime.timedelta(days=1)
        self.run_date = self.load_end_date(self.metadata.padded_start_time.min())
        self.end_date = self.load_end_date(self.metadata.padded_end_time.max())

    
    def load_cost(self):
        '''convert big query runtime table into pandas dataframe
        Includes ? and workflow_id.
SELECT
    (SELECT value from UNNEST(labels) where key = 'cromwell-workflow-id' limit 1) as workflow_id,
    (SELECT value from UNNEST(labels) where key = 'cromwell-sub-workflow-name' limit 1) as sub_workflow_name,
    (SELECT value from UNNEST(labels) where key = 'wdl-task-name' limit 1) as wdl_task_name,
    sum(cost) as cost
FROM `fin-admin-328417.billing_export.gcp_billing_export_v1_01F75B_F18B56_39F532` 
WHERE DATE(_PARTITIONTIME) > "2021-10-01"
    AND project.id = "nygc-comp-p-12d3"
    AND sku.description != "Network Google Ingress from Americas to Americas"
    AND cost > 0
    AND EXISTS (SELECT 1 FROM UNNEST(labels) WHERE key = 'cromwell-workflow-id' and value in ("cromwell-954927a1-cd2c-4991-a771-38b3b36fa549", "cromwell-9c49910d-fa48-494c-971d-36247d4e919e"))
group by sub_workflow_name, wdl_task_name, workflow_id 
order by cost desc
LIMIT 1000
        '''
        log.info('Loading Cost: ' + self.workflow_uuid + '...')
        query_string = '''SELECT
    (SELECT value from UNNEST(labels) where key = 'cromwell-workflow-id' limit 1) as workflow_id,
    (SELECT value from UNNEST(labels) where key = 'cromwell-sub-workflow-name' limit 1) as sub_workflow_name,
    (SELECT value from UNNEST(labels) where key = 'wdl-task-name' limit 1) as wdl_task_name,
    sum(cost) as cost
    FROM `''' + self.billing_export + '''` 
    WHERE DATE(_PARTITIONTIME) >= "''' + self.run_date + '''"
        AND DATE(_PARTITIONTIME) <= "''' +self.end_date + '''"
        AND project.id = "''' + self.gcp_query_project + '''"
        AND sku.description != "Network Google Ingress from Americas to Americas"
        AND cost > 0
        AND EXISTS (SELECT 1 FROM UNNEST(labels) WHERE key = 'cromwell-workflow-id' and value in ("cromwell-''' + self.workflow_uuid + '''"))
    group by sub_workflow_name, wdl_task_name, workflow_id
    order by cost desc
    LIMIT ''' + self.limit + '''
        '''
        print(query_string)
        if self.test:
            print(query_string)
        else:
            self.cost = (
                self.bqclient.query(query_string)
                .result()
                .to_dataframe(bqstorage_client=self.bqstorageclient)
            )
        log.info('Done loading Cost: ' + self.workflow_uuid )


def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-file-prefix',
                        help='Prefix for cost summary file',
                        required=True
                        )
    parser.add_argument('--output-metrics',
                        help='CSV file with outputMetrics. Includes workflow uuid in the name. '
                        'Also includes task start_time and end_time',
                        nargs='*',
                        required=True
                        )
    parser.add_argument('--gcp-project',
                        help='GCP project id if this is not the set project '
                        'assumes you have read only access',
                        default=False,
                        required=False
                        )
    parser.add_argument('--url',
                        help='URL for cromwell server.',
                        required=True
                        )
    parser.add_argument('--limit',
                        help='Limit for bigquery SQL.',
                        default=10000,
                        required=False
                        )
    parser.add_argument('--dry-run',
                        help='Write but skip running the bigquery SQL query.',
                        action='store_true'
                        )
    parser.add_argument('--billing-export',
                        help='Billing export value for the bigquery SQL query.',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    cost = Cost(output_metrics=args['output_metrics'],
                url=args['url'],
                out_file_prefix=args['out_file_prefix'],
                limit=args['limit'],
                gcp_project=args['gcp_project'],
                billing_export=args['billing_export'],
                test=args['dry_run'])


if __name__ == "__main__":
    main()       
    
