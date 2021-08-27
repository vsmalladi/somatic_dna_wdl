import sys
from google.cloud import bigquery
from google.cloud import bigquery_storage
import json
import datetime
import logging as log
import pandas as pd
import Colorer
import make_auth
import argparse

pd.set_option('display.max_columns', 500)
log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class Runtime():
    def __init__(self,
                 output_info,
                 limit=10000,
                 metrics_file=False):
        self.output_info = output_info
        if not metrics_file:
            print( output_info['workflow_uuid'])
            print(output_info["project_data"]['project'])
            self.metrics_file = output_info["project_data"]['project'].replace(' ', '_') + '.' +  output_info['workflow_uuid'] + '_outputMetrics.csv'
        else:
            self.metrics_file = metrics_file
        # find relevant project info
        self.limit = str(limit)
        self.metrics_limit = str(limit * 4)
        # login
        self.credentials, self.gcp_project = make_auth.login()
        self.sub_workflow_uuids = self.load_sub_workflow_uuids()
        self.run_date = self.get_rundate()
        self.get_sample_info()
        # gather big cloud tables
        self.bqclient = bigquery.Client(project=self.gcp_project , 
                                        credentials=self.credentials)
        self.bqstorageclient = bigquery_storage.BigQueryReadClient(credentials=self.credentials)
        self.load_metadata()
        self.load_runtime()
        self.load_metrics()
        hashable = ['instance_name', 'workflow_id', 'attempt']
        self.metadata = self.metadata.drop_duplicates(subset=hashable)
        self.instance_id_map = self.get_max_mems(self.metadata, 
                                                 instance_id_map=self.instance_id_map, 
                                                 execution_status=False,
                                                 metrics_key='mem_used_gb')
        self.non_retry_instance_id_map = self.get_max_mems(self.metadata, 
                                                           instance_id_map=self.non_retry_instance_id_map, 
                                                           execution_status='Done',
                                                           metrics_key='mem_used_gb')
        self.non_retry_metadata = self.metadata[self.metadata.execution_status == 'Done'].copy()
        if self.metadata[self.metadata.execution_status != 'RUNNING'].empty:
            self.metadata.to_csv(self.metrics_file, index=False)
        else:
            self.metadata = self.load_plot_metrics(self.metadata,
                                                   instance_id_map=self.instance_id_map)
            self.non_retry_metadata = self.load_plot_metrics(self.non_retry_metadata,
                                                             instance_id_map=self.non_retry_instance_id_map)
            self.metadata.to_csv(self.metrics_file, index=False)
            self.non_retried_metrics_file = self.metrics_file.replace('.csv', '.non_retried.csv')
            self.non_retry_metadata.to_csv(self.non_retried_metrics_file, index=False)
        
    @staticmethod
    def read(file):
        with open(file) as output_info_file:
            output_info = json.load(output_info_file)
            return output_info
        
    def get_task_core_h(self, grouped, row, ids=['id', 
                                                 'task_call_name', 
                                                 'workflow_name']):
        ''' get runtime in core hours for the task (within a workflow) '''
        return (row.cpu_count * float(row.run_time_m)) / 60
    
    def match_instance_id(self, instance_id_map, row, execution_status='Done'):
        for inst in instance_id_map:
            if (row.workflow_id in instance_id_map[inst]['workflow_ids'] \
                    and row.task_call_name == instance_id_map[inst]['task_call_name']):
                if execution_status:
                    if execution_status == row.execution_status:
                        yield inst
                else:
                    yield inst
            
    
    def get_max_mems(self, metadata, instance_id_map, 
                     execution_status='Done',
                     job_id_col='task_call_name', metrics_key='mem_used_gb'):
        '''Get max mem for a runs with the same sub-workflow_id and task_call_name'''
        runtime_results = metadata.drop_duplicates(['task_call_name', 'workflow_id', 'instance_name'])
        max_values = {}
        for index, row in runtime_results.iterrows():
            instance_ids = []
            for inst in self.match_instance_id(instance_id_map, 
                                               row, 
                                               execution_status=execution_status):
                instance_ids.append(inst)
            instance_ids = list(set(instance_ids))
            try:
                max_value = self.metrics[self.metrics.instance_id.isin(instance_ids)][metrics_key].max()
            except ValueError:
                max_value = 0
            for inst in instance_ids:
                instance_id_map[inst][metrics_key] = max_value
        return instance_id_map
    
    def get_max_mem(self, grouped, row, ids, 
                    instance_id_map,
                    metadata, metrics_key='mem_used_gb'):
        '''Get max mem for a given subset of the runs'''
        runtime_results = metadata.copy()
        for id in ids:
            runtime_results = runtime_results[runtime_results[id] == row[id]].copy()
        max_values = []
        for inst in instance_id_map:
            if (row.workflow_id in instance_id_map[inst]['workflow_ids']) \
                    and (row.task_call_name == instance_id_map[inst]['task_call_name']):
                max_values.append(instance_id_map[inst][metrics_key])
        try:
            return max(max_values)
        except ValueError:
            return 0
        
    def get_flow_runtime(self, grouped, row, ids):
        ''' get max runtime in hours for group.
        For sub workflow and workflows this will also include wait time and
        return the end-to-end wallclock time for a particular id.
        For tasks this will only return the longest runtime.
        '''
        if 'task_call_name' in ids:
            runtime_results = grouped.run_time_m.agg(['max']).reset_index()
            for id in ids:
                runtime_results = runtime_results[runtime_results[id] == row[id]].copy()
            runtime = runtime_results['max'].tolist()[0]
            return runtime / 60.0
        else:
            end_time_results = grouped.end_time.agg(['max']).reset_index()
            for id in ids:
                end_time_results = end_time_results[end_time_results[id] == row[id]].copy()
            end_time = end_time_results['max'].tolist()[0]
            start_time_results = grouped.start_time.agg(['min']).reset_index()
            for id in ids:
                start_time_results = start_time_results[start_time_results[id] == row[id]].copy()
            start_time = start_time_results['min'].tolist()[0]
            runtime = end_time - start_time
        return runtime.total_seconds() / 60.0 / 60.0
        
    def load_plot_metrics(self, metadata, instance_id_map):
        '''add cpu_time, max_mem, wall_clock, cpu_used_percent, mem_efficiency, mem_used_gb'''
        # get wallclock
        metadata['start_time'] = pd.to_datetime(metadata.start_time, infer_datetime_format=True)
        metadata['end_time'] = pd.to_datetime(metadata.end_time, infer_datetime_format=True)
        metadata['run_time'] = metadata['end_time'] - metadata['start_time']
        metadata['run_time_m'] = metadata.apply(lambda row:
                                                     (row.run_time.total_seconds() / 60.0), axis=1)
        # get cpu_time
        metadata['cpu_time_m'] = metadata.apply(lambda row:
                                                (row.run_time_m
                                                 * row.cpu_count),
                                                axis=1)
        # ==============
        #     Tasks
        # ==============
        grouped = metadata.groupby('task_call_name')
        metadata['mean_task_core_h'] = metadata.apply(lambda row:
                                                           (row.cpu_count
                                                            * float(grouped.mean().run_time_m[row.task_call_name])) / 60,
                                                           axis=1)
        # task wallclock time from start to finish
        metadata['mean_task_run_time_h'] = metadata.apply(lambda row:
                                                          float(grouped.mean().run_time_m[row.task_call_name]) / 60,
                                                          axis=1)
        
        # ============================
        #     Tasks w/in subworkflow
        # ============================
        grouped = metadata.groupby(['id', 'task_call_name', 'workflow_name', 'instance_name'])
        # sub task wallclock time from start to finish
        metadata['sample_task_run_time_h'] = metadata.apply(lambda row: 
                                                            self.get_flow_runtime(grouped, row,
                                                                                  ids=['id', 
                                                                                       'instance_name',
                                                                                       'task_call_name',
                                                                                       'workflow_name']),
                                                            axis=1)
        metadata['max_mem_g'] = metadata.apply(lambda row: self.get_max_mem(grouped, row,
                                                                            metadata=metadata,
                                                                            metrics_key='mem_used_gb',
                                                                            instance_id_map=instance_id_map,
                                                                            ids=['id',
                                                                                 'instance_name',
                                                                                 'task_call_name',
                                                                                 'workflow_name']),
                                                                                 axis=1)
        metadata['sample_task_core_h'] = metadata.apply(lambda row:
                                                        self.get_task_core_h(grouped, row,
                                                                             ids=['id',
                                                                                  'instance_name',
                                                                                  'task_call_name',
                                                                                  'workflow_name']),
                                                        axis=1)
        # ============================
        #       sub workflows
        # ============================
        grouped = metadata.groupby(['id', 'workflow_name'])
        sub_metadata = grouped.agg(sample_subworkflow_core_h=pd.NamedAgg(column='sample_task_core_h', aggfunc=sum)).reset_index()
        metadata['sample_subworkflow_core_h'] = metadata.apply(lambda row:
                                                               sub_metadata[(sub_metadata.id == row.id)
                                                                            & (sub_metadata.workflow_name == row.workflow_name)].sample_subworkflow_core_h.tolist()[0], 
                                                               axis=1)
        # sub workflows wallclock time from start to finish
        metadata['sample_subworkflow_run_time_h'] = metadata.apply(lambda row:
                                                                   self.get_flow_runtime(grouped, row,
                                                                                         ids=['id', 'workflow_name']), 
                                                                   axis=1)
        metadata['subworkflow_max_mem_g'] = metadata.apply(lambda row: metadata[(metadata.id == row.id)
                                                                                & (metadata.workflow_name == row.workflow_name)].max_mem_g.max(),
                                                                                axis=1)
        # ============================
        #       workflows
        # ============================
        # full workflow wallclock time from start to finish for an id
        grouped = metadata.groupby('id')
        sub_metadata = grouped.agg(sample_workflow_core_h=pd.NamedAgg(column='sample_task_core_h', aggfunc=sum)).reset_index()
        metadata['sample_workflow_core_h'] = metadata.apply(lambda row: 
                                                            sub_metadata[(sub_metadata.id == row.id)].sample_workflow_core_h.tolist()[0],
                                                            axis=1)
        metadata['sample_workflow_run_time_h'] = metadata.apply(lambda row: 
                                                                self.get_flow_runtime(grouped, row, 
                                                                                      ids=['id']), 
                                                                axis=1)
        metadata['workflow_max_mem_g'] = metadata.apply(lambda row: metadata[(metadata.id == row.id)].max_mem_g.max(),
                                                        axis=1) 
        return metadata
    
    def load_end_date(self, date_record):
        return date_record.strftime(format='%Y-%m-%d')      

    def modify_date(self, run_date):
        return '-'.join(run_date.split('-')[0:3])
    
    def get_sample_info(self):
        self.sample_count = len(self.output_info['project_data']['sampleIds'])
        self.pair_count = len(self.output_info['project_data']['pairId'])
        
    def get_rundate(self):
        return self.modify_date(self.output_info['project_data']['run_date']) 
        
    def get_gcp_project(self):
        return self.output_info['project_data']['options']['monitoring_image'].split('/')[1]

    def format_sub_workflow(self, sub_workflow_uuids):
        sub_workflow_uuid_list = '("' + '", "'.join(sub_workflow_uuids) + '")'
        return sub_workflow_uuid_list
    
    def load_metadata(self):
        '''convert big query metadata table into pandas dataframe
            Includes start_time, end_time, execution_status, cpu_count,
            mem_total_gb, task_call_name, workflow_name
        '''
        log.info('Loading Metadata...')
        
        sub_workflow_uuid_list = self.format_sub_workflow(self.sub_workflow_uuids.sub_workflow_uuid.unique())
        
        query_string = '''
        SELECT * FROM `''' + self.gcp_project + '''.cromwell_monitoring.metadata`
            WHERE DATE(start_time) >= "''' + self.run_date + '''"
            AND workflow_id IN ''' + sub_workflow_uuid_list
            
#             + '''"
#             LIMIT ''' + self.limit + '''
#         '''
        self.metadata =  (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        self.end_date = self.load_end_date(self.metadata.end_time.max())
        sub_workflow_uuid_sample_map = self.sub_workflow_uuids.copy()
        sub_workflow_uuid_sample_map.columns = ['workflow_id', 'id']
        sub_workflow_uuid_sample_map = sub_workflow_uuid_sample_map.drop_duplicates().copy()
        self.metadata['inputs'] = self.metadata.apply(lambda row: str(row.inputs).replace('\n', ' '), axis=1)
        unhashable = ['disk_mounts', 'disk_total_gb', 'disk_types', 'inputs']
        # hashable = [col for col in self.metadata.columns if not col in unhashable]
        hashable = ['instance_name', 'workflow_id', 'attempt']
        self.metadata = self.metadata.drop_duplicates(subset=hashable)
        self.metadata = self.metadata.merge(sub_workflow_uuid_sample_map, on='workflow_id')
        
    def load_runtime(self):
        '''convert big query runtime table into pandas dataframe
        Includes both instance_id and workflow_id.
        Use to map workflow_id to instance_ids
        '''
        log.info('Loading Runtime...')
        query_string = '''
        SELECT * FROM `''' + self.gcp_project + '''.cromwell_monitoring.runtime`
            WHERE DATE(start_time) >= "''' + self.run_date + '''"
            AND DATE(start_time) <= "''' + self.end_date + '''"
            LIMIT ''' + self.limit + '''
        '''
        runtime = (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        sub_workflow_uuids = self.sub_workflow_uuids.sub_workflow_uuid.tolist()
        self.runtime = runtime[runtime.workflow_id.isin(sub_workflow_uuids)].copy()
        self.instance_ids = self.runtime.instance_id.unique().tolist()
        self.instance_id_map = {}
        for instance_id in self.instance_ids:
            self.instance_id_map[instance_id] = {}
            self.instance_id_map[instance_id]['workflow_ids'] = self.runtime[self.runtime.instance_id == instance_id].workflow_id.unique().tolist()
            self.instance_id_map[instance_id]['task_call_name'] = self.runtime[self.runtime.instance_id == instance_id].task_call_name.tolist()[0]
        self.non_retry_instance_id_map = {}
        map = {row.instance_name: row.execution_status for index, row in self.metadata[['instance_name', 'execution_status']].drop_duplicates().iterrows()}
        execution_status_map = {row.instance_id : map[row.instance_name] for index, row in 
                                self.runtime[['instance_name', 'instance_id']].drop_duplicates().iterrows()}
        non_retry_instance_ids = [inst for inst in execution_status_map]
        for instance_id in non_retry_instance_ids:
            self.non_retry_instance_id_map[instance_id] = {}
            self.non_retry_instance_id_map[instance_id]['workflow_ids'] = self.runtime[(self.runtime.instance_id == instance_id)].workflow_id.unique().tolist()
            self.non_retry_instance_id_map[instance_id]['task_call_name'] = self.runtime[(self.runtime.instance_id == instance_id)].task_call_name.tolist()[0]
        
    def format_instance_id(self):
        instance_id_list = '(' + ', '.join([str(i) for i in self.instance_ids]) + ')'
        return instance_id_list
    
    def load_metrics(self):
        '''convert big query metrics table into pandas dataframe
        Includes instance_id, cpu_used_percent, mem_used_gb, disk_used_gb, 
        disk_read_iops, disk_write_iops.
        Use to find max and plot resource usage over time.
        '''
        log.info('Loading Metrics...')
        
        instance_id_list = self.format_instance_id()
        
        query_string = '''
        SELECT * FROM `''' + self.gcp_project + '''.cromwell_monitoring.metrics`
            WHERE DATE(timestamp) >= "''' + self.run_date + '''"
            AND DATE(timestamp) <= "''' + self.end_date + '''"
            AND instance_id IN ''' + instance_id_list 
        self.metrics = (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )

    def load_sub_workflow_uuids(self):
        sub_workflow_uuids_list = []
        for id in self.output_info['sub_workflow_uuids']:
            data = pd.DataFrame({'sub_workflow_uuid' : self.output_info['sub_workflow_uuids'][id]})
            data['id'] = id
            sub_workflow_uuids_list.append(data)
        sub_workflow_uuids = pd.concat(sub_workflow_uuids_list, 
                                       ignore_index=False)
        return sub_workflow_uuids

    
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-info-file',
                        help='JSON file with outputInfo. Includes workflow uuid, '
                        'options JSON dictionary submitted to cromwell, '
                        'project pairing, sample, '
                        'genome build, library, output files, '
                        'subworkflow uuids, id sorted output files'
                        'and interval list information',
                        required=False
                        )
    parser.add_argument('--multi',
                        help='Records for multiple runs are stored in the JSON file with outputInfo.',
                        required=False,
                        action='store_true'
                        )
    parser.add_argument('--name',
                        default='nygc_pipeline',
                        help='Use with the --multi flag to name the output.',
                        required=False
                        )
    parser.add_argument('--manifest',
                        help='Use with the --multi flag to list outputInfo JSON files.',
                        required=False
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    output_infos = []
    manifest_files = []
    if args['multi']:
        with open(args['manifest']) as manifest:
            for line in manifest:
                output_infos.append(Runtime.read(line.strip()))
                manifest_files.append(line.strip())
    else:
        output_info = Runtime.read(args['output_info_file'])
        output_infos = [output_info]
    non_retried_metrics = []
    metrics = []
    for output_info in output_infos:
        runtime = Runtime(limit=10000,
                          output_info=output_info)
        non_retried_metrics.append(runtime.non_retried_metrics_file)
        metrics.append(runtime.metrics_file)
    if args['multi']:
        manifest = pd.DataFrame({'non_retried_metrics' : non_retried_metrics,
                                 'metrics' : metrics,
                                 'output_info' : manifest_files})
        manifest.to_csv(args['name'] + '_MetricsInfoManifest.csv', index=False)
    

if __name__ == "__main__":
    main()       
    