import sys
from google.cloud import bigquery
from google.cloud import bigquery_storage
import json
import datetime
import logging as log
import pandas as pd
import Colorer

pd.set_option('display.max_columns', 500)
log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class Runtime():
    def __init__(self,
                 limit=10000,
                 gcp_project=False,
                 output_info_file=False,
                 metrics_file=False,
                 sample_ids=False,
                 pair_ids=False,
                 run_date=False,
                 sub_workflow_uuids=False):
        assert sub_workflow_uuids or output_info_file, 'Define either output_info_file or sub_workflow_uuids'
        if not output_info_file:
            assert run_date, 'Run date must be defined if task uuids is defined'
            assert pair_ids, 'pairIds must be defined if task uuids is defined'
            assert sample_ids, 'sampleIds must be defined if task uuids is defined'
        if not metrics_file:
            self.metrics_file = output_info_file.replace('_outputInfo.json', '_outputMetrics.csv')
        else:
            self.metrics_file = metrics_file
        # find relevant project info
        self.limit = str(limit)
        self.metrics_limit = str(limit * 4)
        self.gcp_project = self.get_gcp_project(gcp_project, output_info_file)
        self.sub_workflow_uuids = self.load_sub_workflow_uuids(sub_workflow_uuids, output_info_file)
        self.run_date = self.get_rundate(run_date, output_info_file)
        self.get_sample_info(sample_ids, pair_ids,
                             output_info_file)
        # gather big cloud tables
        self.bqclient = bigquery.Client()
        self.bqstorageclient = bigquery_storage.BigQueryReadClient()
        self.load_metadata()
        self.load_runtime()
        self.load_metrics()
        self.instance_id_map = self.get_max_mems(self.metadata, 
                                                 instance_id_map=self.instance_id_map, 
                                                 execution_status=False,
                                                 metrics_key='mem_used_gb')
        self.non_retry_instance_id_map = self.get_max_mems(self.metadata, 
                                                           instance_id_map=self.non_retry_instance_id_map, 
                                                           execution_status='Done',
                                                           metrics_key='mem_used_gb')
        self.non_retry_metadata = self.metadata[self.metadata.execution_status == 'Done'].copy()
        self.metadata = self.load_plot_metrics(self.metadata,
                                               instance_id_map=self.instance_id_map)
        self.non_retry_metadata = self.load_plot_metrics(self.non_retry_metadata,
                                                         instance_id_map=self.non_retry_instance_id_map)
        self.metadata.to_csv(self.metrics_file, index=False)
        self.non_retry_metadata.to_csv(self.metrics_file.replace('.csv', '.non_retried.csv'), index=False)
        
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
            max_value = self.metrics[self.metrics.instance_id.isin(instance_ids)][metrics_key].max()
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
        return max(max_values)
        
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

    def modify_date(self, run_date):
        return '-'.join(run_date.split('-')[0:3])
    
    def get_sample_info(self, sample_ids, pair_ids, 
                        output_info_file):
        if sample_ids:
            self.sample_count = len(sample_ids)
            self.pair_count = len(pair_ids)
        else:
            with open(output_info_file) as output_info_object:
                output_info = json.load(output_info_object)
            self.sample_count = len(output_info['project_data']['sampleIds'])
            self.pair_count = len(output_info['project_data']['pairId'])
        
    def get_rundate(self, run_date, output_info_file):
        if run_date:
            return self.modify_date(run_date)
        else:
            with open(output_info_file) as output_info_object:
                output_info = json.load(output_info_object)
            return self.modify_date(output_info['project_data']['run_date']) 
        
    def get_gcp_project(self, gcp_project, output_info_file):
        if gcp_project:
            return gcp_project
        else:
            with open(output_info_file) as output_info_object:
                output_info = json.load(output_info_object)
            return output_info['project_data']['options']['monitoring_image'].split('/')[1]

    
    def load_metadata(self):
        '''convert big query metadata table into pandas dataframe
            Includes start_time, end_time, execution_status, cpu_count,
            mem_total_gb, task_call_name, workflow_name
        '''
        log.info('Loading Metadata...')
        query_string = '''
        SELECT * FROM `''' + self.gcp_project + '''.cromwell_monitoring.metadata`
            WHERE DATE(start_time) >= "''' + self.run_date + '"'
#             + '''"
#             LIMIT ''' + self.limit + '''
#         '''
        meta = (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        sub_workflow_uuid_sample_map = self.sub_workflow_uuids.copy()
        sub_workflow_uuid_sample_map.columns = ['workflow_id', 'id']
        sub_workflow_uuid_sample_map = sub_workflow_uuid_sample_map.drop_duplicates().copy()
        self.metadata = meta[meta.workflow_id.isin(sub_workflow_uuid_sample_map.workflow_id.tolist())].copy()
        self.metadata['inputs'] = self.metadata.apply(lambda row: str(row.inputs).replace('\n', ' '), axis=1)
        unhashable = ['disk_mounts', 'disk_total_gb', 'disk_types', 'inputs']
        hashable = [col for col in self.metadata.columns if not col in unhashable]
        self.metadata = self.metadata.drop_duplicates(hashable)
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
        execution_status_map = {row.instance_id : map[row.instance_name] for index, row in self.runtime[['instance_name', 'instance_id']].drop_duplicates().iterrows()}
        non_retry_instance_ids = [inst for inst in execution_status_map]
        for instance_id in non_retry_instance_ids:
            self.non_retry_instance_id_map[instance_id] = {}
            self.non_retry_instance_id_map[instance_id]['workflow_ids'] = self.runtime[(self.runtime.instance_id == instance_id)].workflow_id.unique().tolist()
            self.non_retry_instance_id_map[instance_id]['task_call_name'] = self.runtime[(self.runtime.instance_id == instance_id)].task_call_name.tolist()[0]
        
    def load_metrics(self):
        '''convert big query metrics table into pandas dataframe
        Includes instance_id, cpu_used_percent, mem_used_gb, disk_used_gb, 
        disk_read_iops, disk_write_iops.
        Use to find max and plot resource usage over time.
        '''
        log.info('Loading Metrics...')
        query_string = '''
        SELECT * FROM `''' + self.gcp_project + '''.cromwell_monitoring.metrics`
            WHERE DATE(timestamp) >= "''' + self.run_date + '''"
        '''
        metrics = (
            self.bqclient.query(query_string)
            .result()
            .to_dataframe(bqstorage_client=self.bqstorageclient)
        )
        self.metrics = metrics[metrics.instance_id.isin(self.instance_ids)].copy()

    def load_sub_workflow_uuids(self, sub_workflow_uuids, output_info_file):
        if sub_workflow_uuids:
            return sub_workflow_uuids
        else:
            with open(output_info_file) as output_info_object:
                output_info = json.load(output_info_object)
            sub_workflow_uuids_list = []
            for id in output_info['sub_workflow_uuids']:
                data = pd.DataFrame({'sub_workflow_uuid' : output_info['sub_workflow_uuids'][id]})
                data['id'] = id
                sub_workflow_uuids_list.append(data)
            sub_workflow_uuids = pd.concat(sub_workflow_uuids_list)
            return sub_workflow_uuids
           

def main():
    output_info_file = sys.argv[2]
    Runtime(limit=10000,
            output_info_file=output_info_file,
            run_date=False,
            sub_workflow_uuids=False)
    
    
if __name__ == "__main__":
    main()
    