import glob
import pandas as pd
import numpy as np
import sys

def add_cost(kind, row, map, costs):
    '''add various costs from table'''
    options = ['capacity_cost', 'core_cost', 'ram_cost', 'egress_cost']
    cost_keys = [key for key in map if map[key]['cost'] == kind]
    if row.main_workflow_name.lower() == row.sub_workflow_name:
        match = costs[(costs.wdl_task_name == row.wdl_task_name) &
                      (costs.workflow_id == row.costs_workflow_id) &
                      (costs.service_type.isin(cost_keys))].copy()
    else:
        match = costs[(costs.sub_workflow_name == row.sub_workflow_name) &
                      (costs.wdl_task_name == row.wdl_task_name) &
                      (costs.workflow_id == row.costs_workflow_id) &
                      (costs.service_type.isin(cost_keys))].copy()
    if not match.empty:
        # preemptible doesn't matter
        if kind in ['capacity_cost', 'egress_cost']:
            return match.cost.sum()
        elif kind in ['core_cost', 'ram_cost']:
            cost_keys = [key for key in cost_keys 
                         if str(map[key]['preemptible']) == str(row.preemptible)]
            match = match[(match.service_type.isin(cost_keys))].copy()
            if not match.empty:
                return match.cost.sum()
    return np.nan

def add_type(row, map, costs):
    '''add instance type from cost info'''
    options = ['core_cost']
    cost_keys = [key for key in map if map[key]['cost'] == 'core_cost']
    if row.main_workflow_name.lower() == row.sub_workflow_name:
        match = costs[(costs.wdl_task_name == row.wdl_task_name) &
                      (costs.service_type.isin(cost_keys))].copy()
    else:
        match = costs[(costs.sub_workflow_name == row.sub_workflow_name) &
                      (costs.wdl_task_name == row.wdl_task_name) &
                      (costs.service_type.isin(cost_keys))].copy()
    if not match.empty:
        # preemptible doesn't matter
        if str(row.preemptible) != 'nan':
            cost_keys = [key for key in cost_keys 
                         if str(map[key]['preemptible']) == str(row.preemptible)]
            match = match[(match.service_type.isin(cost_keys))].copy()
            if not match.empty:
                hits = [map[i]['type'] for i in match.service_type.unique().tolist()]
                return ','.join(hits)
    return np.nan

def get_cost(value, length):
    if length == 0:
        return np.nan
    return value / length

def load_runtime(file, costs, uuid):
    runtime = pd.read_csv(file)
    map = {'Spot Preemptible Custom Instance Core running in Americas' : {'cost' : 'core_cost',
                                                                    'type' : 'spot_custom',
                                                                    'preemptible' : True},
           'Spot Preemptible Custom Instance Ram running in Americas' : {'cost' : 'ram_cost',
                                                                    'type' : 'spot_custom',
                                                                    'preemptible' : True},
           'SSD backed Local Storage attached to Spot Preemptible VMs': {'cost' : 'capacity_cost',
                                      'type' : 'ssd_local_attached_spot_capacity',
                                      'preemptible' : np.nan},
            'SSD backed PD Capacity': {'cost' : 'capacity_cost',
                                      'type' : 'ssd_pd_capacity',
                                      'preemptible' : np.nan},
           'Custom Instance Core running in Americas': {'cost' : 'core_cost',
                                                        'type' : 'custom',
                                                        'preemptible' : False},
           'Preemptible Custom Instance Core running in Americas': {'cost' : 'core_cost',
                                                                    'type' : 'custom',
                                                                    'preemptible' : True},
           'Custom Instance Ram running in Americas': {'cost' : 'ram_cost',
                                                       'type' : 'custom',
                                                       'preemptible' : False},
           'Preemptible Custom Instance Ram running in Americas' : {'cost' : 'ram_cost',
                                                                    'type' : 'custom',
                                                                    'preemptible' : True},
           'Storage PD Capacity': {'cost' : 'capacity_cost',
                                   'type' : 'storage_pd_capacity',
                                   'preemptible' : np.nan},
           'Preemptible N1 Predefined Instance Core running in Americas': {'cost' : 'core_cost',
                                                                    'type' : 'N1_predefined',
                                                                    'preemptible' : True},
           'Preemptible N1 Predefined Instance Ram running in Americas': {'cost' : 'ram_cost',
                                                                    'type' : 'N1_predefined',
                                                                    'preemptible' : True},
           'Network Internet Egress from Americas to Americas': {'cost' : 'egress_cost',
                                                                 'type' : 'network_egress',
                                                                 'preemptible' : np.nan}
           }
    runtime['main_workflow_id'] = runtime.apply(lambda row: uuid, axis=1)
    runtime['costs_workflow_id'] = runtime.apply(lambda row: 'cromwell-' + row.main_workflow_id, axis=1)
    runtime['sub_workflow_name'] = runtime.apply(lambda row: row.workflow_name.lower(), axis=1)
    runtime['wdl_task_name'] = runtime.apply(lambda row: row.non_alias_task_call_name.split('.')[-1].lower(), axis=1)
    options = ['capacity_cost', 'core_cost', 'ram_cost', 'egress_cost']
    for option in options:
        runtime['total_' + option] = runtime.apply(lambda row: add_cost(option, row, map, costs), axis=1)
        runtime['count_' + option] = runtime.groupby(['main_workflow_id', 
                                                      'sub_workflow_name', 
                                                      'wdl_task_name'])['total_' + option].transform('count')
        runtime[option] = runtime.apply(lambda row: get_cost(row['total_' + option], 
                                                             row['count_' + option]), axis=1)
        remove = ['total_' + option, 'count_' + option]
        remain = [col for col in runtime.columns if not col in remove]
        runtime = runtime[remain]
    runtime['service_type'] = runtime.apply(lambda row: add_type(row, map, costs), axis=1)
    return runtime

#get ${uuid}
uuid = sys.argv[1]
# get ${uuid}_outputCosts.csv file
costs_file = sys.argv[2]
# get ${project_name}.${uuid}_outputMetrics.csv
runtime_file = sys.argv[3]
output_prefix = sys.argv[4]
# process
costs = pd.read_csv(costs_file)
runtime = load_runtime(runtime_file, costs, uuid)
runtime.to_csv(output_prefix + '.outputMetrics.cost.csv', index=False)
total = pd.DataFrame({'workflow_uuid' : [uuid],
                      'cost' : [costs.cost.sum()]})
total.to_csv(output_prefix + '.outputMetrics.total.csv', index=False)

    
# output_dir = '/gpfs/commons/projects/TCGA/gdc-awg/ALCHEMIST/WGS/analysis/v7/gcp/tmp/'
# # Ataq seq
# cost_glob = '/gpfs/internal/sweng/wdl_cost/cost_test/*_outputCosts.csv'
# input_glob = '/gpfs/commons/home/jroche/RND300/outputMetrics_expandMetadata/*'
# output_dir = '/gpfs/commons/projects/TCGA/gdc-awg/ALCHEMIST/WGS/analysis/v7/gcp/tmp/'
# # output = output_dir + 'ataq_seq.outputMetrics.cost.csv'
# # dfs= [] 
# # all_costs = []   
# # for file in glob.glob(cost_glob):
# #     uuid = file.split('/')[-1].replace('_outputCosts.csv', '')
# #     costs = pd.read_csv(file)
# #     runtime_file = glob.glob('/gpfs/commons/home/jroche/RND300/outputMetrics_expandMetadata/*' + uuid +  '_outputMetrics.csv')[0]
# #     runtime = load_runtime(runtime_file, costs, uuid)
# #     dfs.append(runtime)
# #     all_costs.append(costs)
# # all = pd.concat(dfs, ignore_index=True)
# # all.to_csv(output, index=False)
# # all_cost = pd.concat(all_costs, ignore_index=True)
# # all_cost.to_csv(output_dir + 'ataq_seq.cost.csv')
# #  
# # alchemist
# cost_glob = '/gpfs/internal/sweng/wdl_cost/cost_test_1/*_outputCosts.csv'
# output_dir = '/gpfs/commons/projects/TCGA/gdc-awg/ALCHEMIST/WGS/analysis/v7/gcp/tmp/'
# output = output_dir + 'alchemist.outputMetrics.cost.csv'
# dfs = []
# for file in glob.glob(cost_glob):
#     uuid = file.split('/')[-1].replace('_outputCosts.csv', '')
#     costs = pd.read_csv(file)
#     runtime_file = glob.glob('/gpfs/commons/projects/TCGA/gdc-awg/ALCHEMIST/WGS/analysis/v7/gcp/ALCHEMIST.*' + uuid +  '_outputMetrics.csv')[0]
#     runtime = load_runtime(runtime_file, costs, uuid)
#     dfs.append(runtime)
# all = pd.concat(dfs, ignore_index=True)
# all.to_csv(output, index=False)
# 
# 
# # alchemist meta
# output = output_dir + 'alchemist.outputMetrics.cost.csv'
# cost = pd.read_csv(output)
# file = '/gpfs/commons/projects/TCGA/gdc-awg/ALCHEMIST/metadata/manifest_map.csv'
# metadata = pd.read_csv(file)
# type_and_cov_file = '/gpfs/commons/projects/TCGA/gdc-awg/ALCHEMIST/metadata/pairs_with_type.csv'
# type_and_cov = pd.read_csv(type_and_cov_file)
# type_and_cov['is_ffpe'] = type_and_cov.apply(lambda row: row.tumor, axis=1)
#pon

#  runtime[(runtime.wdl_task_name == 'svabawgs') &
#          (runtime.sub_workflow_name == 'svaba')]
#     
