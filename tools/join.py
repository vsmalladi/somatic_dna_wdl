import glob
import pandas as pd
import numpy as np
import sys
import logging as log

try:
    import Colorer
except ImportError:
    pass

def load_runtime(file, costs, uuid):
    runtime = pd.read_csv(file)
    map = {
        'Network Internet Egress from Americas to EMEA':
            {'cost': 'egress_cost',
             'type': 'network_internet_egress',
             'preemptible': np.nan
             },
        'Spot Preemptible Custom Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'spot_custom',
             'preemptible': True
             },
        'Spot Preemptible Custom Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'spot_custom',
             'preemptible': True
             },

        'SSD backed Local Storage':
            {'cost': 'capacity_cost',
             'type': 'ssd_local_capacity',
             'preemptible': np.nan
             },
        'SSD backed Local Storage attached to Spot Preemptible VMs':
            {'cost': 'capacity_cost',
             'type': 'ssd_local_attached_spot_capacity',
             'preemptible': np.nan
             },
        'SSD backed PD Capacity':
            {'cost': 'capacity_cost',
             'type': 'ssd_pd_capacity',
             'preemptible': np.nan
             },
        'Custom Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'custom',
             'preemptible': False
             },
        'Preemptible Custom Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'custom',
             'preemptible': True
             },
        'Custom Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'custom',
             'preemptible': False
             },
        'Preemptible Custom Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'custom',
             'preemptible': True
             },
        'Storage PD Capacity':
            {'cost': 'capacity_cost',
             'type': 'storage_pd_capacity',
             'preemptible': np.nan},
        'Spot Preemptible N1 Predefined Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'spot_N1_predefined',
             'preemptible': True
             },
        'Spot Preemptible N1 Predefined Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'spot_N1_predefined',
             'preemptible': True
             },
        'Spot N1 Predefined Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'spot_N1_predefined',
             'preemptible': False
             },
        'Spot N1 Predefined Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'spot_N1_predefined',
             'preemptible': False
             },
        'N1 Predefined Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'N1_predefined',
             'preemptible': False
             },
        'N1 Predefined Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'N1_predefined',
             'preemptible': False
             },
        'N2 Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'N1_predefined',
             'preemptible': False
             },
        'N2 Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'N1_predefined',
             'preemptible': False
             },
        'Preemptible N1 Predefined Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'N1_predefined',
             'preemptible': True
             },
        'Preemptible N1 Predefined Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'N1_predefined',
             'preemptible': True
             },
        'Spot Preemptible N2 Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'spot_N2',
             'preemptible': True
             },
        'Spot Preemptible N2 Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'spot_N2',
             'preemptible': True
             },
        'Preemptible N2 Instance Core running in Americas':
            {'cost': 'core_cost',
             'type': 'N2',
             'preemptible': True
             },
        'Preemptible N2 Instance Ram running in Americas':
            {'cost': 'ram_cost',
             'type': 'N2',
             'preemptible': True
             },
        'Network Internet Egress from Americas to Americas':
            {'cost': 'egress_cost',
             'type': 'network_internet_egress',
             'preemptible': np.nan
             },
        'Network Inter Region Egress from Americas to Americas':
            {'cost': 'egress_cost',
             'type': 'network_inter_region_egress',
             'preemptible': np.nan
             },
        'Network Inter Zone Egress':
            {'cost': 'egress_cost',
             'type': 'network_inter_zone_egress',
             'preemptible': np.nan
             },
    }

    runtime['main_workflow_id'] = uuid
    runtime['cromwell_workflow_id'] = 'cromwell-' + runtime['main_workflow_id']

    # # note: tasks from root_workflow do not have 'sub_workflow_name'. fill this based on runtime
    # # join.py only takes a single uuid so should only have one unique 'main_workflow_name':
    # sub_workflow_name = runtime['main_workflow_name'].iloc[0]
    # costs['sub_workflow_name'].fillna(sub_workflow_name, inplace=True)
    # A more generalizable approach mapping unique uuids.(unlikely but just to be safe or if multiple uuids allowed)
    unique_roots = runtime.drop_duplicates(['main_workflow_name', 'cromwell_workflow_id']).set_index('cromwell_workflow_id')
    costs['sub_workflow_name'] = np.where(costs['sub_workflow_name'].isna(),
                                          costs['cromwell_workflow_id'].map(unique_roots['main_workflow_name']),
                                          costs['sub_workflow_name']
                                          )

    costs['map'] = costs['service_type'].map(map)
    missing = costs[costs['map'].isna()]['service_type'].unique().tolist()
    # placeholder for items not in map so script doesn't fail. 
    missing_map = {'cost': 'unknown_cost',
                   'type': 'unknown',
                   'preemptible': np.nan
                   }
    costs['map'] = np.where(costs['map'].isna(), missing_map, costs['map'])
    costs['cost_type'] = costs['map'].apply(lambda x: x.get('cost'))
    costs['core_type'] = costs['map'].apply(lambda x: x.get('type'))

    grouper = ['cromwell_workflow_id', 'sub_workflow_name', 'wdl_task_name']
    wide_costs = pd.pivot_table(costs,
                                values='cost',
                                index=grouper,
                                columns='cost_type',
                                aggfunc='sum',
                                fill_value=0
                                ).reset_index()

    runtime = pd.merge(runtime, wide_costs, on=grouper, how='outer') # leaving outer to check validity. can change to left when done

    cost_opts = ['capacity_cost', 'core_cost', 'ram_cost', 'egress_cost', 'unknown_cost']
    cost_cols = [c for c in cost_opts if c in runtime.columns]

    # the costs are still for the full group, divide them evenly
    runtime['group_count'] = runtime.groupby(grouper)['task_call_name'].transform('count')
    for c in cost_cols:
        runtime[c] = runtime[c] / runtime['group_count']
    runtime['total_cost'] = runtime[cost_cols].sum(axis=1)


    # also give estimated total cost based on task runtimes, could also do this for each sub-category if desired
    runtime['group_cost'] = runtime.groupby(grouper)['total_cost'].transform('sum')
    runtime['group_runtime'] = runtime.groupby(grouper)['run_time_m'].transform('sum')
    runtime['task_cost_estimate'] = runtime['group_cost'] * (runtime['run_time_m'] / runtime['group_runtime'])

    # add the service type, treat as set. Almost always going to be a single type
    core_costs = costs[costs['cost_type'] == 'core_cost'].copy()
    service_types = core_costs.groupby(grouper).agg({'core_type': set}).reset_index()
    runtime = pd.merge(runtime, service_types, on=grouper, how='outer')

    if len(missing) > 0:
            log.warning('Following service_types are not yet described in map. Add them to join.py: ')
            log.warning(', '.join(missing))

    runtime.drop(columns=['group_count', 'group_cost', 'group_runtime'], inplace=True)
    return runtime


# get ${uuid}
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
total = pd.DataFrame({'workflow_uuid': [uuid],
                      'cost': [costs.cost.sum()]})
total.to_csv(output_prefix + '.outputMetrics.total.csv', index=False)
