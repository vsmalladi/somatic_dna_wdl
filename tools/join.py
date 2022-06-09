import glob
import pandas as pd
import numpy as np
import sys
import logging as log

try:
    import Colorer
except ImportError:
    pass


def is_preemptible(service):
    if 'preemptible' in service:
        return True
    else:
        return False


def formatter(service, drop=[], truncate=False, suffix=False):
    # takes a string service type, split to list of words, drop some useless words,
    # drop up to truncate word and suffix. can't just keep it do to 'storage'-> 'capacity'
    word_list = service.split()
    clean_words = [a for a in word_list if a not in drop]

    if truncate:
        try:
            truncate_index = clean_words.index(truncate)
        except ValueError:  # storage services are called capacity
            truncate_index = -1
        words = clean_words[:truncate_index]
    if suffix:
        words.append(suffix)

    return '_'.join(words)


def mapper(row):
    # cost_type & preemptible are pretty clearly defined. cost_description not so much
    service_type = row['service_type'].lower()

    if 'storage' in service_type or 'capacity' in service_type:
        cost_type = 'capacity_cost'
        preemptible = np.nan
        cost_description = formatter(service_type,
                                     drop=['backed', 'storage', 'to', 'preemptible', 'vms'],
                                     truncate='capacity',
                                     suffix='capacity'
                                     )
    elif 'egress' in service_type:
        cost_type = 'egress_cost'
        preemptible = np.nan
        cost_description = formatter(service_type,
                                     truncate='egress',
                                     suffix='egress'
                                     )
    elif 'core' in service_type:
        cost_type = 'core_cost'
        preemptible = is_preemptible(service_type)
        cost_description = formatter(service_type,
                                     drop=['preemptible'],
                                     truncate='instance'
                                     )
    elif 'ram' in service_type:
        cost_type = 'ram_cost'
        preemptible = is_preemptible(service_type)
        cost_description = formatter(service_type,
                                     drop=['preemptible'],
                                     truncate='instance'
                                     )
    else:
        cost_type = 'unknown_cost'
        cost_description = 'unknown'
        preemptible = np.nan

    return cost_type, cost_description, preemptible


def load_runtime(file, costs, uuid):
    runtime = pd.read_csv(file)

    runtime['main_workflow_id'] = uuid
    runtime['cromwell_workflow_id'] = 'cromwell-' + runtime['main_workflow_id']

    unique_roots = runtime.drop_duplicates(['main_workflow_name', 'cromwell_workflow_id']).set_index(
        'cromwell_workflow_id')
    costs['sub_workflow_name'] = np.where(costs['sub_workflow_name'].isna(),
                                          costs['cromwell_workflow_id'].map(unique_roots['main_workflow_name']),
                                          costs['sub_workflow_name']
                                          )

    costs[['cost_type', 'cost_desc', 'preemptible']] = costs.apply(mapper, axis=1, result_type="expand")
    grouper = ['cromwell_workflow_id', 'sub_workflow_name', 'wdl_task_name']
    wide_costs = pd.pivot_table(costs,
                                values='cost',
                                index=grouper,
                                columns='cost_type',
                                aggfunc='sum',
                                fill_value=0
                                ).reset_index()

    runtime = pd.merge(runtime, wide_costs, on=grouper, how='outer')

    cost_opts = ['capacity_cost', 'core_cost', 'ram_cost', 'egress_cost', 'unknown_cost']
    cost_cols = [c for c in cost_opts if c in runtime.columns]

    # the costs are still for the full group, divide them evenly
    runtime['group_count'] = runtime.groupby(grouper)['task_call_name'].transform('count')
    for c in cost_cols:
        runtime[c] = runtime[c] / runtime['group_count']
    runtime['total_cost'] = runtime[cost_cols].sum(axis=1)

    runtime['group_cost'] = runtime.groupby(grouper)['total_cost'].transform('sum')
    runtime['group_runtime'] = runtime.groupby(grouper)['actual_runtime_m'].transform('sum')
    runtime['runtime_scaled_total_cost'] = runtime['group_cost'] * (runtime['actual_runtime_m'] / runtime['group_runtime'])

    # add the service type, treat as set. Almost always going to be a single type
    core_costs = costs[costs['cost_type'] == 'core_cost'].copy()
    service_types = core_costs.groupby(grouper).agg({'cost_desc': set}).reset_index()
    runtime = pd.merge(runtime, service_types, on=grouper, how='outer')

    missing = costs[costs['cost_type'] == 'unknown_cost']['service_type'].unique().tolist()
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
