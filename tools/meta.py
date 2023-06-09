'''
Created on Feb 8, 2021

@author: jshelton
'''
import os
import sys
import json
import argparse
import pandas as pd
import jsonschema
from jsonschema import validate
import logging as log
import pprint
import Colorer
import parse

log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)

def git_log(program, soft_fail=False):
    '''
        Function returns the most recent git commit and tag.
    '''
    import subprocess
    potential_git_repo = os.path.abspath(os.path.dirname(program))
    while True:
        if not os.path.isdir(potential_git_repo):
            log.error(')No parent directories of the program contain a .git directory :' + program )
            if soft_fail:
                return "", "not a git repo", "", ""
            else:
                sys.exit(1)
        if os.path.isdir(potential_git_repo + '/.git'):
            break
        else:
            potential_git_repo += '/..'
    git_dir = potential_git_repo + '/.git'
    work_tree_dir = potential_git_repo + '/.git'
    commit = subprocess.check_output(['git', '--git-dir=' + git_dir,
                                      '--work-tree=' + work_tree_dir,
                                      'log', '-1' ]).decode()
    try:
        tag = subprocess.check_output(['git', '--git-dir=' + git_dir,
                                       '--work-tree=' + work_tree_dir,
                                       'describe', '--abbrev=0', '--tags' ]).decode()
        uniq_tag = subprocess.check_output(['git', '--git-dir=' + git_dir,
                                            '--work-tree=' + work_tree_dir,
                                            'describe', '--tags' ]).decode()
    except subprocess.CalledProcessError:
        tag = ''
        uniq_tag = ''
    branch = subprocess.check_output(['git', '--git-dir=' + git_dir,
                                      '--work-tree=' + work_tree_dir,
                                      'branch' ])
    branch = str(branch).split()[1].split('\\n')[0]
    commit = commit.split('\n')[0].split()[-1]
    return commit, tag.split('\n')[0], uniq_tag.split('\n')[0], branch

def read(file):
    with open(file) as project_info_file:
        project_info = json.load(project_info_file)
        return project_info

def load_pairs(file):
    pairs = pd.read_csv(file)
    assert 'tumorId' in pairs.columns and 'normalId' in pairs.columns, 'Error: can not find tumorId and normalId columns in pairs file: ' + ' '.join(pairs.columns)
    pairs['pairId'] = pairs.apply(lambda row: row.tumorId + '--' + row.normalId, axis=1)
    return pairs[['tumorId', 'normalId', 'pairId']], pairs


def load_sample_ids(file):
    '''Return a deduplicated list of sample ids'''
    sample_ids = pd.read_csv(file)
    assert 'sampleId' in sample_ids.columns , 'Error: can not find "sampleId" columns in samples file: ' + ' '.join(sample_ids.columns)
    return sample_ids.sampleId.unique().tolist()


def fill_pair_relationship(row):
    '''Add pair and sample level info to object
    PairRelationship '''
    pair_relationship = {"pairId" : row.pairId}
    pair_relationship['normalId'] = row['normalId']
    pair_relationship['tumorId'] = row['tumorId']
    pair_relationship['normalPrefix'] = row['normalId']
    pair_relationship['tumorPrefix'] = row['tumorId']
    return pair_relationship      


def note_updates(key, new_value, project_info):
    ''' Add/replace value in project_info with value from flag'''
    if new_value:
        if key in project_info:
            if not project_info[key] == new_value:
                log.warning('Note ' + key + ' is already defined and differs from flag value:\n' + str(project_info[key]) + '\n from the value indicated with the flag:\n' + str(new_value) + '\nThe value from the flag supersedes the previous value')
        project_info[key] = new_value
    return project_info


def note_custom_updates(key, alt_project_info, project_info):
    ''' Add/replace value in project_info with value from custom inputs json'''
    match = [match_key for match_key in alt_project_info if match_key.split('.')[-1] == key]
    if len(match) == 1:
        if key in project_info:
            log.warning('Note that the value for ' + key + ' will be taken from --custom-inputs file\n')
        project_info[key] = alt_project_info[match[0]]
    return project_info


def verify_required(key, args, project_info):
    '''Verify a required key'''
    in_project_info = key in project_info and project_info[key]
    in_args = key in args and args[key]
    if not in_project_info and not in_args:
        log.error('Note that the value for ' + key + ' is required to be specified in flags or in the  --project-info file')
        sys.exit(1)
    return True


def fill_in_pair_info(project_info, full_pairs, suffix='.bai'):
    '''
    Read from input pairs CSV file
    REPLACE later with step to ingest standard sample sheet
    '''
    pair_infos = []
    normal_sample_infos = []
    if 'tumorBam' in full_pairs and 'normalBam' in full_pairs:
        for info in project_info['listOfPairRelationships']:
            match = full_pairs[(full_pairs.tumorId == info['tumorId']) &
                               (full_pairs.normalId == info['normalId'])]
            pair_infos.append({'normalId' : info['normalId'],
                               'tumorId' : info['tumorId'],
                               'pairId' : info['pairId'],
                               'tumorFinalBam' : {'bam': match.tumorBam.tolist()[0],
                                                  'bamIndex' : match.tumorBam.tolist()[0].replace('.bam', suffix)},
                               'normalFinalBam' : {'bam': match.normalBam.tolist()[0],
                                                   'bamIndex' : match.normalBam.tolist()[0].replace('.bam', suffix)}
                               })
            normal_sample_infos.append({'sampleId' : info['normalId'],
                                        'finalBam' : {'bam': match.normalBam.tolist()[0],
                                                      'bamIndex' : match.normalBam.tolist()[0].replace('.bam', suffix)}
                               })
        project_info = note_updates(key='pairInfos', new_value=pair_infos, project_info=project_info)
        project_info = note_updates(key='normalSampleBamInfos', new_value=normal_sample_infos, project_info=project_info)
    return project_info
        
    
def populate(args, custom_inputs):
    '''Update the dictionary of project, run and sample related metadata (project_info)
    where a new argument has been specified'''
    project_info = {}
    project_info['options'] = read(args['options'])
    # update wdl pipeline monitor version
    commit, tag, uniq_tag, branch = git_log(__file__)
    project_info['commit'] = commit
    project_info['tag'] = tag
    project_info['branch'] = branch
    project_info['uniq_tag'] = uniq_tag
    # update wdl pipeline version (if stored in a git repo)
    commit_wdl, tag_wdl, uniq_tag_wdl, branch_wdl = git_log(args['wdl_file'],
                                                            soft_fail=True)
    project_info['commit_wdl'] = commit_wdl
    project_info['tag_wdl'] = tag_wdl
    project_info['branch_wdl'] = branch_wdl
    project_info['uniq_tag_wdl'] = uniq_tag_wdl
    # update as needed
    verify_required(key='library', args=args, project_info=project_info)
    project_info = note_updates(key='library', new_value=args['library'], project_info=project_info)
    verify_required(key='genome', args=args, project_info=project_info)
    project_info = note_updates(key='genome', new_value=args['genome'], project_info=project_info)
    verify_required(key='project_name', args=args, project_info=project_info)
    project_info = note_updates(key='project_name', new_value=args['project_name'], project_info=project_info)
    if args['library'] == 'WGS':
        project_info['intervalList'] = 'default'
    else:
        verify_required(key='intervalList', args=args, project_info=project_info)
        project_info = note_updates(key='intervalList', new_value=args['intervalList'], project_info=project_info)
    if args['pairs_file']:
        pairs, full_pairs = load_pairs(args['pairs_file'])
        pair_info = []
        pair_info_relationships = []
        # update basic pairing information with tumor, normalId and pairId from pairs file
        for index, row in pairs.iterrows():
            # add pairRelationship
            current_pair_info_relationship = fill_pair_relationship(row)
            pair_info_relationships.append(current_pair_info_relationship)  
        project_info = note_updates(key='listOfPairRelationships', new_value=pair_info_relationships, project_info=project_info)
        # update BAM objects in project_info from custom inputs JSON (if pairInfos and normalSampleBamInfos  are in the custom json)
        if custom_inputs:
            for alt_project_info_file in custom_inputs:
                alt_project_info = read(alt_project_info_file)
                project_info = note_custom_updates(key='pairInfos', alt_project_info=alt_project_info, project_info=project_info)
                project_info = note_custom_updates(key='normalSampleBamInfos', alt_project_info=alt_project_info, project_info=project_info)
        # update BAM objects in project_info from sample sheet (if pairs file is a sample sheet)
        if not 'pairInfos' in project_info \
                    or not 'normalSampleBamInfos' in project_info:
            if 'tumorBam' in full_pairs.columns and 'normalBam' in full_pairs.columns:
                project_info = fill_in_pair_info(project_info, full_pairs)
        # update pairing objects in project_info
        pair_ids = list(set([info['pairId'] for info in project_info['listOfPairRelationships']]))
        project_info = note_updates(key='pairIds', new_value=pair_ids, project_info=project_info)
        normalIds = list(set([info['normalId'] for info in project_info['listOfPairRelationships']]))
        project_info = note_updates(key='normalIds', new_value=normalIds, project_info=project_info)
        tumorIds = list(set([info['tumorId'] for info in project_info['listOfPairRelationships']]))
        project_info = note_updates(key='tumorIds', new_value=tumorIds, project_info=project_info)
    # fill in list of samples
    try:
        if args['samples_file']:
            sample_ids = load_sample_ids(args['samples_file'])
        else:
            sample_ids = list(set(project_info['normalIds'] + project_info['tumorIds']))
    except KeyError:
        sample_ids = []
    project_info = note_updates(key='sampleIds', new_value=sample_ids, project_info=project_info)
    if custom_inputs:
        for alt_project_info_file in custom_inputs:
                alt_project_info = read(alt_project_info_file)
                project_info = note_custom_updates(key='normalSampleBamInfos', alt_project_info=alt_project_info, project_info=project_info)
                project_info = note_custom_updates(key='sampleInfos', alt_project_info=alt_project_info, project_info=project_info)
    return project_info

def load_json(file):
    with open(file) as input:
        inputs = json.load(input)
        return inputs

def write_json(file_out, object):
    with open(file_out, 'w') as input_info_file:
        json.dump(object, input_info_file, indent=4)

def write_labels_json(args, project_info, labels_files,
                 custom_labels_out):
    '''
    'analysis_project' : The name of the analysis project to associate the workflow with. Workflows and samples will be grouped under this project in Tugboat
    'sample_ids':  # comma separeted in the final output
    sample_analysis_group:  The name of the group of samples being analyzed. This is typically the name of the sample for a singleton group or a comma separated list of sample names for multi-sample groups, or some common identifying attribute for the samples like the participant id. Samples will be organized by sample analysis group name in Tugboat
    category: The category of the project, e.g. Production, Clinical, Research
    customer: The customer identifier associated with the project, typically first initial + last name. Some projects just use NYGC if there's no external customer
    reference_genome: The reference genome used for the samples
    workflow_version: The workflow version being run, e.g. 1.0.0
    workflow_release: The release tag of the workflow version, e.g. 20221213. SWENG will create multiple releases of a workflow version if we change an underlying dependency and we typically tag it with the date. This is not required if it's the same as workflow version
    '''
    custom_labels = {}
    if labels_files != None:
        for labels_file in labels_files:
            custom_labels_dict = load_json(labels_file)
            for key in custom_labels_dict:
                custom_labels[key] = custom_labels_dict[key]
    labels = {'analysis_project': project_info['project_name'],
              'category': 'CompBio',
              'customer': 'NYGC',
              'reference_genome': project_info['genome'],
              'workflow_version': project_info['tag_wdl'],
              'workflow_commit': project_info['commit_wdl'] ,
              'workflow_branch': project_info['branch_wdl'] ,
              'workflow_release': project_info['uniq_tag_wdl']}
    ids = []
    if 'sampleIds' in project_info:
        ids.append(','.join(project_info['sampleIds']))
    if 'pairIds' in project_info:
        ids.append(','.join(project_info['pairIds']))
        ids.append(','.join(project_info['tumorIds']))
        ids.append(','.join(project_info['normalIds']))
        labels['sample_analysis_group'] = ','.join(ids)
    for key in custom_labels:
        labels[key] = custom_labels[key]
    write_json(file_out=custom_labels_out, object=labels)

def write_wdl_json(args, project_info, custom_inputs,
                   custom_inputs_out,
                   local=False):
    parent_dir = os.path.abspath(os.path.dirname(__file__))
    if local:
        pipeline_input = parent_dir + '/../config/pipeline_references_local.json'
        interval_input = parent_dir + '/../config/interval_references_local.json'
        genome_input = parent_dir + '/../config/fasta_references_local.json'
    else:
        pipeline_input = parent_dir + '/../config/pipeline_references.json'
        interval_input = parent_dir + '/../config/interval_references.json'
        genome_input = parent_dir + '/../config/fasta_references.json'
    input = parse.Wdl(args['wdl_file'], 
                      genome_input=genome_input,
                      interval_input=interval_input,
                      pipeline_input=pipeline_input,
                      genome=args['genome'],
                      custom_inputs=custom_inputs,
                      validate=not args['skip_validate'],
                      project_info=project_info,
                      local=args['local'])
    file_out = custom_inputs_out
    write_json(file_out=file_out, object=input.final_inputs)
                
                
def test_schema(json_data):
    '''Test that project metadata fits schema'''
    parent_dir = os.path.abspath(os.path.dirname(__file__))
    schema_file = parent_dir + '/project_schema.json'
    with open(schema_file) as schema_input:
        schema = json.load(schema_input)
    try:
        validate(instance=json_data, schema=schema)
    except jsonschema.exceptions.ValidationError as err:
        log.error(' JSON data is invalid: ')
        help(err)
        log.error(err.message)
        return False
    return True
                

def get_args():
    '''Parse input flags
        Need to add optional task-specific input json
        
    If "production" or "external" inputs are true this triggers the skipping of internal-only files 
        (the steps that run with these also are written to skip tasks that localize these files 
        so any inability to read these files will not negatively affect the run
    '''
    parser = argparse.ArgumentParser(description='Program uses the WDL to determine which variables are required. '
                                                 'Creation of input JSON: \n'  
                                                 'The WDL is used to determine which variables are required. '
                                                 'Required or optional variables are defined from custom inputs JSON. '
                                                 'Any required variable not defined in the custom inputs JSON will be defined from the '
                                                 'reference JSONs in the config directory (as long as the variable names are identical).'
                                                 'The pairing/sample info CSVs (--pairs-file/--samples-file) are used to create pairRelationships and '
                                                 '(if columns named tumorBam and normalBam exist) map BAMs to pairs. '
                                                 'If "production" or "external" inputs are true then validation of NYGC internal-only files is skipped '
                                                 'The pipelines are written to skip tasks that localize these files if "production" or "external" are true '
                                                 'so any inability to read these files will not negatively affect the run.')
    parser.add_argument('--options',
                        help='Options json file (required)',
                            required=True
                            )
    parser.add_argument('--wdl-file',
                        help='WDL workflow. '
                        'An input JSON that matches this '
                        'WDL workflow will be created (required)',
                        required=True
                        )
    parser.add_argument('--file-out',
                        help='Output file for project info',
                        required=True
                        )
    parser.add_argument('--custom-inputs-out',
                        help='Output file for custom workflow inputs',
                        required=True
                        )
    parser.add_argument('--custom-labels-out',
                        help='Output file for workflow labels',
                        required=True
                        )
    parser.add_argument('--library',
                        help='Sequence library type. If not supplied '
                            'define library using --project-data',
                            choices=['WGS',
                                     'Exome'],
                            required=False
                            )
    parser.add_argument('--genome',
                        help='Genome key to use for pipeline. If not supplied '
                        'define genome using --project-data',
                        choices=['Human_GRCh38_full_analysis_set_plus_decoy_hla',
                                 'Human_GRCh38_tcga'],
                        required=False
                        )
    parser.add_argument('--project-name',
                        help='Project name associated with account. If not supplied '
                        'define genome using --project-data',
                        required=False
                        )
    parser.add_argument('--pairs-file',
                        help='CSV file with items that are required to have '
                        '"tumorId", "normalId", "pairId" as columns '
                        'Optionally, include "tumorBam", "normalBam" columns to create '
                        '"pairInfos" and "normalSampleBamInfos" automatically.'
                        ,
                        required=False
                        )
    parser.add_argument('--samples-file',
                        help='Not generally required. '
                        'If steps run only require sampleId and do not use pairing '
                        'information sample info can be populated with a CSV file. '
                        'The CSV file requires a columns named ["sampleId"]. If not supplied '
                        'define samples using --project-data',
                        required=False
                        )
    parser.add_argument('--interval-list',
                        help='File basename for interval list.'
                        'If not supplied the default (the SureSelect '
                        'interval list for your genome) will be used',
                        dest='intervalList',
                        choices=['SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla'],
                        required=False
                        )
    parser.add_argument('--labels',
                        help='Optional comma separated list of JSON files with custom labels. '
                        'By default analysis_project,  sample_analysis_group, sample_ids, category '
                        'customer, reference_genome, workflow_version, workflow_release '
                        'will be filled in using as much available information as possible. '
                        'If the labels json include any of these that value will override the auto-generated value.'
                        'By default the input will be added to the top-level workflow.',
                        required=False,
                        default=False
                        )
    parser.add_argument('--custom-inputs',
                        help='Optional comma separated list of JSON files with custom input variables. '
                        'The name of the variable in the input file must match the '
                        'name of the variable in the WDL workflow. '
                        'It is not required that the input specify the workflow. '
                        'By default the input will be added to the top-level workflow.',
                        required=False,
                        default=False
                        )
    parser.add_argument('--skip-validate',
                        help='Skip the step where input files are validated. '
                        'Otherwise all gs//: URIs will be checked to '
                        'see that a file exists. '
                        'Disable with caution.Cromwell will launch instances and run without checking. '
                        'Test a small pairs file to ensure all references exist and at least some sample input files '
                        'can be read by the current user. ',
                        required=False,
                        action='store_true'
                        )
    parser.add_argument('--local',
                        help='Use local path for all files and run on cluster.',
                        required=False,
                        action='store_true'
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

def main():
    args = get_args()
    if args['custom_inputs']:
        custom_inputs = args['custom_inputs'].split(',')
    else:
        custom_inputs = None
    if args['labels']:
        labels = args['labels'].split(',')
    else:
        labels = None
    # project, run and sample related metadata (project_info)
    project_info = populate(args, custom_inputs)
    passed = test_schema(project_info)
    if passed:
        file_out = args['file_out']
        # print project json
        write_json(file_out=file_out,
                   object=project_info)
        # create/print WDL inputs JSON
        write_wdl_json(args, project_info,
                       custom_inputs=custom_inputs,
                       custom_inputs_out=args['custom_inputs_out'],
                       local=args['local'])
        write_labels_json(args, project_info,
                          labels_files=labels,
                          custom_labels_out=args['custom_labels_out'])

if __name__ == "__main__":
    main()       

