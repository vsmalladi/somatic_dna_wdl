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

def git_log(program):
    '''
        Function returns the most recent git commit and tag.
    '''
    import subprocess
    potential_git_repo = os.path.abspath(os.path.dirname(program))
    while True:
        if not os.path.isdir(potential_git_repo):
            log.error(')No parent directories of the program contain a .git directory :' + program )
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
    assert 'tumor' in pairs.columns and 'normal' in pairs.columns, 'Error: can not find tumor and normal columns in pairs file: ' + ' '.join(pairs.columns)
    pairs['pairId'] = pairs.apply(lambda row: row.tumor + '--' + row.normal, axis=1)
    return pairs[['tumor', 'normal', 'pairId']]

def load_bam(row, pair_info, kind):
    '''temp version to add bam file replace once GCP file manifest format is known'''
    if kind == 'tumor':
        bam = {"bam" : "gs://nygc-comp-s-82ed-input/WGS/CA-103/CA-0103T-D-W.25x.bam",
               "bamIndex" : "gs://nygc-comp-s-82ed-input/WGS/CA-103/CA-0103T-D-W.25x.bai"}
    if kind == 'normal':
        bam = {"bam" : "gs://nygc-comp-s-82ed-input/WGS/CA-103/CA-0103N-D-W.15x.bam",
               "bamIndex" : "gs://nygc-comp-s-82ed-input/WGS/CA-103/CA-0103N-D-W.15x.bai"}
    return bam
 
def load_sample_info(sample_id):
    '''temp version to add FASTQ file replace once GCP file manifest format is known'''
    sample_info = {"sampleId" : sample_id}
    fastq_pair = {'fastqR1': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L001_001.R1.fastq.gz",
                   'fastqR2': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L001_001.R2.fastq.gz",
                   'readgroupId': 'CA-0110N-D-W_AAGGTACA_HY7KNCCXX_L001',
                   'flowcell': 'HY7KNCCXX',
                   'lane': 'L001',
                   'barcode': 'AAGGTACA'}
    fastq_pair2 = {'fastqR1': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L002_001.R1.fastq.gz",
                   'fastqR2': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L002_001.R2.fastq.gz",
                   'readgroupId': 'CA-0110N-D-W_AAGGTACA_HY7KNCCXX_L001',
                   'flowcell': 'HY7KNCCXX',
                   'lane': 'L002',
                   'barcode': 'AAGGTACA'}
    sample_info['listOfFastqPairs']  = [fastq_pair, fastq_pair2]
    return sample_info
        

def fill_pair(row):
    '''Add pair and sample level info to object'''
    pair_info = {"pairId" : row.pairId}
    if load_bam(row, pair_info, kind='tumor') and load_bam(row, pair_info, kind='normal'):
        pair_info['tumorFinalBam'] = load_bam(row, pair_info, kind='tumor')
        pair_info['normalFinalBam'] = load_bam(row, pair_info, kind='normal')
    pair_info['normal'] = row['normal']
    pair_info['tumor'] = row['tumor']
    return pair_info

def fill_sample(sample_id):
    '''Add sample level info to object'''
    sample_info = load_sample_info(sample_id)
    return sample_info

def note_updates(key, args_key, project_info):
    if args_key:
        if key in project_info:
            if not project_info[key] == args_key:
                log.warning('Note that the value for ' + key + ' differs in the --project-info file:\n' + str(project_info[key]) + '\n from the value indicated with the flag:\n' + str(args_key) + '\nThe value from the flag supersedes the file value')
        project_info[key] = args_key
    return project_info

def verify_required(key, args, project_info):
    '''Verify a required key'''
    in_project_info = key in project_info and project_info[key]
    in_args = key in args and args[key]
    if not in_project_info and not in_args:
        log.error('Note that the value for ' + key + ' is required to be specified in flags or in the  --project-info file')
        sys.exit(1)
    return True
    
def repopulate(args):
    '''Ammend the dictionary of project and sample related metadata 
    where a new argument has been specified'''
    if args['project_data']:
        project_info = read(args['project_data'])
    else:
        project_info = {}
    # update wdl pipeline version
    commit, tag, uniq_tag, branch = git_log(__file__)
    project_info['commit'] = commit
    project_info['tag'] = tag
    project_info['branch'] = branch
    project_info['uniq_tag'] = uniq_tag
    # update as needed
    verify_required(key='library', args=args, project_info=project_info)
    project_info = note_updates(key='library', args_key=args['library'], project_info=project_info)
    verify_required(key='genome', args=args, project_info=project_info)
    project_info = note_updates(key='genome', args_key=args['genome'], project_info=project_info)
    verify_required(key='project', args=args, project_info=project_info)
    project_info = note_updates(key='project', args_key=args['project'], project_info=project_info)
    if args['library'] == 'WGS':
        project_info['intervalList'] = 'default'
    else:
        verify_required(key='intervalList', args=args, project_info=project_info)
        project_info = note_updates(key='intervalList', args_key=args['intervalList'], project_info=project_info)
    if args['pairs_file']:
        pairs = load_pairs(args['pairs_file'])
        pair_info = []
        for index, row in pairs.iterrows():
            current_pair_info = fill_pair(row)
            pair_info.append(current_pair_info)
    project_info = note_updates(key='pairInfos', args_key=pair_info, project_info=project_info)
    
    pair_ids = list(set([info['pairId'] for info in project_info['pairInfos']]))
    project_info = note_updates(key='pairId', args_key=pair_ids, project_info=project_info)
    
    normals = list(set([info['normal'] for info in project_info['pairInfos']]))
    project_info = note_updates(key='normals', args_key=normals, project_info=project_info)
    
    tumors = list(set([info['tumor'] for info in project_info['pairInfos']]))
    project_info = note_updates(key='tumors', args_key=tumors, project_info=project_info)
    sample_ids = list(set(project_info['normals'] + project_info['tumors']))
    project_info = note_updates(key='sampleIds', args_key=sample_ids, project_info=project_info)
    sample_info = []
    for sample_id in sample_ids:
        current_sample_id_info = fill_sample(sample_id)
        sample_info.append(current_sample_id_info)
    project_info = note_updates(key='sampleInfos', args_key=sample_info, project_info=project_info)
    return project_info

def write_wdl_json(args, project_info, project_info_file):
    parent_dir = os.path.abspath(os.path.dirname(__file__))
    interval_input = parent_dir + '/../config/interval_references.json'
    genome_input = parent_dir + '/../config/fasta_references.json'
    input = parse.Wdl(args['wdl_file'], 
                      genome_input=genome_input, 
                      interval_input=interval_input, 
                      genome=args['genome'],
                      custom_inputs=args['custom_inputs'],
                      validate=not args['skip_validate'],
                      project_info_file=project_info_file)
    file_out = os.path.basename(args['wdl_file']).replace('.wdl', '') + 'Input.json'
    with open(file_out, 'w') as input_info_file:
        json.dump(input.inputs, input_info_file, indent=4)
            
            
def write_json(args, project_info, file_out):
    if not args['project_data']:
        with open(file_out, 'w') as project_info_file:
            json.dump(project_info, project_info_file, indent=4)
                
                
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
        log.error(err)
        return False
    return True
                

def get_args():
    '''Parse input flags
        Need to add optional task-specific input json
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--library',
                        help='Sequence library type. If not supplied '
                            'define library using --project-data',
                            choices=['WGS',
                                     'Exome'],
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
    parser.add_argument('--genome',
                        help='Genome key to use for pipeline. If not supplied '
                        'define genome using --project-data',
                        choices=['Human_GRCh38_full_analysis_set_plus_decoy_hla'],
                        required=False
                        )
    parser.add_argument('--project',
                        help='Project name associated with account. If not supplied '
                        'define genome using --project-data',
                        required=False
                        )
    parser.add_argument('--pairs-file',
                        help='JSON file with items that are required to have '
                        '"tumor", "normal" sample_ids defined. If not supplied '
                        'define pairing using --project-data',
                        required=False
                        )
    parser.add_argument('--samples-file',
                        help='Not generally required. '
                        'If steps run only require sample_id and do not use pairing '
                        'information sample info can be populated with a CSV file. '
                        'The CSV file requires a columns named ["sampleId"]. If not supplied '
                        'define samples using --project-data',
                        required=False
                        )
    parser.add_argument('--wdl-file',
                        help='WDL workflow. '
                        'To output an input JSON '
                        'that matches a WDL workflow parse the workflow file in as a flag.',
                        required=False
                        )
    parser.add_argument('--project-data',
                        help='Optional JSON file with project pairing, sample, '
                        'genome build, library and interval list information',
                        required=False
                        )
    parser.add_argument('--custom-inputs',
                        help='Optional JSON file with custom input variables. '
                        'The name of the variable in the input file must match the '
                        'name of the variable in the WDL workflow. '
                        'It is not required that the input specify the workflow. '
                        'By default the input will be added to the top-level workflow.',
                        nargs='*',
                        required=False
                        )
    parser.add_argument('--skip-validate',
                        help='Skip the step where input files are validated. '
                        'Otherwise all gs//: URIs will be checked to '
                        'see that a file exists. '
                        'Disable with caution.Cromwell will launch instances and run without checking. '
                        'Test a small pairs file to ensure all references exist and at least some sample input files '
                        'can be read by the current user. ',
                        required=False
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

def main():
    args = get_args()
    project_info = repopulate(args)
    passed = test_schema(project_info)
    if passed:
        file_out = args['project'].replace(' ', '_') + '_projectInfo.json'
        write_json(args, project_info,
                   file_out=file_out)
        if args['wdl_file']:
            write_wdl_json(args, project_info, 
                           project_info_file=file_out)

if __name__ == "__main__":
    main()       

