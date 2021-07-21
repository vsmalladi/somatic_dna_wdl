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
    return pairs[['tumor', 'normal', 'pairId']], pairs

def load_sample_ids(file):
    '''Return a deduplicated list of sample ids'''
    sample_ids = pd.read_csv(file)
    assert 'sampleId' in sample_ids.columns , 'Error: can not find "sampleId" columns in samples file: ' + ' '.join(sample_ids.columns)
    return sample_ids.sampleId.unique().tolist()

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
    if sample_id == 'CA-0103N-D-W':
        sample_info = {"sampleId" : sample_id}
        fastq_pair = {'fastqR1': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L001_001.R1.fastq.gz",
                       'fastqR2': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L001_001.R2.fastq.gz",
                       'readgroupId': 'CA-0110N-D-W_AAGGTACA_HY7KNCCXX_L001',
                       'flowcell': 'HY7KNCCXX',
                       'lane': 'L001',
                       'sampleId' : sample_id,
                       'rgpu': 'HY7KNCCXX.L002.AAGGTACA',
                       'barcode': 'AAGGTACA'}
        fastq_pair2 = {'fastqR1': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L002_001.R1.fastq.gz",
                       'fastqR2': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103N-D-W_CCGAAGTA_HYVLCCCXX_L002_001.R2.fastq.gz",
                       'readgroupId': 'CA-0110N-D-W_AAGGTACA_HY7KNCCXX_L001',
                       'flowcell': 'HY7KNCCXX',
                       'lane': 'L002',
                       'sampleId' : sample_id,
                       'rgpu': 'HY7KNCCXX.L002.AAGGTACA',
                       'barcode': 'AAGGTACA'}
    else:
        sample_info = {"sampleId" : sample_id}
        fastq_pair = {'fastqR1': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103T-D-W_AGATCGCA_HYVLCCCXX_L005_001.R1.fastq.gz",
                       'fastqR2': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103T-D-W_AGATCGCA_HYVLCCCXX_L005_001.R2.fastq.gz",
                       'readgroupId': 'CA-0103T-D-W_AGATCGCA_HY7KNCCXX_L005',
                       'flowcell': 'HY7KNCCXX',
                       'lane': 'L005',
                       'sampleId' : sample_id,
                       'rgpu': 'HYVLCCCXX.L005.AGATCGCA',
                       'barcode': 'AGATCGCA'}
        fastq_pair2 = {'fastqR1': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103T-D-W_AGATCGCA_HYVLCCCXX_L006_001.R1.fastq.gz",
                       'fastqR2': "gs://nygc-comp-s-82ed-input/WGS/CA-103/fastq/CA-0103T-D-W_AGATCGCA_HYVLCCCXX_L006_001.R2.fastq.gz",
                       'readgroupId': 'CA-0103T-D-W_AGATCGCA_HY7KNCCXX_L006',
                       'flowcell': 'HY7KNCCXX',
                       'lane': 'L006',
                       'sampleId' : sample_id,
                       'rgpu': 'HYVLCCCXX.L006.AGATCGCA',
                       'barcode': 'AGATCGCA'}
    sample_info['listOfFastqPairs']  = [fastq_pair, fastq_pair2]
    return sample_info

def fill_pair_relationship(row):
    '''Add pair and sample level info to object
    PairRelationship '''
    pair_relationship = {"pairId" : row.pairId}
    pair_relationship['normal'] = row['normal']
    pair_relationship['tumor'] = row['tumor']
    return pair_relationship      

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
    ''' Add/replace value in project_info with value from flag'''
    if args_key:
        if key in project_info:
            if not project_info[key] == args_key:
                log.warning('Note that the value for ' + key + ' differs in the --project-info file:\n' + str(project_info[key]) + '\n from the value indicated with the flag:\n' + str(args_key) + '\nThe value from the flag supersedes the file value')
        project_info[key] = args_key
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

def fill_in_pair_info(project_info, full_pairs, suffix='.bam.bai'):
    pair_infos = []
    normal_sample_infos = []
    if 'tumor_bam' in full_pairs and 'normal_bam' in full_pairs:
        for info in project_info['listOfPairRelationships']:
            match = full_pairs[(full_pairs.tumor == info['tumor']) &
                               (full_pairs.normal == info['normal'])]
            pair_infos.append({'normal' : info['normal'],
                               'tumor' : info['tumor'],
                               'pairId' : info['pairId'],
                               'tumorFinalBam' : {'bam': match.tumor_bam.tolist()[0],
                                                  'bamIndex' : match.tumor_bam.tolist()[0].replace('.bam', suffix)},
                               'normalFinalBam' : {'bam': match.normal_bam.tolist()[0],
                                                   'bamIndex' : match.normal_bam.tolist()[0].replace('.bam', suffix)}
                               })
            normal_sample_infos.append({'sampleId' : info['normal'],
                                        'finalBam' : {'bam': match.normal_bam.tolist()[0],
                                                      'bamIndex' : match.normal_bam.tolist()[0].replace('.bam', suffix)}
                               })
    else:
        for info in project_info['listOfPairRelationships']:
            pair_infos.append({'normal' : info['normal'],
                               'tumor' : info['tumor'],
                               'pairId' : info['pairId']})
            normal_sample_infos.append({'sampleId' : info['normal']})
    project_info = note_updates(key='pairInfos', args_key=pair_infos, project_info=project_info)
    project_info = note_updates(key='normalSampleBamInfos', args_key=normal_sample_infos, project_info=project_info)
    return project_info

def fill_in_sample_info(project_info):
    sample_infos = []
    normal_sample_infos = []
    for sample_id in project_info['sampleIds']:
        sample_infos.append({'sampleId' : sample_id})
        if sample_id in project_info['normals']:
            normal_sample_infos.append({'sampleId' : sample_id})
    project_info = note_updates(key='sampleInfos', args_key=sample_infos, project_info=project_info)
    project_info = note_updates(key='normalSampleInfos', args_key=normal_sample_infos, project_info=project_info)
    return project_info
        
    
def repopulate(args):
    '''Ammend the dictionary of project and sample related metadata 
    where a new argument has been specified'''
    if args['project_data']:
        project_info = read(args['project_data'])
    else:
        project_info = {}
    project_info['options'] = read(args['options'])
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
    if args['pairs_file'] and not args['project_data']:
        pairs, full_pairs = load_pairs(args['pairs_file'])
        pair_info = []
        pair_info_relationships = []
        for index, row in pairs.iterrows():
            # add test pairInfo
            if args['test_data']:
                current_pair_info = fill_pair(row)
                pair_info.append(current_pair_info)
            # add pairRelationship
            current_pair_info_relationship = fill_pair_relationship(row)
            pair_info_relationships.append(current_pair_info_relationship)   
        # pair-only         
        project_info = note_updates(key='listOfPairRelationships', args_key=pair_info_relationships, project_info=project_info)
        if args['test_data']:
            project_info = note_updates(key='pairInfos', args_key=pair_info, project_info=project_info)
            project_info = note_updates(key='normalSampleBamInfos', args_key=pair_info, project_info=project_info)
        else:
            if args['custom_inputs']:
                for alt_project_info_file in args['custom_inputs']:
                    alt_project_info = read(alt_project_info_file)
                    project_info = note_custom_updates(key='pairInfos', alt_project_info=alt_project_info, project_info=project_info)
                    project_info = note_custom_updates(key='normalSampleBamInfos', alt_project_info=alt_project_info, project_info=project_info)
            if not 'pairInfos' in project_info:
                project_info = fill_in_pair_info(project_info, full_pairs)
            assert 'pairInfos' in project_info, 'pairInfos needed but no entry found for key in --custom-inputs file\n'
        pair_ids = list(set([info['pairId'] for info in project_info['listOfPairRelationships']]))
        project_info = note_updates(key='pairId', args_key=pair_ids, project_info=project_info)
        
        normals = list(set([info['normal'] for info in project_info['listOfPairRelationships']]))
        project_info = note_updates(key='normals', args_key=normals, project_info=project_info)
        
        tumors = list(set([info['tumor'] for info in project_info['listOfPairRelationships']]))
        project_info = note_updates(key='tumors', args_key=tumors, project_info=project_info)
    # fill in list of samples
    if args['samples_file'] and not args['project_data']:
        sample_ids = load_sample_ids(args['samples_file'])
    else:
        sample_ids = list(set(project_info['normals'] + project_info['tumors']))
    project_info = note_updates(key='sampleIds', args_key=sample_ids, project_info=project_info)
    if args['test_data']:
        sample_info = []
        normal_sample_info = []
        for sample_id in sample_ids:
            # add mock fastq files for now (will be replaced later by custom input if made available)
            if sample_id in project_info['normals']:
                normal_sample_info.append(fill_sample(sample_id))
            current_sample_id_info = fill_sample(sample_id)
            sample_info.append(current_sample_id_info)
        project_info = note_updates(key='normalSampleInfos', args_key=normal_sample_info, project_info=project_info)
        project_info = note_updates(key='sampleInfos', args_key=sample_info, project_info=project_info)
    else:
        if args['custom_inputs']:
            for alt_project_info_file in args['custom_inputs']:
                    alt_project_info = read(alt_project_info_file)
                    project_info = note_custom_updates(key='normalSampleInfos', alt_project_info=alt_project_info, project_info=project_info)
                    project_info = note_custom_updates(key='sampleInfos', alt_project_info=alt_project_info, project_info=project_info)
                    if not 'sampleInfos' in project_info:
                        project_info = fill_in_sample_info(project_info)
                    assert 'sampleInfos' in project_info, 'sampleInfos needed but no entry found for key in --custom-inputs file\n'
    return project_info

def write_wdl_json(args, project_info, project_info_file):
    parent_dir = os.path.abspath(os.path.dirname(__file__))
    pipeline_input = parent_dir + '/../config/pipeline_references.json'
    interval_input = parent_dir + '/../config/interval_references.json'
    genome_input = parent_dir + '/../config/fasta_references.json'
    input = parse.Wdl(args['wdl_file'], 
                      genome_input=genome_input,
                      interval_input=interval_input,
                      pipeline_input=pipeline_input,
                      genome=args['genome'],
                      read_length=args['read_length'],
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
        help(err)
        log.error(err.message)
        return False
    return True
                

def get_args():
    '''Parse input flags
        Need to add optional task-specific input json
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--test-data',
                        help='Add FASTQ and BAM objects to run as test data',
                        action='store_true'
                       )
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
    parser.add_argument('--interval-list',
                        help='File basename for interval list.'
                        'If not supplied the default (the SureSelect '
                        'interval list for your genome) will be used',
                        dest='intervalList',
                        choices=['SureSelect_V6plusCOSMIC.target.GRCh38_full_analysis_set_plus_decoy_hla'],
                        required=False
                        )
    parser.add_argument('--project-data',
                        help='Optional JSON file with project pairing, sample, '
                        'genome build, library and interval list information',
                        required=False
                        )
    parser.add_argument('--read-length',
                        help='Required only for steps like BiqSeq2 that use read_length as input',
                        required=False,
                        dest='read_length'
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
                        required=False,
                        action='store_true'
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

