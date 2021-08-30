'''
Created on Feb 8, 2021

@author: jshelton
'''
import os
import sys
import json
import argparse
import pandas as pd
import subprocess
import logging as log
import Colorer
from collections.abc import Iterable
import copy
import glob
import meta

log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)

        
class Gather():
    '''Load rundata from file and confirm project data schema '''
    def __init__(self, uris_file, options_file,
                 project_id,
                 library,
                 genome,
                 pairs_file=False, 
                 samples_file=False):
        self.library = library
        self.genome = genome
        self.project_id = project_id
        self.uris = pd.read_csv(uris_file, header=None, names=['uri'])
        self.options = self.read(options_file)
        if pairs_file:
            self.pairs, self.full_pairs = meta.load_pairs(pairs_file)
            self.sample_ids = False
        elif samples_file:
            self.sample_ids = meta.load_sample_ids(samples_file)
        self.compile()
        self.write(file_out=project_id + '.RunInfos.json')
        
    def read(self, file):
        with open(file) as project_info_file:
            project_info = json.load(project_info_file)
            return project_info
        
    def get_uuid(self, uri):
        uuid = uri.split('/')[4]
        return uuid
    
    def get_id(self, uri):
        '''assumes no periods in sample id'''
        id = uri.split('/')[-1].split('.')[0]
        return id
    
    def get_from_id(self, id, field='pairId'):
        '''Match probable id to id name options'''
        match = self.pairs[(self.pairs.tumor == id) | (self.pairs.normal == id) | (self.pairs.pairId == id)]
        return match[field].tolist()[0]
    
    def get_project_data(self, row):
        project_data = {}
        if self.sample_ids:
            project_data['listOfPairRelationships'] = {}
            project_data['sampleIds'] = [row.id]
            project_data['options'] = self.options
        else:
            project_data['tumors'] = [self.get_from_id(row.id, field='tumor')]
            project_data['normals'] = [self.get_from_id(row.id, field='normal')]
            project_data['pairId'] = [self.get_from_id(row.id, field='pairId')]
            project_data['listOfPairRelationships'] = [{'tumor': project_data['tumors'][0],
                                                       'normal' : project_data['normals'][0],
                                                       'pairId' : project_data['pairId'][0]}]
            project_data['sampleIds'] = [project_data['tumors'][0],
                                         project_data['normals'][0]]
            project_data['options'] = self.options
            project_data['library'] = self.library
            project_data['genome'] = self.genome
            project_data['project'] = self.project_id
            return project_data
        
        
    def compile(self):
        self.uris['workflow_uuid'] = self.uris.apply(lambda row: self.get_uuid(row.uri), axis=1)
        self.uris['id'] = self.uris.apply(lambda row: self.get_id(row.uri), axis=1)
        self.uris['project_data'] = self.uris.apply(lambda row: self.get_project_data(row), axis=1)
        self.uris = self.uris[['workflow_uuid', 'project_data']]
        run_info_json = self.uris.to_json(orient='records')
        self.run_info_json = json.loads(run_info_json)
        
    def write(self, file_out):
        with open(file_out, 'w') as run_infos_file:
            json.dump(self.run_info_json, run_infos_file, indent=4)

def main():
    args = get_args()
    run_data = Gather(uris_file=args['uris'],
                      project_id=args['project_id'],
                      options_file=args['options'],
                      pairs_file=args['pairs_file'],
                      samples_file=args['samples_file'],
                      library=args['library'],
                      genome=args['genome'])


 
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--uris',
                        help='File of URIs that are of interest',
                        required=True
                        )
    parser.add_argument('--options',
                        help='Options JSON file.',
                        required=True
                        )
    parser.add_argument('--project-id',
                        help='Overall project name (research project rather than GCP project).',
                        required=True
                        )
    parser.add_argument('--pairs-file',
                        help='JSON file with items that are required to have '
                        '"tumor", "normal" sample_ids defined. If not supplied '
                        'define pairing using --project-data',
                        required=False
                        )
    parser.add_argument('--library',
                        help='Sequencing library',
                        required=True,
                        choices=['WGS', 'Exome', 'RNA']
                        )
    parser.add_argument('--genome',
                        help='Genome key value',
                        required=True,
                        choices=["Human_GRCh38_full_analysis_set_plus_decoy_hla",
                                 "Human_GRCh38_tcga"]
                        )
    parser.add_argument('--samples-file',
                        help='Not generally required. '
                        'If steps run only require sample_id and do not use pairing '
                        'information sample info can be populated with a CSV file. '
                        'The CSV file requires a columns named ["sampleId"]. If not supplied '
                        'define samples using --project-data',
                        required=False
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

if __name__ == "__main__":
    main()       
