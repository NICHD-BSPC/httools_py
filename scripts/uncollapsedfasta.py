#!/usr/bin/env python


import os
from pathlib import Path
import ruamel.yaml
import argparse
from Bio import SeqIO
import pandas as pd
import re

ap = argparse.ArgumentParser()
ap.add_argument('--config', help='config.yaml file', default='config/config.yaml')
ap.add_argument('--fasta', help='screen fasta file', default=False)
ap.add_argument('--ids', help='screen IDs file', default=False)
ap.add_argument('--integration', help='integration sites file', default=False)
ap.add_argument('--sample', help='sample', default=False)
args = ap.parse_args()

yaml = ruamel.yaml.YAML()
config = yaml.load(open(args.config))
samples = [args.sample] if args.sample else config['sample'].keys()
legacy_mode = config['legacy_mode']

logdir = os.path.dirname(args.config) + '/data/' + config['name'] + '/logs'
Path(logdir).mkdir(parents=True, exist_ok=True)
mydir = os.path.dirname(args.config) + '/data/' + config['name'] + '/uncollapsed_fastas'
Path(mydir).mkdir(parents=True, exist_ok=True)



class SampleConfig:
    """Make sample-specific config"""

    def __init__(self, config, sample, args):
        self.sample = sample
        self.samplevs = sample + '.' + config['genomevs'][config['genome']]
        self.fasta = self.get_fastafn(self.sample, config, args)
        self.ids = self.get_ids(self.samplevs, config, args)
        self.integ = self.get_integ(self.samplevs, config, args)

    def get_fastafn(self, sample, config, args):
        fastanm = args.fasta if args.fasta else os.path.dirname(args.config) + '/data/' + config['name'] + '/fastqscreen/screen_' + sample + '.fa'
        return fastanm

    def get_ids(self, sample, config, args):
        ids = args.ids if args.ids else os.path.dirname(args.config) + '/data/' + config['name'] + '/filblast/id_mappings_' + sample + '.txt'
        return ids

    def get_integ(self, sample, config, args):
        integ = args.integration if args.integration else os.path.dirname(args.config) + '/data/' + config['name'] + '/location/true_integration_' + sample + '.txt'
        return integ

def get_id_dict(samplecfg):
    # reads mappings of ssp_id vs. sequence id,
    # subsets to ssp_ids in integration file
    # returns dictionnary sequence ids with number of duplicates
    ids = pd.read_csv(samplecfg.ids, sep='\t', index_col='Id')
    integ = pd.read_csv(samplecfg.integ, sep='\t')['Id'].tolist()
    ids = ids.loc[integ,'id'].tolist()

    def split_val(val, subpattern, splitpattern):
        val = re.sub(subpattern, '', val)
        val = re.split(splitpattern, val)
        return(val)

    splitted = [split_val(x, "\[|'|\]| ", ',') for x in ids]
    flat = [x for v in splitted for x in v]
    ids_dict = {x: int(split_val(x, 'dupl', '_')[1]) for x in flat}
    return ids_dict

def write_fa(samplecfg, ids_dict):
    basename = os.path.splitext(os.path.split(samplecfg.integ)[1])[0]
    fn = mydir + '/uncollapsed_' + basename + '.fa'
    with open(fn, 'w') as fnfile:
        for seqrec in SeqIO.parse(samplecfg.fasta, 'fasta'):
            n = ids_dict[seqrec.id] if (seqrec.id in ids_dict) else 0
            for x in range(0, n):
                SeqIO.write(seqrec, fnfile, 'fasta')


def main():
    for sample in samples:
        samplecfg = SampleConfig(config, sample, args)
        ids_dict = get_id_dict(samplecfg)
        write_fa(samplecfg, ids_dict)
main()


