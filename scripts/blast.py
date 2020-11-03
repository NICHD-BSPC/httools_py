#!/usr/bin/env python


import os
from pathlib import Path
import ruamel.yaml
import argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline


ap = argparse.ArgumentParser()
ap.add_argument('--config', help='config.yaml file', default='config/config.yaml')
ap.add_argument('--fasta', help='screen fasta file', default=False)
ap.add_argument('--genome', help='genome database fasta file', default=False)
ap.add_argument('--sample', help='sample', default=False)
args = ap.parse_args()

yaml = ruamel.yaml.YAML()
config = yaml.load(open(args.config))
samples = [args.sample] if args.sample else config['sample'].keys()
legacy_mode = config['legacy_mode']

logdir = 'data/' + config['name'] + '/logs'
Path(logdir).mkdir(parents=True, exist_ok=True)
mydir = 'data/' + config['name'] + '/blast'
Path(mydir).mkdir(parents=True, exist_ok=True)


class SampleConfig:
    """Make sample-specific config"""

    def __init__(self, config, sample, args):
        self.sample = sample
        self.genome = self.get_genome(config, args)
        self.view = config['blastview']
        self.evalue = config['blastevalue']
        self.fasta = self.get_fastafn(sample, args)
        self.suffix = self.get_suffix(config, args)

    def get_fastafn(self, sample, args):
        fastanm = args.fasta if args.fasta else 'data/' + config['name'] + '/fastqscreen/screen_' + sample + '.fa'
        return fastanm

    def get_genome(self, config, args):
        genomefa = args.genome if args.genome else os.path.join(os.path.dirname(__file__), '..', config['genomedb'][config['genome']])
        return genomefa

    def get_suffix(self, config, args):
        suffix = os.path.splitext(os.path.basename(args.genome))[0] if args.genome else config['genomevs'][config['genome']]
        return suffix

def blast_func(samplecfg):
    blastlog = ['BLAST']
    # create index if does not exist
    if not (os.path.exists(samplecfg.genome + '.nin') and
            os.path.exists(samplecfg.genome + '.nhr') and
            os.path.exists(samplecfg.genome + '.nsq')):
        makedb = NcbimakeblastdbCommandline(cmd='makeblastdb', dbtype='nucl', input_file=samplecfg.genome)
        stdout, stderr = makedb()
        blastlog.append('\n\nNcbimakeblastdb\n\n')
        blastlog.append(stdout)
        blastlog.append(stderr)
    # blast
    sampleblast = NcbiblastnCommandline(task = 'blastn',
                                        query = samplecfg.fasta,
                                        db = samplecfg.genome,
                                        outfmt = samplecfg.view,
                                        evalue = samplecfg.evalue,
                                        out = mydir + '/blast_' + samplecfg.sample + '.' + samplecfg.suffix + '.txt')
    stdout, stderr = sampleblast()
    blastlog.append('\n\nNcbiblastn\n\n')
    blastlog.append(stdout)
    blastlog.append(stderr)
    return blastlog

def main():
    for sample in samples:
        samplecfg = SampleConfig(config, sample, args)
        logfn = logdir + '/blast_' + config['name'] + '_' + sample + '.' + samplecfg.suffix + '.log.txt'
        with open(logfn, 'w+') as logfile:
            blastlog = blast_func(samplecfg) if os.path.exists(samplecfg.fasta) else '\nfile does not exist: ' + samplecfg.fasta + '\n'
            [logfile.write(line) for line in blastlog]
main()


