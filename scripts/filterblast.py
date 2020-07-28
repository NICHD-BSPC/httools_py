#!/usr/bin/env python


import os
from pathlib import Path
import ruamel.yaml
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO


ap = argparse.ArgumentParser()
ap.add_argument('--config', help='config.yaml file', default='../test/config/config_SN_legacy.yaml')
ap.add_argument('--blast', help='blast file', default=False)
ap.add_argument('--screenSN', help='screen SN file', default=False)
ap.add_argument('--sample', help='sample', default=False)
args = ap.parse_args()

yaml = ruamel.yaml.YAML()
config = yaml.load(open(args.config))
samples = [args.sample] if args.sample else config['sample'].keys()
legacy_mode = config['legacy_mode']

logdir = os.path.dirname(args.config) + '/data/' + config['name'] + '/logs'
Path(logdir).mkdir(parents=True, exist_ok=True)
mydir = os.path.dirname(args.config) + '/data/' + config['name'] + '/filblast'
Path(mydir).mkdir(parents=True, exist_ok=True)


class SampleConfig:
    """Make sample-specific config"""
    def __init__(self, config, sample, args):
        self.sample = sample
        self.lib_design = config['sample'][sample]['lib_design']
        self.max_score_diff = config['max_score_diff']
        self.error = config['blastevalue']
        self.tsd = config['tsd'][config['sample'][sample]['integrase']]
        self.SN = 1 if isinstance(config['sample'][sample]['SN_position'], int) and isinstance(config['sample'][sample]['SN_length'], int) else 0
        self.blastfn = self.get_blastfn(sample, config, args)
        self.SNdict = self.get_SNdict(sample, config, args) if self.SN else 0

    def get_blastfn(self, sample, config, args):
        blastfn = args.blast if args.blast else os.path.dirname(args.config) + '/data/' + config['name'] + '/blast/blast_' + \
            sample + '.' + config['genomevs'][config['genome']] + '.txt'
        return blastfn

    def get_SNdict(self, sample, config, args):
        SNfn = args.screenSN if args.screenSN else os.path.dirname(args.config) + '/data/' + config['name'] + '/fastqscreen/screen_' + sample + '_SN.txt'
        SNdict = {rec.id : rec.seq for rec in SeqIO.parse(SNfn, "fasta")}
        return SNdict


class BlastRow:
    """Screen blast result line"""
    def __init__(self, blast, samplecfg):
        self.id = blast['id']
        self.dupl = int(self.id.split('_')[1].replace('dupl', ''))
        self.indpt = 0 # placeholder
        self.chr = blast['chr']
        self.strand = self.get_pos(blast, samplecfg)[1]
        self.position = self.get_pos(blast, samplecfg)[0]
        self.Id = 'SSP_' + self.chr + '_' + str(self.position) + '_' + self.strand
        self.SNs = self.get_SNs(str(samplecfg.SNdict[self.id])) if samplecfg.SN else 0

    def get_pos(self, blast, samplecfg):
        if samplecfg.lib_design == 'U5':
            strand = '+' if blast['sstart'] <= blast['send'] else '-'
            pos = blast['sstart'] if blast['sstart'] <= blast['send'] else blast['sstart'] - samplecfg.tsd + 1
        elif samplecfg.lib_design == 'U3':
            strand = '-' if blast['sstart'] <= blast['send'] else '+'
            pos = blast['sstart'] if blast['sstart'] <= blast['send'] else blast['sstart'] - samplecfg.tsd + 1
        return pos, strand

    def get_SNs(self, SNstring):
        SNstring = SNstring.replace('[','').replace(']','').replace("'",'')
        SNs = SNstring.split(sep='-')
        return SNs

def check_lib_design(samplecfg):
    if samplecfg.lib_design not in ['U5', 'U3']:
        raise ValueError('lib_design can only be U5 or U3')

def is_uniquelymapped(blastgroup, samplecfg, uniquedict, multidict):
    """Determine uniquely or multi-mapped"""
    # sort by id, bitscore then evalue
    blastgroup = blastgroup.reset_index()
    if len(blastgroup.index) == 1:
        # add to uniquedf
        uniquedict.append(blastgroup.iloc[0,:].to_dict())
    elif ((abs(blastgroup.loc[0, 'bitscore'] - blastgroup.loc[1, 'bitscore']) >= samplecfg.max_score_diff) & \
            (blastgroup.loc[0, 'evalue'] / blastgroup.loc[1, 'evalue'] <= samplecfg.max_score_diff)) & \
            (legacy_mode == 0):
        # add to uniquedf
        uniquedict.append(blastgroup.iloc[0,:].to_dict())
    elif ((abs(blastgroup.loc[0, 'bitscore'] - blastgroup.loc[1:,'bitscore'].max()) >= samplecfg.max_score_diff) & \
          (blastgroup.loc[0, 'evalue'] / blastgroup.loc[1:,'evalue'].max() <= samplecfg.max_score_diff)) & \
            (legacy_mode == 1):
        # add to uniquedf
        uniquedict.append(blastgroup.iloc[0,:].to_dict())
    else:
        # add to multidf the rows within the threshold ranges of top match, add all if legacy_mode
        if legacy_mode == 0:
            minbit = blastgroup.loc[0, 'bitscore'] - samplecfg.max_score_diff
            maxevalue = blastgroup.loc[0, 'evalue'] / samplecfg.max_score_diff
            blastgroup = blastgroup[(blastgroup['evalue'] <= maxevalue) & (blastgroup['bitscore'] >= minbit)]
        [multidict.append(blastgroup.iloc[i,:].to_dict()) for i in range(len(blastgroup))]

    return uniquedict, multidict


def add_rowinfo(blastrow, samplecfg):
    """Determine position, orientation, dupl number, indpt number and Id"""
    info = BlastRow(blastrow, samplecfg)
    blastrow['position'] = info.position
    blastrow['strand'] = info.strand
    blastrow['dupl'] = info.dupl
    blastrow['indpt'] = info.indpt
    blastrow['Id'] = info.Id
    blastrow['SNs'] = info.SNs if samplecfg.SN else 0
    return blastrow

def collapse(df):
    """Collapse ids mapped to identical strand-specific positions
    and count number of duplicates and independents,
    keep track of ids per Id"""
    agg_dict = {'chr':'first',
                'pident':'first',
                'length':'first',
                'mismatch':'first',
                'gapopen':'first',
                'qstart':'first',
                'qend':'first',
                'sstart':'first',
                'send':'first',
                'evalue':'first',
                'bitscore':'first',
                'position':'first',
                'strand':'first',
                'dupl':'sum',
                'indpt':'sum',
                'id':(lambda x: list(x)),
                'SNs':(lambda x: list(x))}
    cdf = df.groupby('Id').agg(agg_dict)
    def SNcount(x):
        x = x.replace('[','').replace(']','').replace("'",'').replace(' ','')
        SNs = x.split(sep=',')
        return len(set(SNs))
    cdf['indpt'] = cdf['SNs'].apply(lambda x: SNcount(str(x)))
    return cdf

class MakeLog:
    """Count number of insertions to generate the log file
    total_seq = number of sequence IDs
    lines = number of blast matches
    start_over1 = matches starting at bp >1
    high_evalue = matches with evalues > threshold
    unique_SSP = total number of uniquely mapped strand-specific positions
    unique_indpt = independent integration events uniquely mapped
    unique_dupl = number of sequence reads uniquely mapped
    multi_SSP = total number of multi mapped strand-specific positions
    multi_indpt = independent integration events multi mapped
    multi_dupl = number of sequence reads multi mapped
    NOTE: the number of indpt and dupl for a given sequence ID in multimapped are counted at each possible position
    """
    # format is not matching the log file from Perl script
    def __init__(self, samplecfg, rawblast, uniquedf, multidf):
        self.total_seq = len(rawblast['id'].unique())
        self.lines = len(rawblast.index)
        self.start_over1 = len(rawblast[rawblast['qstart'] > 1].index)
        self.high_evalue = len(rawblast[rawblast['evalue'] <= samplecfg.error].index)
        self.unique_ssp = self.per_chro(uniquedf, 'chr', 'dupl', len) if len(uniquedf.index) > 0 else 0
        self.unique_indpt = self.per_chro(uniquedf, 'chr', 'indpt', sum) if len(uniquedf.index) > 0 else 0
        self.unique_dupl = self.per_chro(uniquedf, 'chr', 'dupl', sum) if len(uniquedf.index) > 0 else 0
        self.multi_ssp = self.per_chro(multidf, 'chr', 'dupl', len) if len(multidf.index) > 0 else 0
        self.multi_indpt = self.per_chro(multidf, 'chr', 'indpt', sum) if len(multidf.index) > 0 else 0
        self.multi_dupl = self.per_chro(multidf, 'chr', 'dupl', sum) if len(multidf.index) > 0 else 0

    def per_chro(self, df, chro_col, nb_col, funct):
        # sum number from nb_col per chromosome in chro_col
        chrodict = {}
        chrodict['total'] = funct(df[nb_col])
        for chro in config[config['chro_listvs'][config['genome']]]:
            chrodict[chro] = funct(df[df[chro_col] == chro][nb_col])
        return chrodict


def main():
    py_cols = {'chr':'chr',
               'position':'position',
               'strand':'strand',
               'indpt':'independent_events',
               'dupl':'duplicate_reads'}
    legacy_cols = {'chr':'chr',
                   'pident':'identity',
                   'length':'length',
                   'mismatch':'mismatch',
                   'gapopen':'gap',
                   'qstart':'start',
                   'qend':'end',
                   'sstart':'chrStart',
                   'send':'chrEnd',
                   'evalue':'e',
                   'bitscore':'bitScore',
                   'chr':'chr',
                   'position':'topStrandCoordinate',
                   'strand':'strand',
                   'indpt':'targetSeq',
                   'dupl':'duplicate'}

    for sample in samples:
        logdf = pd.DataFrame(columns = samples,
                           index = ['total_seq', 'lines', 'start_over1', 'high_evalue',
                                    'unique_ssp', 'unique_indpt', 'unique_dupl',
                                    'multi_ssp', 'multi_indpt', 'multi_dupl',])

        samplecfg = SampleConfig(config, sample, args)
        samplename = samplecfg.sample + '.' + config['genomevs'][config['genome']]
        check_lib_design(samplecfg)
        col_names = ['id', 'chr', 'pident', 'length', 'mismatch', 'gapopen',
                               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        if os.path.exists(samplecfg.blastfn):
            rawblast = pd.read_csv(samplecfg.blastfn, sep='\t', header=None,
                        names=col_names)
            # discard any blast match starting at > 1 or with evalue > config['blastevalue']
            blast = rawblast[(rawblast['qstart'] == 1) & (rawblast['evalue'] <= samplecfg.error)]
            blast = blast.sort_values(['id', 'bitscore', 'evalue'], ascending=[True, False, True])
            # determine if uniquely mapped or multimapped
            uniquedict = []
            multidict = []
            dataframes = [group for _, group in blast.groupby('id')]
            for blastgroup in dataframes:
                uniquedict, multidict = is_uniquelymapped(blastgroup, samplecfg, uniquedict, multidict)
            # add info about number of dupl, ...
            for ud in uniquedict:
                add_rowinfo(ud, samplecfg)
            for md in multidict:
                add_rowinfo(md, samplecfg)
            uniquedf = pd.DataFrame.from_dict(uniquedict)
            multidf = pd.DataFrame.from_dict(multidict)
            # collapse per strand-specific position
            uniquedf = collapse(uniquedf) if len(uniquedf.index) > 0 \
                else pd.DataFrame(columns = legacy_cols.keys())
            multidf = collapse(multidf) if len(multidf.index) > 0 \
                else pd.DataFrame(columns = legacy_cols.keys())
            # generate log
            samplelog = MakeLog(samplecfg, rawblast, uniquedf, multidf)
            # print out mapping Id seq_id files
            uniquefn = mydir + '/id_mappings_' + samplename + '.txt'
            multifn = mydir + '/id_mappings_multimatch_' + samplename + '.txt'
            idcols = ['id', 'SNs'] if samplecfg.SN else ['id']
            if all(i in uniquedf.columns  for i in idcols):
                uniquedf[idcols].to_csv(uniquefn, sep='\t')
            if all(i in multidf.columns  for i in idcols):
                multidf[idcols].to_csv(multifn, sep='\t')
            # filter out chromosomes not in the full or short chro list config[config['chro_listvs'][config['genome']]
            chro_list = config[config['chro_listvs'][config['genome']]]
            uniquedf = uniquedf[uniquedf['chr'].isin(chro_list)]
            multidf = multidf[multidf['chr'].isin(chro_list)]
            # replace 0-indpt by na if not Serial Number
            if samplecfg.SN == 0:
                samplelog.unique_indpt = 'na'
                samplelog.multi_indpt = 'na'
                uniquedf['indpt'] = 'na'
                multidf['indpt'] = 'na'
            logdf[sample] = vars(samplelog).values()
            # sort by chr and position
            uniquedf = uniquedf.sort_values(by=['chr', 'position'])
            multidf = multidf.sort_values(by=['chr', 'position'])
            # reformat to match perl scripts
            output_cols = legacy_cols if legacy_mode else py_cols
            uniquedf = uniquedf[output_cols.keys()].rename(columns = output_cols)
            multidf = multidf[output_cols.keys()].rename(columns = output_cols)
            uniquefn = mydir + '/integration_' + samplename + '.txt'
            uniquedf.to_csv(uniquefn, sep='\t')
            multifn = mydir + '/integration_multimatch_' + samplename + '.txt'
            multidf.to_csv(multifn, sep='\t')
        else:
            logdf[sample] = 'file does not exist: ' + samplecfg.blastfn

        # log file
        logfn = logdir + '/integration_' + config['name'] + '_' + samplename + '.log.txt'
        logdf.to_csv(logfn, sep='\t')


main()


