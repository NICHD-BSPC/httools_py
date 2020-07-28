#!/usr/bin/env python


import os
from pathlib import Path
import ruamel.yaml
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import pairwise2
import csv
import gzip
import binascii
import re
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator


ap = argparse.ArgumentParser()
ap.add_argument('--config', help='config.yaml file', default='config/config.yaml')
ap.add_argument('--sample', help='sample', default=False)
args = ap.parse_args()

yaml = ruamel.yaml.YAML()
config = yaml.load(open(args.config))
samples = [args.sample] if args.sample else config['sample'].keys()
legacy_mode = config['legacy_mode']

logdir = os.path.dirname(args.config) + '/data/' + config['name'] + '/logs'
Path(logdir).mkdir(parents=True, exist_ok=True)
mydir = os.path.dirname(args.config)  +'/data/' + config['name'] + '/fastqscreen'
Path(mydir).mkdir(parents=True, exist_ok=True)


class SampleConfig:
    """Make sample-specific config"""

    def __init__(self, config, sample):
        self.sample = sample
        self.seq = config['sample'][sample]['sequence']
        self.barcode_start = config['sample'][sample]['barcode_start']
        self.barcode = self.seq[config['sample'][sample]['barcode_start']-1:config['sample'][sample]['barcode_length']]
        self.integrase = config['sample'][sample]['integrase']
        self.lib_design = config['sample'][sample]['lib_design']
        self.length_to_match = config['length_to_match']
        self.min_length = config['min_length']
        self.allowed_mismatches = config['allowed_mismatches']
        self.linker = config['linker']
        self.ltrcircle = config['ltrcircle'][self.lib_design]
        self.plasmid = config['plasmid'][self.lib_design]
        self.primary_re = config['primary_re'][self.lib_design]
        self.primary_incomplete = config['primary_incomplete'][self.lib_design]
        self.second_re = config['second_re'][self.lib_design]
        self.second_incomplete = config['second_incomplete'][self.lib_design]
        self.dist_to_second_incomplete = config['dist_to_second_incomplete'][self.lib_design]
        self.pbs = config['pbs'][self.lib_design]
        self.tsd = config['tsd'][self.integrase]
        if isinstance(config['sample'][sample]['SN_position'], int) and isinstance(config['sample'][sample]['SN_length'], int):
            self.sn_position = config['sample'][sample]['SN_position']
            self.sn_length = config['sample'][sample]['SN_length']
        else:
            self.sn_position = 0
            self.sn_length = 0
        self.end_of_ltr = self.seq[-(self.length_to_match + self.sn_length) :]


class ScreenFastq:
    """Parsing and filtering of fastq record"""

    def __init__(self, seq, seqid, samplecfg):
        self.seq = seq
        self.id = seqid
        self.samplecfg = samplecfg
        self.barcode = self.has_barcode(self.seq, samplecfg)
        if self.barcode:
            self.end_ltr = self.has_end_ltr(self.seq, samplecfg)
            if self.end_ltr:
                self.plasmid = self.has_plasmid(self.seq, samplecfg)
                self.ltrcircle = self.has_ltrcircle(self.seq, samplecfg)
                if (samplecfg.sn_length):
                    self.second_incomplete = self.has_second_incomplete(self.seq, samplecfg) \
                        if ((self.plasmid + self.ltrcircle) == 0) else 0
                else:
                    self.second_incomplete = self.has_second_incomplete(self.seq, samplecfg)
                self.primary_incomplete = self.has_primary_incomplete(self.seq, samplecfg) \
                    if ((self.plasmid + self.ltrcircle + self.second_incomplete) == 0) else 0
                self.pbs = self.has_pbs(self.seq, samplecfg)
                self.tmplinker = self.has_linker(self.seq, samplecfg)
                self.trimmed = self.trim_seq(self.seq, samplecfg)
                self.short = self.is_short(self.trimmed, samplecfg) if sum( \
                                [self.plasmid, self.primary_incomplete,self.second_incomplete, self.ltrcircle]) == 0 \
                                else 0
                self.filters = sum([self.plasmid, self.primary_incomplete,self.second_incomplete, self.ltrcircle, self.short])
                self.linker = self.has_linker(self.seq, samplecfg) if ((self.filters - self.short) == 0) else 0
                self.pass_good = 1 if (self.filters == 0 and self.pbs == 0) else 0
                self.pass_pbs = 1 if (self.filters == 0 and self.pbs == 1) else 0
                self.sn = self.get_sn(self.seq, samplecfg)
    # NOTE: barcode and end_ltr will have the most screened out sequences, so continue the screen
    # only if they pass
    # linker is only counted if this is a good sequence, whether or not is it short

    def skipped(self, criteria):
        """ Determine whether the sequence should be filtered on the samplecfg criteria
        """
        return 1 if criteria == 'na' else 0

    def exactmatch(self, seq, match, start):
        """Look for exact match of 'match' within seq, starting at index 'start'-1 of seq
        """
        if len(seq) < start + len(match):
            score = 0
        else:
            score = 1 if seq[(start-1):start-1+len(match)] == match else 0
        return score

    def fuzzymatch(self, seq, match, start, maxmis):
        """Calculate match score of 'match' within the subset of seq starting at index 'start'-1 and
        of the same length as 'match'
        Compare with max score minus max number of allowed mismatches
        maxmis: max mismatches allowed
        """
        self.subseq = seq[start-1:(start-1+len(match))]
        if len(self.subseq) < len(match):
            score = 0
        else:
            score = 1 if pairwise2.align.localxd(self.subseq, match, -1,-1,0,0, score_only=True) \
               >= len(match) - maxmis else 0
        return score


    def SN_fuzzymatch(self, seq, match, start, maxmis, sn_position, sn_length):
        """Same as fuzzymatch, but with >
        """
        self.subseq1 = seq[start-1:sn_position-1]
        self.subseq2 = seq[sn_position-1+sn_length:(start-1+len(match))]
        self.match1 = match[0:sn_position-start]
        self.match2 = match[sn_position-start+sn_length:]
        if len(self.subseq1)+len(self.subseq2)+sn_length < len(match):
            score = 0
        else:
            self.score1 = pairwise2.align.localxd(self.subseq1, self.match1, -1,-1,0,0, score_only=True)
            self.score2 = pairwise2.align.localxd(self.subseq2, self.match2, -1,-1,0,0, score_only=True)
            score = 1 if self.score1 + self.score2 >= len(match) - maxmis else 0
        return score

    def slidingmatch(self, seq, match, maxmis):
        """Calculate match score of 'match' anywhere within seq
        Compare with max score minus number of allowed mismatches
        maxmis: max mismatches allowed
        Return position of match (index + 1)
        """
        # pairwise2.align does bad when aligning 2 sequences of different size without gap penalty. 
        # But don't want gap penalty in the match score, so do alignment first with gap penalty to get
        # the region of the match, then do second alignment without gap penalty of the subseq to get
        # the score
        # only used for LTR circle, which are not many. 1st pass screen with scores only to speed things up
        if pairwise2.align.localxs(seq, match,-1,-1, score_only=True) >= len(match) - maxmis:
            self.s, self.m, self.score, self.start, self.end = pairwise2.align.localxs(seq, match,-1,-1)[0]
            self.subseq = seq[self.start:self.end]
            if pairwise2.align.localxd(self.subseq, match, -1,-1,0,0, score_only=True):
                return self.start+1 if pairwise2.align.localxd(self.subseq, match, -1,-1,0,0, score_only=True) \
                    >= self.end - self.start - maxmis else 0
            else:
                return 0
        else:
            return 0

    def has_barcode(self, seq, samplecfg):
        return self.exactmatch(seq, samplecfg.barcode, samplecfg.barcode_start)

    def has_end_ltr(self, seq, samplecfg):
        # perl script screen_illumina_Tf1_sequence-1.0.pl has if >= allowedmismatches is screened out.
        # should have been > because as it is it only allows for 1 mismatch
        maxmis = samplecfg.sn_length + samplecfg.allowed_mismatches -1 if legacy_mode \
            else samplecfg.sn_length + samplecfg.allowed_mismatches
        score = self.SN_fuzzymatch(seq,
                          samplecfg.end_of_ltr,
                          len(samplecfg.seq)-len(samplecfg.end_of_ltr) + 1,
                          maxmis,
                          samplecfg.sn_position,
                          samplecfg.sn_length) if samplecfg.sn_length else \
                self.fuzzymatch(seq,
                          samplecfg.end_of_ltr,
                          len(samplecfg.seq)-len(samplecfg.end_of_ltr) + 1,
                          maxmis)
        return score

    def has_plasmid(self, seq, samplecfg):
        if self.skipped(samplecfg.plasmid):
            return 0
        else:
            return self.fuzzymatch(seq,
                          samplecfg.plasmid,
                          len(samplecfg.seq) + 1,
                          samplecfg.allowed_mismatches)

    def has_primary_incomplete(self, seq, samplecfg):
        if self.skipped(samplecfg.primary_incomplete):
            return 0
        else:
            # match immediately after end of LTR or 1bp after
            return self.exactmatch(seq,
                          samplecfg.primary_incomplete,
                          len(samplecfg.seq) + 1) + \
               self.exactmatch(seq,
                          samplecfg.primary_incomplete,
                          len(samplecfg.seq) + 2)

    def has_second_incomplete(self, seq, samplecfg):
        if self.skipped(samplecfg.second_incomplete):
            return 0
        else:
            # perl version screens out when distance <= $allowedmismatch, this allows only for 1 mismatch NO, ALLOWS FOR 2
            # distance should be calculated on the available sequence read length
            maxmis = samplecfg.allowed_mismatches #-1 if legacy_mode \
                #else samplecfg.allowed_mismatches
            trimmed_second_incomplete = samplecfg.second_incomplete[
                0:(len(seq) - len(samplecfg.seq) - samplecfg.dist_to_second_incomplete+1)
                ]
            return self.fuzzymatch(seq,
                          trimmed_second_incomplete,
                          len(samplecfg.seq) + samplecfg.dist_to_second_incomplete,
                          maxmis)


    def has_pbs(self, seq, samplecfg):
        if self.skipped(samplecfg.pbs):
            return 0
        else:
            return self.fuzzymatch(seq,
                          samplecfg.pbs,
                          len(samplecfg.seq) + 1,
                          1)

    def has_ltrcircle(self, seq, samplecfg):
        if self.skipped(samplecfg.ltrcircle):
            return 0
        else:
            return 1 if self.slidingmatch(seq,
                            samplecfg.ltrcircle,
                            samplecfg.allowed_mismatches) > 0 else 0

    def has_linker(self, seq, samplecfg):
        if self.skipped(samplecfg.linker):
            return 0
        else:
            # linker should be exact match anywhere, look for exact seq in read
            return (samplecfg.linker in seq)

    def is_short(self, seq, samplecfg):
        return 1 if len(seq) <= samplecfg.min_length else 0

    def trim_seq(self, seq, samplecfg):
        trimmedseq = seq[len(samplecfg.seq):]
        if self.skipped(samplecfg.linker):
            return(trimmedseq)
        else:
            trimmedseq = trimmedseq[:trimmedseq.find(samplecfg.linker)+1] if samplecfg.linker in trimmedseq \
                else trimmedseq
            return(trimmedseq)

    def get_sn(self, seq, samplecfg):
        return str(seq[samplecfg.sn_position -1 : samplecfg.sn_position -1 + samplecfg.sn_length])




def _is_gzipped(fn):
    """
    Filename-independent method of checking if a file is gzipped or not. Uses
    the magic number.
    xref https://stackoverflow.com/a/47080739
    """
    with open(fn, 'rb') as f:
        return binascii.hexlify(f.read(2)) == b'1f8b'

def openfile(tmp, mode):
    """
    Returns an open file handle; auto-detects gzipped files.
    """
    if _is_gzipped(tmp):
        return gzip.open(tmp, mode)
    else:
        return open(tmp, mode)


def collapse_seqs(df, fn, idx='seq', ID='id', SN=False):
    """ Collapse dataframe by sequences, count and list the group ids
    Write to fasta file
    Return number of unique sequences
    """
    gpdf = df.groupby([idx]).count()
    gpdf['ids'] = pd.Series(df.groupby([idx])[ID].aggregate(lambda ids: ids.tolist()))
    if SN:
        gpdf['SN'] = pd.Series(df.groupby([idx])[SN].aggregate(lambda SNs: SNs.unique().tolist()))
        gpdf['nb_SN'] = pd.Series(df.groupby([idx])[SN].aggregate(lambda SNs: len(SNs.unique().tolist())))
        with open(fn, 'w') as fnfile:
            fnfile.writelines('>seq' + str(gpdf.index.get_loc(seq) +1) +
                              '_dupl' + str(gpdf.loc[seq, 'id']) +
                              '_indpt' + str(gpdf.loc[seq, 'nb_SN']) +
                              '\n' + seq +
                              '\n' for seq in gpdf.index)
        with open(os.path.splitext(fn)[0] + '_SN.txt', 'w') as snfile:
            snfile.writelines('>seq' + str(gpdf.index.get_loc(seq) +1) +
                              '_dupl' + str(gpdf.loc[seq, 'id']) +
                              '_indpt' + str(gpdf.loc[seq, 'nb_SN']) +
                              '\n' + str(gpdf.loc[seq, 'SN']) +
                              '\n' for seq in gpdf.index)
        with open(os.path.splitext(fn)[0] + '_ids.txt', 'w') as idfile:
            idfile.writelines('>seq' + str(gpdf.index.get_loc(seq) +1) +
                              '_dupl' + str(gpdf.loc[seq, 'id']) +
                              '_indpt' + str(gpdf.loc[seq, 'nb_SN']) +
                              '\n' + str(gpdf.loc[seq, 'ids']) +
                              '\n' for seq in gpdf.index)
        return len(gpdf), str(gpdf['nb_SN'].sum())
    else:
        with open(fn, 'w') as fnfile:
            fnfile.writelines('>seq' + str(gpdf.index.get_loc(seq) +1) +
                              '_dupl' + str(gpdf.loc[seq, 'id']) +
                              '\n' + seq +
                              '\n' for seq in gpdf.index)
        with open(os.path.splitext(fn)[0] + '_ids.txt', 'w') as idfile:
            idfile.writelines('>seq' + str(gpdf.index.get_loc(seq) +1) +
                              '_dupl' + str(gpdf.loc[seq, 'id']) +
                              '\n' + str(gpdf.loc[seq, 'ids']) +
                              '\n' for seq in gpdf.index)
        return [len(gpdf), 'na']

def screen_sample(config, samplecfg, counter, mydir, logdir):
    """ Filter fastq for sequences corresponding to sample
    Output fastqs and logs
    Return number of seq in fastq
    """
    nseq = 0
    seqdf = []
    pbsdf = []
    cols = counter.columns
    counterd = dict.fromkeys(samples)
    for sample in samples:
        counterd[sample] = dict.fromkeys(cols, 0)

    for fqfn in config['fastq']:
        print('\nProcessing ' + fqfn + '\n')
        fqpath = os.path.join(os.path.dirname(args.config), fqfn)
        with openfile(fqpath, 'rt') as fqfile:
            for seqid, seq, qual in FastqGeneralIterator(fqfile):
                nseq += 1
                # screen seq
                currentseq = ScreenFastq(seq, seqid, samplecfg)
                #print(nseq)
                # add to counter
                for k in cols:
                    counterd[samplecfg.sample][k] += getattr(currentseq, k, 0)
                # add seq to dict if passed all filters
                if hasattr(currentseq, 'pass_good') and currentseq.pass_good:
                    seqdf.append({'id':currentseq.id, 'seq':str(currentseq.trimmed), 'sn':currentseq.sn})
                if hasattr(currentseq, 'pass_pbs') and currentseq.pass_pbs:
                    pbsdf.append({'id':currentseq.id, 'seq':str(currentseq.trimmed), 'sn':currentseq.sn})

    seqdf = pd.DataFrame(seqdf)
    pbsdf = pd.DataFrame(pbsdf)

    # convert counterd in to counter (the dataframe)
    counter = pd.DataFrame.from_dict(counterd, orient='index')

    # group df by sequence and write to fasta, number of unique seq returned
    SN = 'sn' if samplecfg.sn_length > 0 else False
    counter.loc[samplecfg.sample, 'unique_good'], counter.loc[samplecfg.sample, 'indpt_good'] = collapse_seqs(
        seqdf, mydir + '/screen_' + str(samplecfg.sample) + '.fa', SN = SN) if len(seqdf) > 0 else [0, 0]
    counter.loc[samplecfg.sample, 'unique_pbs'], counter.loc[samplecfg.sample, 'indpt_pbs'] = collapse_seqs(
        pbsdf, mydir + '/screen_PBS_' + str(samplecfg.sample) + '.fa', SN = SN) if len(pbsdf) > 0 else [0, 0]

    return [counter, nseq, len(seqdf)]

def check_sn_config(samplecfg):
    snseq = samplecfg.seq[samplecfg.sn_position -1 : samplecfg.sn_position -1 + samplecfg.sn_length]
    if re.search(r'[aAtTcCgGnN]', snseq):
        raise ValueError('Invalid format for Serial Number ' + snseq + ' in ' + samplecfg.sample,
                         '\nCannot contain aAtTcCgGnN\n')

def no_read_error(samplecfg, counter, nseq, nobc):
    """ Return an error if no sequence passed all the filters
    Write out an error file containing the number of sequences that passed each filter,
    for easier debugging
    """
    errmessage = 'Zero read passing the sequence-specific filters for ' + samplecfg.sample + \
                 '. A common reason for this is an error in sequences specified in the config yaml file\n'
    counterfn = logdir + '/fastq_screen_' + config['name'] + '_' + samplecfg.sample + '.error.txt'
    with open(counterfn, 'w+') as logfile:
        logfile.write('# Error: ' + errmessage)
        logfile.write('Raw sequences: ' + str(nseq) + '\n' + \
                      'Reads without barcode: ' + str(nobc) + '\n\n')
    counter.to_csv(counterfn, sep='\t', mode='a')
    raise ValueError(errmessage)


def main():
    # initialize counter dataframe
    screen_criteria = ['barcode', 'end_ltr', 'plasmid', 'primary_incomplete', 'second_incomplete',
                       'ltrcircle', 'linker', 'short', 'pass_good', 'unique_good', 'pass_pbs', 'unique_pbs']
    skiplist = ['plasmid', 'primary_incomplete', 'second_incomplete', 'pbs', 'ltrcircle', 'linker']
    py_cols = {'Barcode start':'Barcode start',
                   'Barcode length':'Barcode length',
                   'SN start':'SN start',
                   'SN length':'SN length',
                   'Transposon end':'Transposon end',
                   'barcode':'With barcode',
                   'Nonspecific':'Nonspecific',
                   'plasmid':'plasmid',
                   'second_incomplete': 'NA',
                   'ltrcircle':'LTRcircle',
                   'primary_incomplete': 'NA',
                   'With LTR':'With LTR',
                   'linker':'With linker',
                   'short':'Short',
                   'pass_good':'Good sequences',
                   'unique_good':'Exclude dup',
                   'indpt_good':'Exclude SN',
                   'dup':'dup',
                   'pass_pbs':'PBS events Good sequences',
                   'indpt_pbs':'PBS events Exclude dup',
                   'unique_pbs':'PBS events Exclude SN',
                   'PBS events dups':'PBS events dups'
                   }
    legacy_SN_cols = {'Barcode start':'Barcode start',
                   'Barcode length':'Barcode length',
                   'SN start':'SN start',
                   'SN length':'SN length',
                   'Transposon end':'Beginning seq',
                   'barcode':'With barcode',
                   'Nonspecific':'Nonspecific',
                   'plasmid':'plasmid',
                   'ltrcircle':'LTRcircle',
                   'second_incomplete': 'NA',
                   'primary_incomplete': 'NA',
                   'With LTR':'With LTR',
                   'linker':'With linker',
                   'short':'Short',
                   'pass_good':'Good sequences',
                   'unique_good':'Exclude SN',
                   'pass_pbs':'PBS events Good sequences',
                   'unique_pbs':'PBS events Exclude SN',
                   }
    legacy_nonSN_cols = {'Barcode start':'Barcode start',
                   'Barcode length':'Barcode length',
                   'Transposon end':'Transposon end',
                   'barcode':'With barcode',
                   'Nonspecific':'Nonspecific',
                   'plasmid':'plasmid',
                   'second_incomplete': 'NA',
                   'ltrcircle':'LTRcircle',
                   'primary_incomplete': 'NA',
                   'With LTR':'With LTR',
                   'linker':'With linker',
                   'short':'Short',
                   'pass_good':'Good sequences',
                   'unique_good':'Exclude dup',
                   'dup':'dup',
                   'pass_pbs':'PBS events Good sequences',
                   'unique_pbs':'PBS events Exclude dup',
                   'PBS events dups':'PBS events dups'
                   }


    # get sample-specific config
    for sample in samples:
        counter = pd.DataFrame(index=samples,
                           columns=screen_criteria).fillna(0)
        samplecfg = SampleConfig(config, sample)
        check_sn_config(samplecfg)
        counter, nseq, ngood = screen_sample(config, samplecfg, counter, mydir, logdir)
        # columns to match the old version
        counter.loc[sample, 'Barcode start'] = str(config['sample'][sample]['barcode_start'])
        counter.loc[sample, 'Barcode length'] = str(config['sample'][sample]['barcode_length'])
        counter.loc[sample, 'SN start'] = str(config['sample'][sample]['SN_position']) if samplecfg.sn_length > 0 else 'na'
        counter.loc[sample, 'SN length'] = str(config['sample'][sample]['SN_length']) if samplecfg.sn_length > 0 else 'na'
        counter.loc[sample, 'Transposon end'] = config['sample'][sample]['sequence']
        counter.loc[sample, 'Nonspecific'] = str(counter.loc[sample, 'barcode'] - counter.loc[sample, 'end_ltr'])
        counter.loc[sample, 'With LTR'] = str(counter.loc[sample, 'pass_good'] + counter.loc[sample, 'short'] + counter.loc[sample, 'pass_pbs'])
        counter.loc[sample, 'dup'] = str(counter.loc[sample, 'pass_good'] - counter.loc[sample, 'unique_good'])
        counter.loc[sample, 'PBS events dups'] = str(counter.loc[sample, 'pass_pbs'] - counter.loc[sample, 'unique_pbs'])
        if legacy_mode:
            output_cols = legacy_SN_cols if samplecfg.sn_length > 0 else legacy_nonSN_cols
            for cols in ['Nonspecific', 'plasmid', 'second_incomplete', 'primary_incomplete', 'ltrcircle', 'With LTR', 
                     'linker', 'short', 'pass_good', 'pass_pbs']:
                counter[[cols]] = np.where((counter[[cols]] == 0) | (counter[[cols]] == '0'), '', counter[[cols]])
        else:
            output_cols = py_cols
        output_cols['primary_incomplete'] = samplecfg.primary_re +' incomplete'
        output_cols['second_incomplete'] = samplecfg.second_re +' incomplete'

        [output_cols.pop(k, None) for k in skiplist if getattr(samplecfg,k) == 'na'] 

        counter = counter[output_cols.keys()].rename(columns = output_cols)
        counter = counter.rename_axis('Name')

        # add total seq and no barcode count and write counter to log
        nobc = nseq - counter['With barcode'].sum()
        if legacy_mode:
            counter[['With barcode']] = np.where((counter[['With barcode']] == 0) | (counter[['With barcode']] == '0'), '', counter[['With barcode']])
        if nobc == 0 and legacy_mode:
            nobc = ''
        counterfn = logdir + '/fastq_screen_' + config['name'] + '_' + sample + '.log.txt'

        # return an error if there are no sequences in seqdf
        if ngood == 0:
            no_read_error(samplecfg, counter, nseq, nobc)

        with open(counterfn, 'w+') as logfile:
            if legacy_mode and samplecfg.sn_length > 0:
                logfile.write('Raw sequences: ' + str(nseq) + '\n' + \
                      'Reads without barcode: ' + str(nobc) + ', reads\n\n')
            else:
                logfile.write('Raw sequences: ' + str(nseq) + '\n' + \
                      'Reads without barcode: ' + str(nobc) + '\n\n')
        counter.to_csv(counterfn, sep='\t', mode='a')

main()
