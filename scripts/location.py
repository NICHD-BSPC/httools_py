#!/usr/bin/env python


import os
from pathlib import Path
import ruamel.yaml
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import pybedtools

ap = argparse.ArgumentParser()
ap.add_argument('--config', help='config.yaml file', default='../config/config.yaml')
ap.add_argument('--integration', help='integration file', default=False)
ap.add_argument('--cds', help='CDS coordinate file', default=False)
ap.add_argument('--sample', help='sample', default=False)
args = ap.parse_args()

yaml = ruamel.yaml.YAML()
config = yaml.load(open(args.config))
samples = [args.sample] if args.sample else config['sample'].keys()
legacy_mode = config['legacy_mode']

logdir = 'data/' + config['name'] + '/logs'
Path(logdir).mkdir(parents=True, exist_ok=True)
mydir = 'data/' + config['name'] + '/location'
Path(mydir).mkdir(parents=True, exist_ok=True)
hmdir = 'data/' + config['name'] + '/location/homologous_recombination'
Path(hmdir).mkdir(parents=True, exist_ok=True)
exdir = 'data/' + config['name'] + '/location/excluded'
Path(exdir).mkdir(parents=True, exist_ok=True)

class SampleConfig:
    """Make sample-specific config"""
    def __init__(self, config, sample, args):
        self.sample = sample
        self.lib_design = config['sample'][sample]['lib_design']
        self.tsd = config['tsd'][config['sample'][sample]['integrase']]
        self.SN = 1 if isinstance(config['sample'][sample]['SN_position'], int) and isinstance(config['sample'][sample]['SN_length'], int) else 0
        self.integfn= self.get_integfn(sample, config, args)
        self.integdf = self.get_integ(self.integfn)[0] if os.path.exists(self.integfn) else 0
        self.integbed = self.get_integ(self.integfn)[1] if os.path.exists(self.integfn) else 0
        self.chromsizes = self.get_chromsizes(config)
        self.cdsfn = args.cds if args.cds else self.get_path(config['genomecds'][config['genome']])
        self.cds = self.get_bed(self.cdsfn)
        self.ltr5 = 0 if self.skipped(config['preexist_ltr'][self.lib_design]['ltr5']) \
            else self.get_bed(self.get_path(config['preexist_ltr'][self.lib_design]['ltr5']))
        self.ltr3 = 0 if self.skipped(config['preexist_ltr'][self.lib_design]['ltr3']) \
            else self.get_bed(self.get_path(config['preexist_ltr'][self.lib_design]['ltr3']))
        self.sololtr = 0 if self.skipped(config['preexist_ltr'][self.lib_design]['sololtr']) \
            else self.get_bed(self.get_path(config['preexist_ltr'][self.lib_design]['sololtr']))
        self.exclude = config['exclude']

    def skipped(self, criteria):
        """ Determine whether the filtering step should be skipped
        """
        return 1 if criteria == 'none' else 0

    def get_integfn(self, sample, config, args):
        integfn = args.integration if args.integration else 'data/' + config['name'] + \
            '/filblast/integration_' + sample + '.' + config['genomevs'][config['genome']] + '.txt'
        return integfn

    def get_integ(self, integfn):
        """ Return dataframe and pybedtools object of integration file"""
        df = pd.read_csv(integfn, sep='\t', index_col=0)
        df = df.rename(columns={'topStrandCoordinate':'position', 'duplicate':'dupl', 'targetSeq':'indpt'}) if legacy_mode else \
            df.rename(columns={'duplicate_reads':'dupl', 'independent_events':'indpt'})
        col_names = ['chr', 'position', 'position', 'id', 'dupl', 'strand']
        df['id'] = df.index
        bed = pybedtools.BedTool.from_dataframe(df[col_names])
        return df, bed

    def get_chromsizes(self, config):
        chromsizes = {}
        for rec in SeqIO.parse(self.get_path(config['genomedb'][config['genome']]), 'fasta'):
            chromsizes[rec.id] = len(rec.seq)
        return chromsizes

    def get_path(self, fn):
        return os.path.join(os.path.dirname(__file__), '..', fn)

    def get_bed(self, fn):
        # need to add mock CDS for start and end of chromosomes, for the first and last intergenic regions
        bed = pybedtools.BedTool(fn)
        extbed = pd.DataFrame(columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
        for chro in self.chromsizes.keys():
            extbed.loc[chro + '_start'] = [chro, 1, 1, chro+'_start', 0, 'na']
            extbed.loc[chro + '_end'] = [chro, self.chromsizes[chro]+1, self.chromsizes[chro]+1, chro+'_end', 0, 'na']
        extbed = pybedtools.BedTool.from_dataframe(extbed)
        return extbed.cat(bed, postmerge=False).sort()


def get_intergtype(uporfstrd, dnorfstrd):
    ori = [uporfstrd, dnorfstrd]
    if (ori == ['+', '+'] or ori == ['-', '-']):
        return 'Tandem'
    elif (ori == ['+', '-']):
        return 'Convergent'
    elif (ori == ['-', '+']):
        return 'Divergent'
    elif (ori[0] == 'na' or ori[1] == 'na'):
        return ''


class Location:
    """Position integration relative to adjacent ORFs/CDSs
    id, chr, position, strand, dupl, indpt refer to the integration
    uporf refers to the upstream ORF (relative to genome)
    dnorf refers to the downstream ORF (relative to genome)
    inorf refers to integration within ORF
    loc refers to the location of integration relative to nearest ORF (UP, DOWN, IN_start, IN_end)
    intergenic refers to the type of orientation of ORFs in the intergenic region (Divergent
        convergent, tandem, or -- if in ORF)
    distance refers to the distance of integration to the nearest ORF
    """
    def __init__(self, integ, samplecfg, tsdshiftCDS, cdsdf, legacy_mode):
        self.ID = integ['id']
        self.chro = integ['chr']
        self.position = integ['position']
        self.strand = integ['strand']
        self.dupl = integ['dupl']
        self.indpt = integ['indpt']
        self.integbed = self.integ_to_bed(self.ID, self.chro, self.position, self.strand)
        # location determined relative to the shifted CDS, then the coordinates are shifted back to original
        self.uporf = self.map_interg(self.integbed, tsdshiftCDS, cdsdf, fu=True)
        self.dnorf = self.map_interg(self.integbed, tsdshiftCDS, cdsdf, fd=True)
        self.closestorf = self.map_interg(self.integbed, tsdshiftCDS, cdsdf)
        # self.closestorf has original coordinates. Need to shift to determine location then shift back
        self.inorf_start = self.map_inorf(
            self.integbed, self.closestorf, cdsdf, samplecfg.tsd, fu=True) if self.closestorf['dist'][0] == 0 else 0
        self.inorf_end = self.map_inorf(
            self.integbed, self.closestorf, cdsdf, samplecfg.tsd, fd=True) if self.closestorf['dist'][0] == 0 else 0
        self.loc = self.get_loc(self.uporf, self.dnorf, self.closestorf, self.inorf_start, self.inorf_end)
        self.intergenic = '--' if self.closestorf['dist'][0] == 0 \
                            else get_intergtype(self.uporf['orfstrand'][0], self.dnorf['orfstrand'][0])
        self.distance = min(abs(self.inorf_start['dist'][0]), abs(self.inorf_end['dist'][0])) \
                            if self.closestorf['dist'][0] == 0 \
                            else min(abs(self.uporf['dist'][0]), abs(self.dnorf['dist'][0]))

    def integ_to_bed(self, ID, chro, position, strand):
        return pybedtools.BedTool.from_dataframe(
            pd.DataFrame([{'chrom':chro, 'chromStart':position, 'chromEnd':position,
                           'name':ID, 'score':0, 'strand':strand}]))

    def map_interg(self, integbed, tsdshiftedCDS, cdsdf, fu=False, fd=False):
        dfcols = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'orfchr', 'orfstart', 'orfend', 'orfname', 'orfscore', 'orfstrand',
                     'dist']
        # if ties in closest, take only the first occurance
        df = integbed.closest(tsdshiftedCDS, D = "ref",
                              fu=fu, fd=fd,
                              k=1, nonamecheck=True,
                              t='first').to_dataframe(header=None, names=dfcols)
        # revert TSD shift for the end coordinate of ORF now that the distance has been calculated
        df = self.tsd_reverse(df, cdsdf)
        return df[['orfname', 'orfstart', 'orfend', 'orfstrand','dist', 'chrom']]

    def map_inorf(self, integbed, orf, cdsdf, tsd, fu=False, fd=False):
        # create 2 intervals of 1bp per ORF for start and end, shift the end interval coordinate
        # NOTE: added -1 and +1 to start and end because the distances were shifted by 1bp compared to
        # the original Perl HTtools
        df = pd.DataFrame({'chrom':[orf['chrom'][0], orf['chrom'][0]],
                           'start': [orf['orfstart'][0]-1, orf['orfend'][0]+1-tsd-1] if legacy_mode else \
                                [orf['orfstart'][0]-1, orf['orfend'][0]+1-tsd+1],
                           'end' : [orf['orfstart'][0]-1, orf['orfend'][0]+1-tsd-1] if legacy_mode else \
                                [orf['orfstart'][0]-1, orf['orfend'][0]+1-tsd+1],
                           'name': [orf['orfname'][0], orf['orfname'][0]],
                           'score': [0, 0],
                           'strand':[orf['orfstrand'][0], orf['orfstrand'][0]]})
        orfbed = pybedtools.BedTool.from_dataframe(df)
        return self.map_interg(integbed, orfbed, cdsdf, fu=fu, fd=fd)

    def get_loc(self, uporf, dnorf, closestorf, inorf_start, inorf_end):
        if (closestorf['dist'][0] == 0 and abs(inorf_start['dist'][0]) <= abs(inorf_end['dist'][0])):
            direction = 'IN_start' if closestorf['orfstrand'][0] == '+' else 'IN_end'
        elif (closestorf['dist'][0] == 0 and abs(inorf_start['dist'][0]) > abs(inorf_end['dist'][0])):
            direction = 'IN_end' if closestorf['orfstrand'][0] == '+' else 'IN_start'
        elif (abs(uporf['dist'][0]) <= abs(dnorf['dist'][0])):
            direction = 'UP' if uporf['orfstrand'][0] == '-' else 'DOWN'
        elif (abs(uporf['dist'][0]) > abs(dnorf['dist'][0])):
            direction = 'DOWN' if dnorf['orfstrand'][0] == '-' else 'UP'
        # special case if start or end of chromosome, to match perl version
        if (uporf['orfstrand'][0] == 'na'):
            direction = 'DOWN' if dnorf['orfstrand'][0] == '-' else 'UP'
        if (dnorf['orfstrand'][0] == 'na'):
            direction = 'UP' if uporf['orfstrand'][0] == '-' else 'DOWN'
        if (dnorf['orfstrand'][0] == 'na' and uporf['orfstrand'][0] == 'na'):
            direction = 'IN_end'
        return(direction)

    def tsd_reverse(self, shifteddf, cdsdf):
        # reverse the shift of end coordinate of CDS (relative to genome) by TSD + 1
        shifteddf['orfend'] = shifteddf.apply(lambda x: cdsdf[cdsdf['name'] == x['orfname']]['end'], axis=1)
        shifteddf['orfstart'] = shifteddf.apply(lambda x: cdsdf[cdsdf['name'] == x['orfname']]['start'], axis=1)
        return(shifteddf)

def tsd_shift(cdsdf, tsd):
    # shift the end coordinate of CDS (relative to genome) by TSD + 1 to account for TSD
    # in assigning ORF vs. intergenic region to integration position
    df = cdsdf.copy()
    df['end'] = df['end'] - tsd - 1 if legacy_mode else \
        df['end'] - tsd +1
    # shift only to start coordinate if the interval is < tsd
    df['end'] = np.where(df['end'] < df['start'], df['start'], df['end'])
    return(pybedtools.BedTool.from_dataframe(df))


def exclude_to_bed(exclude):
    """ Takes list of positions to exclude and return BedTool object"""
    df = pd.DataFrame(exclude, columns=['name'])
    df[['chr', 'start', 'strand']] = df.name.str.split('_', expand=True)
    df['score'] = 1
    df = df[['chr', 'start', 'start', 'name', 'score', 'strand']]
    return(pybedtools.BedTool.from_dataframe(df))


def filter_homrecomb(samplecfg, config, fn):
    """ Parse integrations into homologous recombination, excluded and true integrations"""
    drop_cols = ['id'] if samplecfg.SN else ['id', 'indpt']
    # get positions to exclude
    exclude = exclude_to_bed(samplecfg.exclude) if samplecfg.exclude != ['none'] else False
    # true integrations
    beds = [x for x in [samplecfg.ltr5, samplecfg.ltr3, samplecfg.sololtr, exclude] if x ]
    if len(beds) > 0:
        trueint = samplecfg.integdf.loc[
            samplecfg.integbed.intersect(
                beds, v=True, F=1, nonamecheck=True
            ).to_dataframe()['name'], :]
    else:
        trueint = samplecfg.integdf
    trueint.drop(columns=drop_cols).to_csv(mydir + '/true_integration_' + fn + '.txt', sep='\t')
    # in LTR 5'
    idx = samplecfg.integbed.intersect(samplecfg.ltr5, wa=True, F=1, nonamecheck=True) if samplecfg.ltr5 else 0
    if idx:
        ltr5 = samplecfg.integdf.loc[idx.to_dataframe()['name'], :]
        ltr5.drop(columns=drop_cols).to_csv(mydir + '/homologous_recombination/ltr5_integration_' + fn + '.txt', sep='\t')
    else:
        ltr5 = pd.DataFrame(columns=samplecfg.integdf.columns)
    # in LTR3'
    idx = samplecfg.integbed.intersect(samplecfg.ltr3, wa=True, F=1, nonamecheck=True) if samplecfg.ltr3 else 0

    if idx:
        ltr3 = samplecfg.integdf.loc[idx.to_dataframe()['name'], :]
        ltr3.drop(columns=drop_cols).to_csv(mydir + '/homologous_recombination/ltr3_integration_' + fn + '.txt', sep='\t')
    else:
        ltr3 = pd.DataFrame(columns=samplecfg.integdf.columns)
    # in solo-LTR
    idx = samplecfg.integbed.intersect(samplecfg.sololtr, wa=True, F=1, nonamecheck=True) if samplecfg.sololtr else 0

    if idx:
        sololtr = samplecfg.integdf.loc[idx.to_dataframe()['name'], :]
        sololtr.drop(columns=drop_cols).to_csv(mydir + '/homologous_recombination/sololtr_integration_' + fn + '.txt', sep='\t')
    else:
        sololtr = pd.DataFrame(columns=samplecfg.integdf.columns)
    # in exclude
    if exclude:
        idx = samplecfg.integbed.intersect(exclude, wa=True, F=1, nonamecheck=True)
        if idx:
            excludepos = samplecfg.integdf.loc[idx.to_dataframe()['name'], :]
            excludepos.drop(columns=drop_cols).to_csv(mydir + '/excluded/excluded_integration_' + fn + '.txt', sep='\t')
        else:
            excludepos = pd.DataFrame(columns=samplecfg.integdf.columns)
    else:
        excludepos = pd.DataFrame(columns=samplecfg.integdf.columns)
    return trueint, ltr5, ltr3, sololtr, excludepos



def make_location(integ, samplecfg, tsdshiftCDS, cdsdf, legacy_mode):
    """ Run the integration row in Location class and
    return the location-formated row"""
    print('integloc')
    integloc = Location(integ, samplecfg, tsdshiftCDS, cdsdf, legacy_mode)
    print('upORF')
    upORF = integloc.inorf_start if integloc.closestorf['dist'][0] == 0 else integloc.uporf
    print('dnORF')
    dnORF = integloc.inorf_end if integloc.closestorf['dist'][0] == 0 else integloc.dnorf
    row = {'ID': integloc.ID,
           'chr': integloc.chro,
           'position': integloc.position,
           'strand': integloc.strand,
           'indpt' : integloc.indpt if samplecfg.SN else np.nan,
           'dupl' : integloc.dupl,
           'uporf' : upORF['orfname'][0],
           'upstart' : upORF['orfstart'][0],
           'upend' : upORF['orfend'][0],
           'upstrand' : upORF['orfstrand'][0],
           'updist' : upORF['dist'][0],
           'dnorf' : dnORF['orfname'][0],
           'dnstart' : dnORF['orfstart'][0],
           'dnend' : dnORF['orfend'][0],
           'dnstrand' : dnORF['orfstrand'][0],
           'dndist' : dnORF['dist'][0],
           'location' : integloc.loc,
           'intergenic_type' : integloc.intergenic,
           'distance' : integloc.distance}
    return row


def collapse(df):
    """ Collapse by upstream ORF and aggregate columns"""
    agg_dict = {'uporf':'first',
                'upstart':'first',
                'upend':'first',
                'upstrand':'first',
                'dnorf':'first',
                'dnstart':'first',
                'dnend':'first',
                'dnstrand':'first',
                'intergenic_type':'first',
                'SSP':'count',
                'indpt':'sum',
                'dupl':'sum'}
    df['SSP'] = 1
    cdf = df.groupby('uporf').agg(agg_dict)
    return cdf

def make_intervals(location, samplecfg):
    """ Parse location file into ORF vs. intergenics,
    sum independent and duplicate number of integrations per interval,
    merge with the full list of intergenics or ORFs"""
    # calculate sum of integrations per interval
    orf = location.loc[location['location'].isin(['IN_start', 'IN_end'])]
    orf = collapse(orf)[['uporf', 'upstart', 'upend', 'upstrand', 'SSP', 'indpt', 'dupl']]
    intergenic = location.loc[location['location'].isin(['UP', 'DOWN'])]
    intergenic = collapse(intergenic)
    # make full list of ORFs and intergenic regions
    fullorf = samplecfg.cds.to_dataframe()[['name', 'chrom', 'start', 'end', 'strand']]
    fullintergenic = pd.concat(
        [fullorf.loc[:len(fullorf)-2,], fullorf.loc[1:,].reset_index().add_suffix('_dn')], axis=1)
    fullintergenic['intergenic'] = fullintergenic.apply(lambda x: get_intergtype(x['strand'], x['strand_dn']), axis=1)
    # merge full list and integration list
    orf = pd.merge(fullorf, orf[['SSP', 'indpt', 'dupl']],
                   left_on='name',
                   right_on=orf.index,
                   how='left').fillna(0)
    orf[['SSP', 'indpt', 'dupl']] = orf[['SSP', 'indpt', 'dupl']].astype(int)
    intergenic = pd.merge(fullintergenic, intergenic[['SSP', 'indpt', 'dupl']],
                   left_on='name',
                   right_on=intergenic.index,
                   how='left').fillna(0)
    intergenic[['SSP', 'indpt', 'dupl']] = intergenic[['SSP', 'indpt', 'dupl']].astype(int)
    intergenic = intergenic.drop(columns=['index_dn'])
    return orf, intergenic

def count_loc(df, featcol, featname, countcol):
    """ Count the integrations mapped to a given feature
    featcol: column name where the feature is given
    featname: name of the feature (UP, DOWN, ...)
    countcol: column name to count
    """
    df['SSP'] = 1
    df = df[df[featcol].isin(featname)][countcol].sum() if featcol else df[countcol].sum()
    return df

def make_counter(location, ltr5, ltr3, sololtr, excludepos, samplecfg, counter):
    countcols = ['SSP', 'indpt', 'dupl'] if samplecfg.SN else ['SSP', 'dupl']
    for countcol in countcols:
        counter.loc['unfiltered_integration', samplecfg.sample+'_'+countcol] = count_loc(
                                samplecfg.integdf, False, False, countcol)
        counter.loc['true_integration', samplecfg.sample+'_'+countcol] = count_loc(
                                location, False, False, countcol)
        counter.loc['LTR5_integration', samplecfg.sample+'_'+countcol] = count_loc(
                                ltr5, False, False, countcol)
        counter.loc['LTR3_integration', samplecfg.sample+'_'+countcol] = count_loc(
                                ltr3, False, False, countcol)
        counter.loc['solo_integration', samplecfg.sample+'_'+countcol] = count_loc(
                                sololtr, False, False, countcol)
        counter.loc['excluded_position', samplecfg.sample+'_'+countcol] = count_loc(
                                excludepos, False, False, countcol)
        counter.loc['total intergenic', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'location', ['UP', 'DOWN'], countcol)
        counter.loc['upstream intergenic', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'location', ['UP'], countcol)
        counter.loc['downstream intergenic', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'location', ['DOWN'], countcol)
        counter.loc['tandem intergenic', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'intergenic_type', ['Tandem'], countcol)
        counter.loc['divergent intergenic', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'intergenic_type', ['Divergent'], countcol)
        counter.loc['convergent intergenic', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'intergenic_type', ['Convergent'], countcol)
        counter.loc['total in ORF', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'location', ['IN_start', 'IN_end'], countcol)
        counter.loc['start of ORF', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'location', ['IN_start'], countcol)
        counter.loc['end of ORF', samplecfg.sample+'_'+countcol] = count_loc(
                                location, 'location', ['IN_end'], countcol)
    return counter


def main():
    for sample in samples:
        counter_cols = list()
        for t in ['SSP', 'indpt', 'dupl']:
            counter_cols.append(sample + '_' + t)
        counter = pd.DataFrame(columns = counter_cols)

        # sample-specific config
        samplecfg = SampleConfig(config, sample, args)
        fn = samplecfg.sample + '.' +  config['genomevs'][config['genome']]
        # distances and location calculated after taking into account the TSD 
        # at the downstream coordinate of intervals. Need to shift end coordinate of CDSs
        cdsdf = samplecfg.cds.to_dataframe()
        tsdshiftCDS = tsd_shift(cdsdf, samplecfg.tsd)

        if os.path.exists(samplecfg.integfn):
            # filtering of positions matching homologous recombination
            trueint, ltr5, ltr3, sololtr, excludepos = filter_homrecomb(samplecfg, config, fn)
            # maping relative to ORFs
            truedicts = trueint.to_dict('records')
            locdicts = [make_location(x, samplecfg, tsdshiftCDS, cdsdf, legacy_mode) for x in truedicts]
            location = pd.DataFrame.from_dict(locdicts)
            orf, intergenic = make_intervals(location, samplecfg)
            # count integrations mapped to the different categories
            counter = make_counter(location, ltr5, ltr3, sololtr, excludepos, samplecfg, counter)
            if (samplecfg.SN == 0):
                location = location.drop(columns=['indpt'])
                orf = orf.drop(columns=['indpt'])
                intergenic = intergenic.drop(columns=['indpt'])
                counter = counter.drop(columns=[sample+'_indpt'])
            # write output files
            location.drop(columns=['SSP']).to_csv(mydir + '/location_' + fn + '.txt', sep='\t', index=False)
            orf.to_csv(mydir + '/ORF_' + fn + '.txt', sep='\t', index=False)
            intergenic.to_csv(mydir + '/intergenic_' + fn + '.txt', sep='\t', index=False)
        else:
            counter.loc['Warning', sample+'_SSP'] = 'file does not exist: ' + samplecfg.integfn

        counter.to_csv(logdir + '/location_' + config['name'] + '_' + fn + '.log.txt', sep='\t')

main()

