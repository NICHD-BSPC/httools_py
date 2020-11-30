#!/usr/bin/env python


import os
from pathlib import Path
import ruamel.yaml
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import pybedtools as pbt

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
        bed = pbt.BedTool.from_dataframe(df[col_names])
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
        bed = pbt.BedTool(fn)
        extbed = pd.DataFrame(columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
        for chro in self.chromsizes.keys():
            extbed.loc[chro + '_start'] = [chro, 1, 1, chro+'_start', 0, 'na']
            extbed.loc[chro + '_end'] = [chro, self.chromsizes[chro]+1, self.chromsizes[chro]+1, chro+'_end', 0, 'na']
        extbed = pbt.BedTool.from_dataframe(extbed)
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
    def __init__(self, integ, samplecfg, tsdshiftCDS):
        #self.ID = integ['id']
        #self.chro = integ['chr']
        #self.position = integ['position']
        #self.strand = integ['strand']
        #self.dupl = integ['dupl']
        #self.indpt = integ['indpt']
        #self.integbed = self.integ_to_bed(self.ID, self.chro, self.position, self.strand)
        # location determined relative to the shifted CDS, then the coordinates are shifted back to original
        self.uporf = self.map_interg(integ, tsdshiftCDS, fu=True)
        self.dnorf = self.map_interg(integ, tsdshiftCDS, fd=True)
        self.closestorf = self.map_interg(integ, tsdshiftCDS)
        self.inorf_start = self.map_inorf(integ, tsdshiftCDS, fu=True)
        self.inorf_end = self.map_inorf(integ, tsdshiftCDS, fd=True)
        self.loc = self.get_loc(self.uporf, self.dnorf, self.closestorf, self.inorf_start, self.inorf_end)
        self.intergenic = self.get_intergtype(self.closestorf, self.uporf, self.dnorf)

        self.distance = self.get_distance(self.inorf_start, self.inorf_end, self.closestorf, self.uporf, self.dnorf)

    def map_interg(self, integ, tsdshiftedCDS, fu=False, fd=False):
        #dfcols = [0'chrom', 1'start', 2'end', 3'name', 4'score', 5'strand',
        #             6'orfchr', 7'orfstart', 8'orfend', 9'orfname', 10'orfscore', 11'orfstrand',
        #             12'dist']
        #           chr1    12693   12693   SSP_chr1_12693_+    0   +
        #           chr1    12158   12985   SPAC212.08c         na  +
        #           0
        # if ties in closest, take only the first occurance
        bed = integ.closest(tsdshiftedCDS, D = "ref",
                              fu=fu, fd=fd,
                              k=1, nonamecheck=True,
                              t='first')
        bed = bed.to_dataframe(header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                       'orf_chrom', 'orf_start', 'orf_end', 'orf_name', 'orf_score', 'orf_strand', 'distance'])
        return bed.set_index('name')

    def map_inorf(self, integ, orf, fu=False, fd=False):
        # create 2 intervals of 1bp per ORF for start and end, for each orf
        # NOTE: added -1 and +1 to start and end because the distances were shifted by 1bp compared to
        # the original Perl HTtools
        with open('inters.txt', 'w') as fout:
            for y in orf:
                # next part of line is ORFs made up from start coordinates
                start = str(y.start-1) if y.start > 1 else str(y.start)
                line = [
                    y.chrom,
                    start, start,
                    y.name,
                    y.score,
                    y.strand
                ]
                fout.write('\t'.join(line) + '\n')
                # next part of line is ORFs made up from end coordinates
                line = [
                    y.chrom,
                    str(y.end+1), str(y.end+1),
                    y.name,
                    y.score,
                    y.strand
                ]
                fout.write('\t'.join(line) + '\n')
        orfbed = pbt.BedTool('inters.txt').sort()
        return self.map_interg(integ, orfbed, fu=fu, fd=fd)

    def get_loc(self, uporf, dnorf, closestorf, inorf_start, inorf_end):
        locdf = pd.DataFrame(columns = ['direction'], index = uporf.index)
        locdf['direction'] = np.where(
            # in ORF, closer to start, + strand ORF
            (closestorf['distance'] == 0) & (abs(inorf_start['distance']) <= abs(inorf_end['distance'])) \
                & (closestorf['orf_strand'] == '+'),
            'IN_start', locdf['direction'])
        locdf['direction'] = np.where(
            # in ORF, closer to start, - strand ORF
            (closestorf['distance'] == 0) & (abs(inorf_start['distance']) <= abs(inorf_end['distance'])) \
                & (closestorf['orf_strand'] == '-'),
            'IN_end', locdf['direction'])
        locdf['direction'] = np.where(
            # in ORF, closer to end, + strand ORF
            (closestorf['distance'] == 0) & (abs(inorf_start['distance']) > abs(inorf_end['distance'])) \
                & (closestorf['orf_strand'] == '+'),
            'IN_end', locdf['direction'])
        locdf['direction'] = np.where(
            # in ORF, closer to end, - strand ORF
            (closestorf['distance'] == 0) & (abs(inorf_start['distance']) > abs(inorf_end['distance'])) \
                & (closestorf['orf_strand'] == '-'),
            'IN_start', locdf['direction'])
        locdf['direction'] = np.where(
            # outside ORF, closer to uporf, + strand ORF
            (closestorf['distance'] != 0) & (abs(uporf['distance']) <= abs(dnorf['distance'])) \
                & (uporf['orf_strand'] == '+'),
            'DOWN', locdf['direction'])
        locdf['direction'] = np.where(
            # outside ORF, closer to uporf, - strand ORF
            (closestorf['distance'] != 0) & (abs(uporf['distance']) <= abs(dnorf['distance'])) \
                & (uporf['orf_strand'] == '-'),
            'UP', locdf['direction'])
        locdf['direction'] = np.where(
            # outside ORF, closer to dnorf, + strand ORF
            (closestorf['distance'] != 0) & (abs(uporf['distance']) > abs(dnorf['distance'])) \
                & (dnorf['orf_strand'] == '+'),
            'UP', locdf['direction'])
        locdf['direction'] = np.where(
            # outside ORF, closer to dnorf, - strand ORF
            (closestorf['distance'] != 0) & (abs(uporf['distance']) > abs(dnorf['distance'])) \
                & (dnorf['orf_strand'] == '-'),
            'DOWN', locdf['direction'])
#        if (int(closestorf[12]) == 0 and abs(int(inorf_start[12])) <= abs(int(inorf_end[12]))):
 #           direction = 'IN_start' if closestorf[11] == '+' else 'IN_end'
 #       elif (int(closestorf[12]) == 0 and abs(int(inorf_start[12])) > abs(int(inorf_end[12]))):
 #           direction = 'IN_end' if closestorf[11] == '+' else 'IN_start'
#        elif (abs(int(uporf[12])) <= abs(int(dnorf[12]))):
#            direction = 'UP' if uporf[11] == '-' else 'DOWN'
 #       elif (abs(int(uporf[12])) > abs(int(dnorf[12]))):
 #           direction = 'DOWN' if dnorf[11] == '-' else 'UP'

        # special case if start or end of chromosome, to match perl version
        locdf['direction'] = np.where(
            # uporf strand is na, dnorf strand is -
            (uporf['orf_strand'] == 'na') & (dnorf['orf_strand'] == '-'),
            'DOWN', locdf['direction'])
        locdf['direction'] = np.where(
            # uporf strand is na, dnorf strand is +
            (uporf['orf_strand'] == 'na') & (dnorf['orf_strand'] == '+'),
            'UP', locdf['direction'])
        locdf['direction'] = np.where(
            # dnorf strand is na, uporf strand is -
            (dnorf['orf_strand'] == 'na') & (uporf['orf_strand'] == '-'),
            'UP', locdf['direction'])
        locdf['direction'] = np.where(
            # dnorf strand is na, uporf strand is +
            (dnorf['orf_strand'] == 'na') & (uporf['orf_strand'] == '+'),
            'DOWN', locdf['direction'])
        locdf['direction'] = np.where(
            # dnorf and uporf strands are na
            (dnorf['orf_strand'] == 'na') & (uporf['orf_strand'] == 'na'),
            'IN_end', locdf['direction'])
#        if (uporf[11] == 'na'):
#            direction = 'DOWN' if dnorf[11] == '-' else 'UP'
#        if (dnorf[11] == 'na'):
#            direction = 'UP' if uporf[11] == '-' else 'DOWN'
#        if (dnorf[11] == 'na' and uporf[11] == 'na'):
#            direction = 'IN_end'
        return(locdf)

    def get_intergtype(self, closestorf, uporf, dnorf):
        inter = pd.DataFrame(columns = ['type'], index = uporf.index)
        inter['type'] = np.where((uporf['orf_strand'] == '+') & (dnorf['orf_strand'] == '+'),
                                 'Tandem', inter['type'])
        inter['type'] = np.where((uporf['orf_strand'] == '-') & (dnorf['orf_strand'] == '-'),
                                 'Tandem', inter['type'])
        inter['type'] = np.where((uporf['orf_strand'] == '+') & (dnorf['orf_strand'] == '-'),
                                 'Convergent', inter['type'])
        inter['type'] = np.where((uporf['orf_strand'] == '-') & (dnorf['orf_strand'] == '+'),
                                 'Divergent', inter['type'])
        inter['type'] = np.where((uporf['orf_strand'] == 'na') | (dnorf['orf_strand'] == 'na'),
                                 '', inter['type'])
        inter['type'] = np.where((closestorf['distance'] == 0),
                                 '--', inter['type'])
        return inter

    def get_distance(self, inorf_start, inorf_end, closestorf, uporf, dnorf):
        dist = pd.DataFrame(columns = ['dist'], index = uporf.index)
        dist['dist'] = np.where(
            (closestorf['distance'] == 0) & (abs(inorf_start['distance']) <= (inorf_end['distance'])),
            abs(inorf_start['distance']), dist['dist'])
        dist['dist'] = np.where(
            (closestorf['distance'] == 0) & (abs(inorf_start['distance']) > (inorf_end['distance'])),
            abs(inorf_end['distance']), dist['dist'])
        dist['dist'] = np.where(
            (closestorf['distance'] != 0) & (abs(uporf['distance']) <= (dnorf['distance'])),
            abs(uporf['distance']), dist['dist'])
        dist['dist'] = np.where(
            (closestorf['distance'] != 0) & (abs(uporf['distance']) > (dnorf['distance'])),
            abs(dnorf['distance']), dist['dist'])
        return dist



def tsd_reverse(shifteddf, cdsdf):
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
    return(pbt.BedTool.from_dataframe(df))

def integ_to_bed(integ):
    integ.loc[:,'score'] = 0
    cols = ['chr', 'position', 'position', 'id', 'score', 'strand']
    return pbt.BedTool.from_dataframe(integ[cols])

def exclude_to_bed(exclude):
    """ Takes list of positions to exclude and return BedTool object"""
    df = pd.DataFrame(exclude, columns=['name'])
    df[['chr', 'start', 'strand']] = df.name.str.split('_', expand=True)
    df['score'] = 1
    df = df[['chr', 'start', 'start', 'name', 'score', 'strand']]
    return(pbt.BedTool.from_dataframe(df))


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



def make_location(integ, integdf, samplecfg, tsdshiftCDS, cdsdf):
    """ Run the integration row in Location class and
    return the location-formated row"""
    integloc = Location(integ, samplecfg, tsdshiftCDS)
    integdf = integdf.loc[:,['id', 'chr', 'position', 'strand', 'indpt', 'dupl']]
    integdf['indpt'] = integdf['indpt'] if samplecfg.SN else np.nan
    integdf['uporf'] = np.where(integloc.closestorf['distance'] == 0,
                                integloc.inorf_start['orf_name'],
                                integloc.uporf['orf_name'])
    # ORF coordinates are still shifted, need to shift back
    integdf['upstart'] = [int(cdsdf[cdsdf['name'] == x].loc[:,'start']) for x in integdf['uporf']]
    integdf['upend'] = [int(cdsdf[cdsdf['name'] == x].loc[:,'end']) for x in integdf['uporf']]
    integdf['upstrand'] = np.where(integloc.closestorf['distance'] == 0,
                                integloc.inorf_start['orf_strand'],
                                integloc.uporf['orf_strand'])
    integdf['updist'] = np.where(integloc.closestorf['distance'] == 0,
                                integloc.inorf_start['distance'],
                                integloc.uporf['distance'])
    integdf['dnorf'] = np.where(integloc.closestorf['distance'] == 0,
                                integloc.inorf_end['orf_name'],
                                integloc.dnorf['orf_name'])
    integdf['dnstart'] = [int(cdsdf[cdsdf['name'] == x].loc[:,'start']) for x in integdf['dnorf']]
    integdf['dnend'] = [int(cdsdf[cdsdf['name'] == x].loc[:,'end']) for x in integdf['dnorf']]
    integdf['dnstrand'] = np.where(integloc.closestorf['distance'] == 0,
                                integloc.inorf_end['orf_strand'],
                                integloc.dnorf['orf_strand'])
    integdf['dndist'] = np.where(integloc.closestorf['distance'] == 0,
                                integloc.inorf_end['distance'],
                                integloc.dnorf['distance'])
    integdf['location'] = integloc.loc['direction']
    integdf['intergenic_type'] = integloc.intergenic['type']
    integdf['distance'] = integloc.distance['dist']

#    row = {'ID': integloc.ID,
#           'chr': integloc.chro,
#           'position': integloc.position,
#           'strand': integloc.strand,
#           'indpt' : integloc.indpt if samplecfg.SN else np.nan,
#           'dupl' : integloc.dupl,
#           # ORF coordinates are still shifted, need to shift back
#           'uporf' : upORF[9],
#           'upstart' : int(cdsdf[cdsdf['name'] == upORF[9]]['start']),
#           'upend' : int(cdsdf[cdsdf['name'] == upORF[9]]['end']),
#           'upstrand' : upORF[11],
#           'updist' : upORF[12],
#           'dnorf' : dnORF[9],
#           'dnstart' : int(cdsdf[cdsdf['name'] == dnORF[9]]['start']),
#           'dnend' : int(cdsdf[cdsdf['name'] == dnORF[9]]['end']),
#           'dnstrand' : dnORF[11],
#           'dndist' : dnORF[12],
#           'location' : integloc.loc,
#           'intergenic_type' : integloc.intergenic,
#           'distance' : integloc.distance}
    return integdf


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


def check_cds(cds):
    # check gene names are unique and convert to dataframe
    cdsdf = cds.to_dataframe()
    if len(cdsdf['name']) > len(cdsdf['name'].unique()):
        raise ValueError('Invalid CDS file. Gene names must be unique:\n' +
                         '\n'.join(cdsdf[cdsdf['name'].duplicated()]['name']) + '\n')
    else:
        return cdsdf

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
        cdsdf = check_cds(samplecfg.cds)
        tsdshiftCDS = tsd_shift(cdsdf, samplecfg.tsd)

        if os.path.exists(samplecfg.integfn):
            # filtering of positions matching homologous recombination
            trueint, ltr5, ltr3, sololtr, excludepos = filter_homrecomb(samplecfg, config, fn)
            # maping relative to ORFs
            truebed = integ_to_bed(trueint)
            location = make_location(truebed, trueint, samplecfg, tsdshiftCDS, cdsdf)
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

