import sys
import os
from textwrap import dedent
import ruamel.yaml
import tempfile
import pandas as pd


fn = config['fn'] if config.get('fn') else 'config.yaml'
yaml = ruamel.yaml.YAML()
config = yaml.load(open(fn))

# ----------------------------------------------------------------------------
# RULES
# ----------------------------------------------------------------------------


final_targets = [
        expand('data/{name}/logs/fastq_screen_{name}_{sample}.log.txt',
               name=config['name'], sample=config['sample'].keys()),
        expand('data/{name}/logs/{logtype}_{name}_{sample}.{version}.log.txt',
            logtype=['blast', 'integration', 'location'], name=config['name'],
               sample=config['sample'].keys(), version=config['genomevs'][config['genome']]),
        expand('data/{name}/fastqscreen/screen_{sample}.fa',
               sample=config['sample'].keys(), name=config['name']),
        expand('data/{name}/blast/blast_{sample}.{version}.txt', name=config['name'],
                   sample=config['sample'].keys(), version=config['genomevs'][config['genome']]),
        expand('data/{name}/filblast/integration_{sample}.{version}.txt', name=config['name'],
                   sample=config['sample'].keys(), version=config['genomevs'][config['genome']]),
        expand('data/{name}/location/{fntype}_{sample}.{version}.txt', name=config['name'],
                   fntype=['true_integration', 'location', 'intergenic', 'ORF'],
                   sample=config['sample'].keys(), version=config['genomevs'][config['genome']]),
        'data/{name}/results.html'.format(name=config['name'])
]

# conditional uncollapsed fasta target
if (config['generate_uncollapsed'] == True):
    unclp = expand('data/{name}/uncollapsed_fastas/uncollapsed_true_integration_{sample}.{version}.fa', name=config['name'],
                   sample=config['sample'].keys(), version=config['genomevs'][config['genome']])
    final_targets.append(unclp)


rule targets:
    """
    Final targets to create
    """
    input:
        final_targets


rule fastqscreen:
    """
    Filtering of fastq files
    """
    input:
        fq=[fq for fq in config['fastq']],
        script=os.path.join(workflow.basedir, "scripts/fastqscreen.py")
    params:
        p1=config['name'],
        p2=config['fastq'],
        p3=config['sample'],
        p4=config['legacy_mode'],
        p5=config['length_to_match'],
        p6=config['min_length'],
        p7=config['allowed_mismatches'],
        p8=config['linker'],
        p9=config['ltrcircle'],
        p10=config['plasmid'],
        p11=config['primary_re'],
        p12=config['primary_incomplete'],
        p13=config['second_re'],
        p14=config['second_incomplete'],
        p15=config['dist_to_second_incomplete'],
        p16=config['pbs'],
        p17=config['tsd']
    output:
        log='data/{name}/logs/fastq_screen_{name}_{sample}.log.txt',
        fns='data/{name}/fastqscreen/screen_{sample}.fa',
        fnid='data/{name}/fastqscreen/screen_{sample}_ids.txt'
    run:
        shell('python {input.script} --config {fn} --sample {wildcards.sample}')

rule blast:
    """
    Mapping of filtered fasta files
    """
    input:
        fns='data/{name}/fastqscreen/screen_{sample}.fa',
        log='data/{name}/logs/fastq_screen_{name}_{sample}.log.txt',
        script=os.path.join(workflow.basedir, "scripts/blast.py")
    params:
        p1=config['genome'],
        p2=config['blastview'],
        p3=config['blastevalue'],
        p4=config['genomedb'],
        p5=config['genomevs'],
    output:
        fns='data/{name}/blast/blast_{sample}.{version}.txt',
        log='data/{name}/logs/blast_{name}_{sample}.{version}.log.txt',
    run:
        shell('python {input.script} --config {fn} --sample {wildcards.sample}')

rule filterblast:
    """
    Filtering of blast results
    """
    input:
        fns='data/{name}/blast/blast_{sample}.{version}.txt', 
        log='data/{name}/logs/blast_{name}_{sample}.{version}.log.txt',
        script=os.path.join(workflow.basedir, "scripts/filterblast.py")
    params:
        p1=config['max_score_diff'],
        p2=config['preexist_ltr'],
        p3=config['chro_listvs'],
        p4=config['full_chro_list'],
        p5=config['short_chro_list']
    output:
        fns=expand('data/{{name}}/filblast/{fntype}_{{sample}}.{{version}}.txt',
            fntype=['integration', 'id_mappings']),
        log='data/{name}/logs/integration_{name}_{sample}.{version}.log.txt',
    run:
        shell('python {input.script} --config {fn} --sample {wildcards.sample}')

rule location:
    """
    Mapping integrations relative to ORFs
    """
    input:
        fns='data/{name}/filblast/integration_{sample}.{version}.txt',
        log='data/{name}/logs/integration_{name}_{sample}.{version}.log.txt',
        script=os.path.join(workflow.basedir, "scripts/location.py")
    params:
        p1=config['genomecds'],
        p2=config['exclude']
    output:
        fns=expand('data/{{name}}/location/{fntype}_{{sample}}.{{version}}.txt',
                   fntype=['true_integration', 'location', 'intergenic', 'ORF']),
        log='data/{name}/logs/location_{name}_{sample}.{version}.log.txt',
    run:
        shell('python {input.script} --config {fn} --sample {wildcards.sample}')

rule cat_logs:
    """
    Concatenate sample logs for filterblast and location into experiment logs
    """
    input:
        logint=expand('data/{{name}}/logs/integration_{{name}}_{sample}.{version}.log.txt',
                sample=config['sample'].keys(), version=config['genomevs'][config['genome']]),
        logloc=expand('data/{{name}}/logs/location_{{name}}_{sample}.{version}.log.txt',
                sample=config['sample'].keys(), version=config['genomevs'][config['genome']]),
        logfas=expand('data/{{name}}/logs/fastq_screen_{{name}}_{sample}.log.txt',
                sample=config['sample'].keys(), version=config['genomevs'][config['genome']])
    output:
        fnint='data/{name}/logs/integration_{name}_log.txt',
        fnloc='data/{name}/logs/location_{name}_log.txt',
        fnfas='data/{name}/logs/fastq_screen_{name}_log.txt'
    run:
        df = {}
        for fn in input.logint:
            df[fn] = pd.read_csv(fn, index_col=[0], sep='\t', header=None)
        finaldf = pd.concat(df, axis=1, join='inner').to_csv(output.fnint, sep='\t', header=None)
        df = {}
        for fn in input.logloc:
            df[fn] = pd.read_csv(fn, index_col=[0], sep='\t', header=None)
        finaldf = pd.concat(df, axis=1, join='inner').to_csv(output.fnloc, sep='\t', header=None)
        df = {}
        for fn in input.logfas:
            df[fn] = pd.read_csv(fn, index_col=None, sep='\t', skiprows=3, header=0, dtype=object)
            loghead = pd.read_csv(fn, index_col=None, sep='\t', nrows=2, header=None)
        finaldf = pd.concat(df, axis=0, join='inner')
        loghead.to_csv(output.fnfas, sep='\t', header=None, index=None)
        finaldf.to_csv(output.fnfas, sep='\t', index=None, mode='a')

rule downstream:
    """
    Report and ORFmap
    """
    input:
        fns=expand('data/{{name}}/location/{fntype}_{sample}.{version}.txt',
                   fntype=['true_integration', 'location', 'intergenic', 'ORF'],
                   sample=config['sample'].keys(), version=config['genomevs'][config['genome']]),
        log=expand('data/{{name}}/logs/{logtype}_{{name}}_log.txt',
                    logtype=['fastq_screen', 'integration', 'location']),
        script=os.path.join(workflow.basedir, "scripts/results.Rmd")
    params:
        p1=config['orf_map_interval'],
        p2=config['avg_orf_length'],
        p3=config['orf_map_window']
    output:
        fn='data/{name}/results.html',
    run:
        cmd = "rmarkdown::render('" + input.script + "',output_file='" + os.path.abspath(output.fn) + \
            "', params=list(fn='" \
            + os.path.abspath(fn) + "'))"
        shell('R -e "{cmd}"')

rule uncollapsed:
    """
    Generate fasta file(s) of trimmed sequence reads corresponding
    to the integration file
    Sequences are trimmed after the end of the LTR and replicated
    as many times as there was duplicate sequences.
    """
    input:
        fn1='data/{name}/fastqscreen/screen_{sample}.fa',
        fn2='data/{name}/filblast/id_mappings_{sample}.{version}.txt',
        fn3='data/{name}/location/true_integration_{sample}.{version}.txt',
        script=os.path.join(workflow.basedir, "scripts/uncollapsedfasta.py")
    params:
        p1=config['generate_uncollapsed']
    output:
        fns='data/{name}/uncollapsed_fastas/uncollapsed_true_integration_{sample}.{version}.fa',
    shell:
        'python {input.script} --config {fn} --sample {wildcards.sample}'


