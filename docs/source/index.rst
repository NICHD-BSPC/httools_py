=============
Documentation
=============

HTtools
=======

HTtools maps and quantifies retrotransposon integrations.
Sequence reads matching the end of the retrotransposon LTR are trimmed
of the LTR end and screened against sequences indicating non
integration-specific fragments (PBS sequence, LTR-circle, donor plasmid,
retrotransposon internal sequences and restriction enzyme recognition
site). Reads are then mapped to genome and sorted according
to uniquely mapped, mapping to multiple positions and mapping to
homologous recombination sites. Uniquely mapped reads are positionned
relative to ORF.

The full processing pipeline is included in a Snakefile and can be
configured using a config yaml file. Upon success, an html summary file
is generated. It contains interactive histograms of integration genomic
distribution, distribution relative to ORFs, heatmap of the most targeted
regions and links to the output TSV files.

Modification of the config.yaml file and rerun of the pipeline will only
trigger the rerun of affected rules.

Each individual step scripts can also be run independently by command
line.

A set of *S. pombe* genome references and parameters are included in HTtools by
default. Other genomes can be added to the ``database`` directory and specifying
the paths in the config.yaml file. See section :ref:`Adding custom genome, annotation
and homologous recombination references <adding-custom-genome-annotation-and-homologous-recombination-references>`.



Four versions of *S. pombe* genome are currently included: 

- 2008 release

- 2012 release

- 2008 release plus the sequence of the expression plasmid pHL2882

- 2012 release plus the sequence of the expression plasmid pHL2882

Experimental setups with and without Serial Number are supported, as
well as integrase-frameshift and native integrase.

Parameters that are typical of the Tf1 retrotransposon are used in the
``test/config/config_*.yaml`` configuration files.

For long-term reproducibility, it is recommended to set up a copy of the HTtools
directory per project and to use an environment manager such as conda 
(`bioconda <https://bioconda.github.io/>`__)
to install the dependencies required. We recommend creating such an environment
in each project directory. This allows the versions of all installed packages to
be maintained independently of any other projects or software on the system.


HTtools Github repository
-----------------------------

The Github repository can be found at https://github.com/NICHD-BSPC/httools.

Either download the zip file or clone the repository with ``git`` by
navigating to the location where you want to clone HTtools and

::

   git clone git@github.com:NICHD-BSPC/httools.git


Software installation
---------------------

Installation of all dependencies is handled by conda, ensuring
reproducibility, streamlined setup, and no need for root administrator
privileges.

Use `bioconda <https://bioconda.github.io/>`__ to automatically install
software into the working directory without needing admin rights on the
machine.

If you have the Anaconda Python distribution, you already have conda.
Otherwise, install `Miniconda <https://conda.io/miniconda.html>`__.

1. Creating a new conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a top-level environment with the workflow requirements in the
HTtools directory. It needs to be activated any time you’ll be working
with this workflow. It is recommended to use the ``-p`` option to create
the environment in the current directory. This makes it easier to track
the environment used for each project. HTtools has been implemented and tested
in Linux and MacOS environments. Other plateforms are not currently supported.

Navigate to your HTtools directory.

In a Linux system, create the environment:

::

   conda env create -p env/ --file requirements-linux.yaml

In a MacOS system, create the environment:

::

   conda env create -p env/ --file requirements-macos.yaml

2. Activating the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to your HTtools directory and activate with:

::

   conda activate env/

Eventually when you’re done, you can “deactivate”, which removes the
environment location from your $PATH until the next time you activate
it.

::

   conda deactivate

3. Running the tests
~~~~~~~~~~~~~~~~~~~~

The shell script test_snakemake.sh will run HTtools on Serial Number and
non Serial Number test datasets, and compare the output files to the
expected outputs. This ensure reproducibility by verifying that the
current tools and packages versions are giving the expected results.
This needs to be done only once after installation.

::

   bash test_snakemake.sh

If everything went fine, the following message will be displayed

   Test passed successfully

If the results generated in the current environment differ from the
expected results, the following error message will be displayed

   Error(s) found, see test/diff-test-screen-full.log

The differences found between the expected and the newly
generated test results are listed in the file ``test/diff-test-screen-full.log``
for debugging.


Configuring the workflow
------------------------

Create a configuration file ``config.yaml`` in the HTtools project directory.
Two templates ``config-Tf1.yaml`` and ``config-Hermes.yaml`` are included
as examples, containing default parameters typical of Tf1 and Hermes respectively.

Adjust the parameters to match your experiment design.

All sample information and workflow configurations are specified in the
``config.yaml`` file.

The following fields need to be adjusted to individual runs:

-  ``name`` experiment name. Should be unique in the project directory to avoid 
   overwritting of results. All results will be stored in a directory labelled ``name``

-  ``fastq`` list of path(s) to the fastq file(s). Path(s) can be either
   absolute or relative to the config file. Can be a .gz file

-  ``sample`` block: The sample block must be copied for each sample in
   the fastq files (typically for each barcode). It must start with a
   unique name and contains the fields:

   -  ``barcode_start`` position in the sequence reads where the barcode
      starts
   -  ``barcode_length`` lenght of the barcode in nulceotides
   -  ``sequence`` expected sequence, from the barcode (included) to the
      end of the LTR. Note: if a Serial Number is included, it must be
      indicated with ``x``\ s
   -  ``integrase`` whether the integrase was native (``wt``) or
      frameshift (``fs``)
   -  ``lib_design`` whether the sequence reads originate from the
      ``U5`` or the ``U3`` end of the retrotransposon
   -  ``SN_position`` (optional) start position of the Serial Number,
      indicate ‘none’ if no SN was used
   -  ``SN_length`` (optional) length of the Serial Number, indicate
      ‘none’ if no SN was used

.. code-block:: yaml

    sample:
        # sample block ----------------------------------------------
        BC3498full:
            barcode_start: 1
            barcode_length: 4
            sequence: CTCACCGCAGTTGATGCATAGGAAGCCxxxxxxxxCAAACTGCGTAGCTAACA
            integrase: wt
            lib_design: U5
            SN_position: 28
            SN_length: 8
        # sample block ----------------------------------------------


-  ``genome`` genome built. Current available options are:

   -  ``1``: 2008 release
   -  ``2``: 2012 release
   -  ``3``: 2008 release plus the sequence of the expression plasmid
      pHL2882
   -  ``4``: 2012 release plus the sequence of the expression plasmid
      pHL2882

-  ``generate_uncollapsed`` whether to output (``True``) or not
   (``False``) fastas of trimmed sequence reads corresponding to the
   integration files. Sequences are trimmed after the end of the LTR and
   are replicated as many times as there were duplicate sequence reads.

-  ``exclude``  positions to exclude, in the format
   chromosome_coordinate_orientation, i.e. ``chr1_240580_-``

   Those positions will be screened out from the true_integrations
   and written in ``data/{name}/location/excluded/`` for reference.

   Indicate 'none' if no position to exclude


Advanced parameters include legacy_mode (see section :ref:`legacy_mode changes <legacy-mode-changes>`
for details), reference sequences used for screening, blast parameters, and are also specified in
the ``config.yaml`` file. Those parameters do not typically need to be modified between experiements
as long as the experimental design remains identical. See the section :ref:`Default advanced 
parameters <default-advanced-parameters>`, as well as the example files located in ``test/config/``
for more details.

Indicate `none` in a filtering step parameter to skip this filtering step.


Running the workflow
--------------------

The workflow performs the following tasks:

-  screening of fastq files for non specific sequence reads
-  mapping of the screened reads to the reference genome using ``blast``
-  filtering of the blast results for uniquely mapped reads
-  positioning of the insertions relative to ORFs and quantification
-  plotting of results and creation of summary html file
-  (optional) creation of fasta files containing reads that correspond
   to the integration files


Since HTtools is based on Snakemake, the entire workflow can be executed on a single machine,
submitted to a cluster, or run on cloud platforms (see the `Snakemake <https://snakemake.readthedocs.io/>`__
documentation for details on these execution methods).

Running on a local system
~~~~~~~~~~~~~~~~~~~~~~~~

HTtools should be run on a system with at least (X CPU, Y RAM) due to the computational complexity.

With the environment activated, navigate to the HTtools directory and run the workflow:

::

   snakemake --config fn=config.yaml -R `snakemake --config fn=config.yaml --list-params-changes` --cores=1

Notes:

-  ``--config fn=config.yaml`` indicates the location of your configuration file. This is assuming a file named
    ``config.yaml`` in the HTtools directory. This is a requirement argument.
-  the command ``snakemake --config fn=config.yaml --list-params-changes`` lists the files affected by any parameter
   changes done in the ``config.yaml`` file since the last snakemake execution. ``-R`` triggers the rules that produce
   those files, effectively re-processing and updating any result file dependent of the changed parameters.
-  ``--cores=1`` sets the number of cores used by the workflow to 1. ``--cores=1`` will work on any system; optionally adjust the
   number of cores according to your system's specifications for optimized speed.
-  log and error messages are indicated within the ``Snakefile.log``

Upon success, results can be found in the directory ``data/{name}`` where ``name`` is the experiment name provided in
the ``config.yaml`` file. See section :ref:`Output files of interest <output-files-of-interest>`
for details.

An error is raised and the workflow is aborted when a sample does not return any read.
This is generally due to an error in the sequences specified in the ``config.yaml`` file.
A modified fastqscreen log file ``data/logs/fastq_screen_{name}_{sample}.error.txt`` is generated and contains
the number of reads passing / blocked by each of the sequence filters for debugging.

Running on a SLURM cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~

Optionally HTtools can be run on a cluster. A wrapper file is included for running on a SLURM cluster. Other
types of clusters are not currently supported.

::

    sbatch --cpus-per-task=4 scripts/WRAPPER_SLURM config.yaml

Notes:

-  adjust the number of ``--cpus-per-task`` to your system's specifications.
-  when running parallel jobs, log and error messages are not indicated within the ``Snakefile.log``
file but can rather be found in ``logs/{rule.name}.{jobID}.e``.


Running individual scripts
--------------------------

Alternatively, scripts for the individual steps can be run
independently. See individual scripts code for usage.

This can be useful for example to position the multimatch integrations
relative to ORFs. In this example, a multimatch integration file is
processed through the location step. From the HTtools directory:

::

   python scripts/location.py --integration path/to/data/{name}/filterblast/integration_multimatch_file.txt
   --config path/to/config.yaml

then the output
``path/to/data/{name}/location/location_multimatch_file.txt`` can be
processed through the ORFmap step. From the HTtools directory:

::

   R -e "rmarkdown::render('scripts/results.Rmd',output_file='../wanted/path/to/results.html', params=list(configfn='../path/to/config.yaml'))"

(please note the ``../`` in the output and params arguments, the paths
must be relative to the results.Rmd file)


Adding custom genome, annotation and homologous recombination references
------------------------------------------------------------------------

The pipeline contains by default a set of S. pombe releases. Adding new references can be done
by following the steps below.

Create a custom genome database from a reference fasta file using the tool ``makeblastdb``
from the NCBI BLAST+ tool suite ([Camacho_et_al.,2009]_). ``makeblastdb`` is included in the environment.

::

    makeblastdb -in {genome.fasta} -out {genome.fas} -dbtype nucl -logfile logfile.txt

A BED6-formated file can be used as custom annotation file. BED6 contains the columns
chrom, chromStart, chromEnd, name, score, strand. The score is not used by the pipeline and can be
set to any value.

Copy the created *.nhr, *.nin, *.nsq files, fasta and annotation files
to the directory ``HTtools/database/{new_database_name}``.

Update the paths to ``genomedb`` and ``genomecds`` in the Advanced parameters section of the
YAML configuration file accordingly.

A custom version of a retrotransposon preexisting insertions can be used to detect possible homologous
recombination. Prepare BED6-formated files corresponding to the 3' terminal repeat outmost coordinate (U5),
to the 5' terminal repeat outmost coordinates (U3) and to single repeats (solo-LTR) that originated from
excision of a retrotransposon. Note that the outmost coordinate corresponds to the 3' extremity if the library
was sequenced from U5 or to the 5' extremity if the library was sequenced from U3. Copy those files to
the directory ``HTtools/database/{new_preexisiting_coordinates_name}``.

Update the paths to ``preexist_ltr`` in the Advanced parameters section of the
YAML configuration file accordingly.

.. [Camacho_et_al.,2009] Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+: architecture and applications.
   BMC Bioinformatics 10, 421 (2009). https://doi.org/10.1186/1471-2105-10-421


Output files of interest
------------------------

Output files of interest:

1) ``data/{sample}/results.html``: summary report containing interactive figures and links to all
   result files.
2) ``data/{name}/filterblast/integration_{sample}.txt``: contains the list of integration positions
   with the number of associated sequence reads. If the experiment set
   up includes Serial Number, the last 2 columns indicate the number of
   independent integration events and the number of sequence reads
   respectively.
3) ``data/{name}/location/true_integration_{sample}.txt``: integrations minus the positions matching
   homologous recombination sites and optionnaly the positions to exclude.
4) ``data/{name}/location/homol-recomb_{sample}.txt``: potential homologous recombination events
   filtered out from integration_{sample}.txt
5) ``data/{name}/location/ORF_{sample}.txt``: lists the ORFs and the corresponding number of
   integrations.
6) ``data/{name}/location/intergenic_{sample}.txt``: lists the intergenic regions and the
   corresponding number of integrations.
7) ``data/{name}/location/location_{sample}.txt``: integration positions with assignment to ORF
   or intergenic region.
8) ``data/{name}/ORFmap/ORFmap_{sample}.txt``: table summarizing the % of integration within
   intervals upstream, downstream and within ORFs.
9) ``data/{name}/logs/log_*.txt``: summary of sequence read and integration numbers.

legacy_mode changes
-------------------

When ``legacy_mode`` is set to ``True`` in the config.yaml, the pipeline
follows the behavior of the HTtools perl scripts suite [Esnault_et_al_2019]_ on which
HTtools_py was based.

.. [Esnault_et_al_2019] Esnault C., Lee M., Ham C, Levin L. Transposable element insertions
   in fission yeast drive adaptation to environmental stress. https://doi.org/10.1101/gr.239699.118


fastqscreen
~~~~~~~~~~~

The perl scripts ``screen_illumina_Tf1_sequence-1.0.pl`` and
``screen_illumina_Tf1_SN_sequence-2.0.pl`` screened out sequences with
>= 2 mismatches to end of LTR, or non-specific sequences. This should
have been > 2 to allow up to 2 mismatches. ``legacy_mode=False`` allows
up to 2 mismatches.

filblast
~~~~~~~~

To determine whether the read is multimapped or uniquely mapped, the
perl version compares all the matches, and assign to multi only if all
the matches are within the threshold. It seems more appropriate to at
first only looks at the top 2. If the best match is within the evalue
threshold of the second best, then assign to multimatch any sequence
within the threshold of the top match. ``legacy_mode=False`` follows the
later.

location
~~~~~~~~

The upstream distances to nearest ORF were off by 6 nucleotides in the
perl scripts. Distances to downstream were correct.
``legacy_mode=False`` fixes this issue.

ORFmap
~~~~~~

The perl script ``ORF_map_v2-nonSN.pl`` was counting the header line as
an integration SSP, thus increasing the total number of SSP by 1. Fixed
with ``legacy_mode=False``.

Notes
-----

.. _fastqscreen-1:

fastqscreen
~~~~~~~~~~~

The sequences characteristic of SpeI incomplete are located ~70 bp from
the begining of the sequence reads. The SpeI incomplete sequence would
partially fall outside of the sequence read when the sequencing length
was 100bp. Longer reads (150bp) are prefered for this reason, although
the 100bp still allow SpeI incomplete correct assignment in most cases.

Sequence distances calculations are using different packages between
perl and python scripts. Out of 10,000 reads, tests showed 100% of
identical assignment between the original perl script and the updated
python version for non SN reads. 0.01% reads were assigned to SpeI
incomplete in in python but not in perl out of 10,000 reads with SN
(100bp reads).

The filtering was sequential in the perl version, and was processed
slightly differently between the SN and non SN version. I.e. SpeI
incomplete is only counted if the sequence was neither categorized as
plasmid, nor ltrcircle in the SN version. The non SN version counts any
SpeI incomplete. This may change the numbers within the filtering
categories but does not affect whether a read is filtered out. This
behavior is conserved in python when ``legacy_mode=False``.

Default advanced parameters
---------------------------

.. code-block:: yaml

    # -----------------------------------------------------------
    # Advanced parameters
    # -----------------------------------------------------------
    # Those parameters do not typically need to be modified.
    # Filters against linker, ltrcircle, plasmid, primary_incomplete, 
    # second_incomplete and pbs are optional. Indicate 'none' to skip
    # those filters.
    legacy_mode: False                          # whether to enable legacy_mode
    length_to_match: 34                         # number of nucleotides to match to reference filtering sequences during fastq screening
    min_length: 14                              # minimum length for trimmed reads to be processed
    allowed_mismatches: 2                       # number of mismatches allowed when screening fastqs for reference filtering sequences
    linker: TAGTCCCTTAAGCGGAG                   # sequence of linker to be filtered out
    ltrcircle:
      U5: TGTCAGCAATACTAGCAGCATGGCTGATACACTA    # sequence of terminal-repeat circle to be filtered out for U5 libraries
      U3: TGTTAGCTACGCAGTTACCATAAACTAAATTCCT    # sequence of terminal-repeat circle to be filtered out for U3 libraries
    plasmid:
      U5: GAAGTAAATGAAATAACGATCAACTTCATATCAA    # sequence of donor plasmid to be filtered out for U5 libraries
      U3: none                                  # sequence of donor plasmid to be filtered out for U3 libraries
    primary_re:
      U5: MseI                                  # name of primary restriction enzyme for U5 libraries
      U3: MseI                                  # name of primary restriction enzyme for U3 libraries
    primary_incomplete:
      U5: TTAA                                  # sequence of primary restriction site to be filtere out for U5 libraries
      U3: TTAA                                  # sequence of primary restriction site to be filtere out for U3 libraries
    second_re:
      U5: SpeI                                  # name of secondary restriction enzyme for U5 libraries
      U3: BspHI                                 # name of secondary restriction enzyme for U3 libraries
    second_incomplete:
      U5: AATTCTTTTCGAGAAAAAGGAATTATTGACTAGT    # sequence of secondary restriction site to be filtere out for U5 libraries
      U3: TTACATTGCACAAGATAAAAATATATCATCATGA    # sequence of secondary restriction site to be filtere out for U3 libraries
    dist_to_second_incomplete:
      U5: 28                                    # distance between end of transposable element end and position of secondary_incomplete sequence for U5 libraries
      U3: 22                                    # distance between end of transposable element end and position of secondary_incomplete sequence for U3 libraries
    pbs:
      U5: ATAACTGAACT                           # primer binding site (PBS) sequence to be filtered out for U5 libraries
      U3: TTGCCCTCCCC                           # primer binding site (PBS) sequence to be filtered out for U3 libraries
    tsd:
      wt: 5                                     # length of target site duplication (TSD) in wild-type integrase context
      infs: 0                                   # length of target site duplication (TSD) in integrase-frameshift context
    blastview: 6                                # view parameter for blast results. Modifying this might interfere with subsequent screening steps
    blastevalue: 0.05                           # evalue threshold for blast
    max_score_diff: 0.0001                      # evalue ratio threshold to assign a read to uniquely mapped or multi-location mapped
    orf_map_interval: 100                       # length of outside-ORF intervals in histograms of distribution relative to ORF
    avg_orf_length: 1500                        # average ORF length; will determine the number of within-ORF intervals
    orf_map_window: 5000                        # span of distance to plot in histograms of distribution relative to ORF
    genomedb:                                   # path to the different versions of genome fasta reference files
      1: database/2007/chr123.fas
      2: database/2012_ASM294v2/chr123.fas
      3: database/2007_with_pHL2882/chr123pHL2882.fas
      4: database/2012_ASM294v2_pHL2882/chr123pHL2882.fas
    genomevs:                                   # name of the different versions of genome fasta reference files
      1: v07str
      2: v12str
      3: v07pHL
      4: v12pHL
    preexist_ltr:                               # path to the annotation files for homologous recombination screening, for U5 and U3 libraries
      U5:
        ltr5: database/LTR_2012_ASM294v2/Tf2_5_LTR.txt
        ltr3: database/LTR_2012_ASM294v2/Tf2_3_LTR.txt
        sololtr: database/LTR_2012_ASM294v2/solo_LTR.txt
      U3:
        ltr5: database/LTR_2012_ASM294v2/Tf2_5_LTR-U3.txt
        ltr3: database/LTR_2012_ASM294v2/Tf2_3_LTR-U3.txt
        sololtr: database/LTR_2012_ASM294v2/solo_LTR-U3.txt
    genomecds:                                  # path to the different versions of annotation reference files
      1: database/2007/cds.txt
      2: database/2012_ASM294v2/cds.txt
      3: database/2007_with_pHL2882/cds.txt
      4: database/2012_ASM294v2_pHL2882/cds.txt
    # List of chromosomes of interest
    # The integration log file will give for infomration purpose the count within each chromosome 
    # in the reference genome, but only the chromosomes from the list below will be included 
    # in the output files integration, intergenic, ORF, location, ORFmap, logoDNA
    chro_listvs:
      1: short_chro_list
      2: full_chro_list
      3: short_chro_list
      4: full_chro_list
    full_chro_list:
      - chr1
      - chr2
      - chr3
      - AB325691
    short_chro_list:
      - chr1
      - chr2
      - chr3



Change log
----------

2020-10-29

- Generalize steps to run HTtools on different plateforms

2020-07-14

- Added capability to run jobs in parallel on HPC
- Screen out a list of positions to exclude
- Plot correlation heatmaps

2020-06-29

-  Full rewrite in python
-  Added a results summary html output
-  Added interactive heat maps of the most targeted intergenic regions
   and most targeted ORFs

2019-10-17

-  Version httools.2.0
-  Fix for U3 workflow the ‘BspHI incomplete screen’, orientation of
   integration, recombination events coordinates, and removed the
   ‘plasmid’ screen
-  Added workflows for integrase-independent experiments (IN-indpt) for
   U5 and U3
-  Filter out sequence reads matching LTR circles
-  Screen out multimatch insertions based on blast e-values rather than
   blast bit scores

2019-08-07

-  Added screen from U3 transposon end
-  Allowing compressed fastq.gz as input file
