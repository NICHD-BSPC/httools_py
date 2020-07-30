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
is generated. It contains interactive graphs of integration genomic
distribution, ORF-maps, heatmap of the most targeted regions and links
to the output TSV files.

Modification of the config.yaml file and rerun of the pipeline will only
trigger the rerun of affected rules.

Each individual step scripts can also be run independently by command
line.

*S. pombe* genome references and parameters are included in HTtools by
default. Other genomes can be added to the ``database`` directory and specifying
the paths in the config.yaml file.

Four versions of *S. pombe* genome are currently included: 

- 2008 release

- 2012 release

- 2008 release plus the sequence of the expression plasmid pHL2882

- 2012 release plus the sequence of the expression plasmid pHL2882

Experimental setups with and without Serial Number are supported, as
well as integrase-frameshift and native integrase.

Parameters that are typical of the Tf1 retrotransposon are used in the
``test/config/config_*.yaml`` configuration files.


HTtools Github repository
-----------------------------

The Github repository can be found at https://github.com/NICHD-BSPC/httools_py.

Either download the zip file or clone the repository with ``git`` by
navigating to the location where you want to clone HTtools and

::

   git clone git@github.com:NICHD-BSPC/httools_py.git


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
the environment used for each project.

Navigate to your HTtools directory.

If you are on Linux, create the environment:

::

   conda env create -p env/ --file requirements-linux.yaml

If you are on macOS, create the environment:

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

Note: if you are running on Biowulf, run this in an interactive session.

::

   bash test_snakemake.sh

If everything went fine, the following message will be displayed

   Test passed successfully

If the results generated in the current environment differ from the
expected results, the following error message will be displayed

   Error(s) found, see test/diff-test-screen-full.log

Contact the Bioinformatics core for assistance.

Configuring the workflow
------------------------

All sample information and workflow configurations are indicated in the
config.yaml file.

The following fields need to be adjusted to individual runs:

-  ``name`` experiment name. Should be unique since all results will be
   stored in a directory labelled with that name

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
      indicate ‘na’ if no SN was used
   -  ``SN_length`` (optional) length of the Serial Number, indicate
      ‘na’ if no SN was used

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
   chromosome_coordinate_orientation, i.e. chr1_240580_-

   Those positions will be screened out from the true_integrations
   and written in ``data/{name}/location/excluded/`` for reference.

   Indicate 'na' if no position to exclude


Advanced parameters include legacy_mode, reference sequences used for
screening, blast parameters, and are also specified in the config.yaml
file. Those parameters do not typically need to be modified. See the
exemple files ``test/config/config_nonSN.yaml`` and
``test/config/config_SN.yaml`` for details.

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

With the environment activated, navigate to the ``httools_py`` directory
and run the workflow:

::

   bash path/to/start_HTtools.sh path/to/config.yaml

For exemple, from the HTtools directory, with a config file called
``config_nonSN.yaml`` and located in HTtools test/config directory:

::

   bash start_HTtools.sh test/config/config_nonSN.yaml

The config.yaml file is a required argument.

Upon success, results can be found in your config.yaml file directory
under ``data/{name}`` where ``name`` is the experiment name provided in
the config.yaml file. See section :ref:`Output files of interest<Output files of interest>`
for details.

An error is raised and the workflow is aborted when a sample does not return any read.
This is generally due to an error in the sequences specified in the ``config.yaml`` file.
A modified fastqscreen log file ``data/logs/fastq_screen_{name}_{sample}.error.txt`` is generated and contains
the number of reads passing / blocked by each of the sequence filters for debugging.

The workflow is set up to automatically run up to 4 parallel jobs on
Biowulf HPC. Adjust ``start_HTtools.sh`` to your `$HOSTNAME` if you are running on a different HPC.
When running parallel jobs, the error messages are not indicated within the ``Snakefile.log``
file but can rather be found in ``logs/{rule.name}.{jobID}.e``.


Running individual scripts
--------------------------

Alternatively, scripts for the individual steps can be run
independently. See individual scripts code for usage.

This can be useful for exemple to position the multimatch integrations
relative to ORFs. In this exemple, an multimatch integration file is
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
follows the behavior of the original perl scripts.

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

Change log
----------

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
