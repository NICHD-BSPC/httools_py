#!/bin/bash

# -----------------------------------------------------------------

# SN version legacy_mode
snakemake -s Snakefile --config fn='test/config/config_SN_legacy.yaml' --cores=1

# step 1: filtering of fastqs
diff -rb <( tail -n +3 test/config/data/test_fullfastqSN/logs/fastq_screen_test_fullfastqSN_log.txt) <( tail -n +3 test/expected/SNlegacy/SRR7068454full_screen_log.txt) \
    > test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

# step 2: blast
diff -rbB <(cut -f2-12 test/config/data/test_fullfastqSN/blast/blast_BC3512full.v12str.txt) <(cut -f2-12 test/expected/SNlegacy/blast_BC3512fullv12str.txt) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

# step 3: filter blast results
diff -rbB <(cut -f2,9,13,14,15,16 test/config/data/test_fullfastqSN/filblast/integration_BC3498full.v12str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16,17 test/expected/SNlegacy/integration_BC3498fullv12str.txt | tail -n +2) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log
diff -rbB <(cut -f2,9,13,14,15,16 test/config/data/test_fullfastqSN/filblast/integration_multimatch_BC3498full.v12str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16,17 test/expected/SNlegacy/integration_multimatch_BC3498fullv12str.txt | tail -n +2) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

# step 4a: filter homologous recombination positions
diff -rwB <(cut -f2,9,13,14,15,16 test/config/data/test_fullfastqSN/location/true_integration_BC3498full.v12str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16,17 test/expected/SNlegacy/true_integration_BC3498fullv12str.txt | tail -n +2) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log
diff -rwB <(cut -f2,9,13,14,15,16 test/config/data/test_fullfastqSN/location/homologous_recombination/sololtr_integration_BC3498full.v12str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16,17 test/expected/SNlegacy/homol-recomb-solo-LTR_BC3498fullv12str.txt | tail -n +1) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

# step 4b: position integrations relative to ORFs
diff -rwB <(awk '!/AB325691/' test/config/data/test_fullfastqSN/location/location_BC3498full.v12str.txt | tail -n +2 | cut -f2,3,4,5,6,10,11,16,17 ) \
    <(awk '!/AB325691/' test/expected/SNlegacy/location_true_BC3498fullv12str.txt | tail -n +2 | cut -f2,3,4,5,6,10,11,16,17)  >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

diff -rwB <(awk '!/_start.*_end/' test/config/data/test_fullfastqSN/location/intergenic_BC3498full.v12str.txt | awk '!/_end.*_start/' | tail -n +2 | cut -f12,13,14 ) \
    <(awk '!/start\t|end */' test/expected/SNlegacy/intergenic_true_BC3498fullv12str.txt | cut -f10,11,12 ) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

diff -rwB <(awk '!/_start|_end/' test/config/data/test_fullfastqSN/location/ORF_BC3498full.v12str.txt | tail -n +2 | cut -f1,3,4,5,6,7,8 ) \
    <(awk '!/start\t|end */' test/expected/SNlegacy/ORF_true_BC3498fullv12str.txt ) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

# step 5: ORFmap
diff -rwB <(awk '/BC3498full/' test/config/data/test_fullfastqSN/ORFmap/ORFmap.tsv | cut -f6,7,8 ) \
    <(cut -f2,3,4 test/expected/SNlegacy/ORFmap_BC3498fullv12str.txt | tail -n +2) >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log

# step 6: uncollapsed fastas
diff -rb  test/config/data/test_fullfastqSN/uncollapsed_fastas/uncollapsed_true_integration_BC3498full.v12str.fa test/expected/SNlegacy/uncollapsed_true_integration_BC3498fullv12str.fa \
    >> test/diff-test-legacy-screen-fullSN.log 2>> test/diff-test-legacy-screen-fullSN.log


# test if log has error
if [ -s test/diff-test-legacy-screen-fullSN.log ]
then
    echo "Error(s) found, see test/diff-test-legacy-screen-fullSN.log"
    head test/diff-test-legacy-screen-fullSN.log
    exit 1
else
    echo "Test SN legacy_mode passed successfully"
fi

# -----------------------------------------------------------------

# non SN version legacy_mode
snakemake -s Snakefile --config fn='test/config/config_nonSN_legacy.yaml' --cores=1

# step 1: filter fastq
diff -rb <( tail -n +3 test/config/data/test_fullfastq/logs/fastq_screen_test_fullfastq_log.txt) <( tail -n +3 test/expected/nonSNlegacy/SRR5305121full_screen_log.txt) \
    > test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log

# step 2: blast
diff -rbB <(cut -f2-12 test/config/data/test_fullfastq/blast/blast_BCfullnonSN.v07str.txt) <(cut -f2-12 test/expected/nonSNlegacy/blast_BCfullnonSNv07str.txt) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log

# step 3: filter blast results
diff -rbB <(cut -f2,9,13,14,16 test/config/data/test_fullfastq/filblast/integration_BCfullnonSN.v07str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16 test/expected/nonSNlegacy/integration_BCfullnonSNv07str.txt | tail -n +2) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log
diff -rbB <(cut -f2,9,13,14,16 test/config/data/test_fullfastq/filblast/integration_multimatch_BCfullnonSN.v07str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16 test/expected/nonSNlegacy/integration_multimatch_BCfullnonSNv07str.txt | tail -n +2) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log

# step 4a: filter homologous recombination positions
diff -rwB <(cut -f2,9,13,14,15,16 test/config/data/test_fullfastq/location/true_integration_BCfullnonSN.v07str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16,17 test/expected/nonSNlegacy/true_integration_BCfullnonSNv07str.txt | tail -n +2) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log
diff -rwB <(cut -f2,9,13,14,15,16 test/config/data/test_fullfastq/location/homologous_recombination/sololtr_integration_BCfullnonSN.v07str.txt | tail -n +2) \
    <(cut -f2,9,14,15,16,17 test/expected/nonSNlegacy/homol-recomb-solo-LTR_BCfullnonSNv07str.txt | tail -n +1) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log

# step 4b: position integrations relative to ORFs
diff -rwB <(cut -f2,3,4,5,6,10,15,16,17 test/config/data/test_fullfastq/location/location_BCfullnonSN.v07str.txt | tail -n +2) \
    <(cut -f2,3,4,5,6,10,15,16,17 test/expected/nonSNlegacy/location_true_BCfullnonSNv07str.txt | tail -n +2) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log
diff -rwB <(awk '!/_start.*_end/' test/config/data/test_fullfastq/location/intergenic_BCfullnonSN.v07str.txt | awk '!/_end.*_start/' | tail -n +2 | cut -f12,13) \
    <(awk '!/start\t|end\t/' test/expected/nonSNlegacy/intergenic_true_BCfullnonSNv07str.txt | cut -f10,11) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log
diff -rwB <(awk '!/_start|_end/' test/config/data/test_fullfastq/location/ORF_BCfullnonSN.v07str.txt | tail -n +2 | cut -f1,3,4,5,6,7) \
    <(awk '!/start\t|end\t/' test/expected/nonSNlegacy/ORF_true_BCfullnonSNv07str.txt) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log

# step 5: ORFmap
diff -rwB <(cut -f6,8 test/config/data/test_fullfastq/ORFmap/ORFmap.tsv | tail -n +2) \
    <(cut -f2,3 test/expected/nonSNlegacy/ORFmap_BCfullnonSNv07str.txt | tail -n +2) >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log

# step 6: uncollapsed fastas
diff -rb  test/config/data/test_fullfastq/uncollapsed_fastas/uncollapsed_true_integration_BCfullnonSN.v07str.fa test/expected/nonSNlegacy/uncollapsed_true_integration_BCfullnonSNv07str.fa \
    >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-legacy-screen-full.log

# test if log has error
if [ -s test/diff-test-legacy-screen-full.log ]
then
    echo "Error(s) found, see test/diff-test-legacy-screen-full.log"
    head test/diff-test-legacy-screen-full.log
    exit 1
else
    echo "Test non SN legacy_mode passed successfully"
fi


# -----------------------------------------------------------------

# SN version non legacy_mode
snakemake -s Snakefile --config fn='test/config/config_SN.yaml' --cores=1

# step 1: filtering of fastqs
diff -rb <( tail -n +3 test/config/data/test_nonlegacySN/logs/fastq_screen_test_nonlegacySN_log.txt) <( tail -n +3 test/expected/SN/SRR7068454full_screen_log.txt) \
    > test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log

# step 2: blast
diff -rbB test/config/data/test_nonlegacySN/blast/blast_BC3512full.v12str.txt test/expected/SN/blast_BC3512fullv12str.txt >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log

# step 3: filter blast results
diff -rbB test/config/data/test_nonlegacySN/filblast/integration_BC3498full.v12str.txt \
    test/expected/SN/integration_BC3498fullv12str.txt >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log
diff -rbB test/config/data/test_nonlegacySN/filblast/integration_multimatch_BC3498full.v12str.txt \
     test/expected/SN/integration_multimatch_BC3498fullv12str.txt >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log

# step 4a: filter homologous recombination positions
diff -rwB test/config/data/test_nonlegacySN/location/true_integration_BC3498full.v12str.txt \
    test/expected/SN/true_integration_BC3498fullv12str.txt >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log
diff -rwB test/config/data/test_nonlegacySN/location/homologous_recombination/sololtr_integration_BC3498full.v12str.txt \
    test/expected/SN/homol-recomb-solo-LTR_BC3498fullv12str.txt >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log

# step 4b: position integrations relative to ORFs
diff -rwB test/config/data/test_nonlegacySN/location/location_BC3498full.v12str.txt \
     test/expected/SN/location_true_BC3498fullv12str.txt  >> test/diff-test-legacy-fullSN.log 2>> test/diff-test-screen-fullSN.log

diff -rwB test/config/data/test_nonlegacySN/location/intergenic_BC3498full.v12str.txt \
    test/expected/SN/intergenic_true_BC3498fullv12str.txt  >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log

diff -rwB test/config/data/test_nonlegacySN/location/ORF_BC3498full.v12str.txt \
    test/expected/SN/ORF_true_BC3498fullv12str.txt >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log

# step 5: ORFmap
diff -rwB test/config/data/test_nonlegacySN/ORFmap/ORFmap.tsv \
     test/expected/SN/ORFmap_BC3498fullv12str.txt >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log

# step 6: uncollapsed fastas
diff -rb  test/config/data/test_nonlegacySN/uncollapsed_fastas/uncollapsed_true_integration_BC3498full.v12str.fa test/expected/SN/uncollapsed_true_integration_BC3498fullv12str.fa \
    >> test/diff-test-screen-fullSN.log 2>> test/diff-test-screen-fullSN.log


# test if log has error
if [ -s test/diff-test-screen-fullSN.log ]
then
    echo "Error(s) found, see test/diff-test-screen-fullSN.log"
    head test/diff-test-screen-fullSN.log
    exit 1
else
    echo "Test SN non legacy_mode passed successfully"
fi

# -----------------------------------------------------------------


# non SN version non legacy_mode
snakemake -s Snakefile --config fn='test/config/config_nonSN.yaml' --cores=1

# step 1: filter fastq
diff -rb <( tail -n +3 test/config/data/test_nonlegacy/logs/fastq_screen_test_nonlegacy_log.txt) <( tail -n +3 test/expected/nonSN/SRR5305121full_screen_log.txt) \
    > test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log

# step 2: blast
diff -rbB test/config/data/test_nonlegacy/blast/blast_BCfullnonSN.v07str.txt test/expected/nonSN/blast_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log

# step 3: filter blast results
diff -rbB test/config/data/test_nonlegacy/filblast/integration_BCfullnonSN.v07str.txt \
    test/expected/nonSN/integration_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log
diff -rbB test/config/data/test_nonlegacy/filblast/integration_multimatch_BCfullnonSN.v07str.txt \
    test/expected/nonSN/integration_multimatch_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log

# step 4a: filter homologous recombination positions
diff -rwB test/config/data/test_nonlegacy/location/true_integration_BCfullnonSN.v07str.txt \
    test/expected/nonSN/true_integration_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log
diff -rwB test/config/data/test_nonlegacy/location/homologous_recombination/sololtr_integration_BCfullnonSN.v07str.txt \
    test/expected/nonSN/homol-recomb-solo-LTR_BCfullnonSNv07str.txt >> test/diff-test-legacy-screen-full.log 2>> test/diff-test-screen-full.log

# step 4b: position integrations relative to ORFs
diff -rwB test/config/data/test_nonlegacy/location/location_BCfullnonSN.v07str.txt \
    test/expected/nonSN/location_true_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log
diff -rwB test/config/data/test_nonlegacy/location/intergenic_BCfullnonSN.v07str.txt \
    test/expected/nonSN/intergenic_true_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log
diff -rwB test/config/data/test_nonlegacy/location/ORF_BCfullnonSN.v07str.txt \
    test/expected/nonSN/ORF_true_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log

# step 5: ORFmap
diff -rwB test/config/data/test_nonlegacy/ORFmap/ORFmap.tsv \
    test/expected/nonSN/ORFmap_BCfullnonSNv07str.txt >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log

# step 6: uncollapsed fastas
diff -rb  test/config/data/test_nonlegacy/uncollapsed_fastas/uncollapsed_true_integration_BCfullnonSN.v07str.fa test/expected/nonSN/uncollapsed_true_integration_BCfullnonSNv07str.fa \
    >> test/diff-test-screen-full.log 2>> test/diff-test-screen-full.log

# test if log has error
if [ -s test/diff-test-screen-full.log ]
then
    echo "Error(s) found, see test/diff-test-screen-full.log"
    head test/diff-test-screen-full.log
    exit 1
else
    echo "Test non SN non legacy_mode passed successfully"
fi


