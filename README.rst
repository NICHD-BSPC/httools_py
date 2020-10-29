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

Documentation
-------------

Documentation can be found at http://NICHD-BSPC.github.io/httools.
