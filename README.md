# TEAnnotator
## Purpose
A pipeline to provide annotations of transposable element (TE) content in a genome, combining the inputs from a variety of annotators (e.g. RepeatMasker, RepeatModeler) into one non-redundant fasta file.
## Requirements
Written in python3.

Requires:

[RepeatMasker](http://repeatmasker.org/RMDownload.html) v4.0.6 or later (and its dependencies)

[RepeatModeler](http://repeatmasker.org/RepeatModeler.html) v1.0.8 or later (and its dependencies)
# Usage
Basic usage is:
```bash
TEAnnotator.py -i input.fas
```
TEAnnotator accepts the following additional arguments (all of which have defaults already set):

-c (number of cores)

-s (taxon model to use for RepeatMasker)

-u (cutoff value for RepeatMasker)

-m (minimum length of repeat annotation)

--lowcomplex (include low complexity repeats, which are screened out by default) 
