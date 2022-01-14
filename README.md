## MeShClust2

The newest version of MeShClust (v3.0) can be obtained from https://github.com/BioinformaticsToolsmith/Identity.git

This repository is no longer supported.

Release version - 2.3.0

### Requirements
g++ 4.9.1 or later, requires Homebrew on Mac OS X
Compilation using g++ (homebrew) and CMake on Mac OS X see [this link](https://stackoverflow.com/questions/29057437/compile-openmp-programs-with-gcc-compiler-on-os-x-yosemite)

### Linux/Unix compilation
``` sh
mkdir bin && cd bin
cmake ..
make
```

### Citation
If you find this tool helpful, please cite:

[James, Benjamin T. et al. (2018), MeShClust2: Application of alignment-free identity scores in clustering long DNA sequences. bioRxiv, 451278.](https://doi.org/10.1101/451278)

### Usage

    Usage: meshclust2 --id 0.x [OPTIONS] *.fasta

    --id          The most important parameter, --id, controls the identity cutoff of the sequences.
                  Needs to be between 0 and 1.
                  If it is not specified, an identity of 0.9 is used.

    --kmer        decides the size of the kmers. It is by default automatically decided by average sequence
                  length, but if provided, MeShClust can speed up a little by not having to find the largest
                  sequence length. Increasing kmer size can increase accuracy, but increases memory consumption.

    --dump       Run until the classifier is trained, and then dump the weights to the file,
                 default 'weights.txt'. Can be used with --recover to recover the weights
                 instead of re-training.

    --recover    Recover weights for the classifier trained by a previous run which used --dump to dump
                 the weights.

    --list       Instead of specifying files as extra arguments, provide a text file with
                 a list of files. Can use pipes or process substitutions such as "--list <(ls *.fasta) "

    --no-train-list    Same as --list, but these files are not passed to the classifier,
                       e.g. unassembled genomes

    --mut-type   {single, both, nonsingle-typical, nonsingle-all, all-but-reversion, all-but-translocation}
                 changes the mutation generation algorithm. By default, "both" is used, utilizing
                 single point and block mutations. On higher identity data sets, "single", which includes only single point mutations,
                 is preferable. The option "nonsingle-typical" uses only block mutations,
                 disallowing single point mutations. Other options include "all", which includes single,
                 block, and nontypical mutations translocation and reversion.

    --feat       determines the combinations of features to be used. By default, "slow" allows 11
                 combinations to be selected from. "fast" removes 2 slower features from "slow"
                 which include logarithm based features.

    --single-file  Using this option, (no value is needed), each file is treated as a single sequence.
                   If multiple sequences in a file are encountered, they are joined with 50 Ns,
                   and the k-mers are not counted in that region.
                   However, to be most accurate, it is advised to not use these sequences in the
                   training step (for mutations) and instead 1) train using un-joined sequences and
                   use --dump to dump to a file, and 2) use --recover with --single-file for the
                   file list.

    --sample     selects the total number of sequences used for both training and testing.
                 2000 is the default value. That is, --sample 2000 provides 2000 training
                 pairs and 2000 testing pairs.

    --num-templates   selects the number of "template" sequences from which to mutate.
                 For example, if 300 (the default) templates are requested, and the number of
                 "samples" is requested to be 2000 (the default), 300 sequences will be read in
                 and mutated 2000/300 times each to create 2000 semi-synthetic pairs.

    --min-feat   (default 4) sets the minimum feature pairs to be used. If set to 2, at least 2 feature pairs
                 will be used. Recall that features include pairwise combinations of the "feat" option.

    --max-feat   (default 4) sets the maximum feature pairs to be used. Diminishing returns appears quickly,
                 so a very large maximum (>10) is not advised.

    --min-id     (default 0.35) sets the lower bound for mutation identity scores to be calculated.
                 Shouldn't need to be set normally, as lower identites take much longer,
                 especially with single mutations only.

    --datatype   (8,16,32,64) Decides the integer size of the histograms. If not provided,
                 all sequences are read in and counted to ensure the largest k-mer does not
                 overflow. If the provided k-mer is too small, it will overflow.

    --threads    sets the number of threads to be used. By default OpenMP uses the number of available cores
                 on your machine, but this parameter overwrites that.

    --output     specifies the output file, in CD-HIT's CLSTR format, described below:
                 A '>Cluster ' followed by an increasing index designates a cluster.
                 Otherwise, the sequence is printed out.
                 A '*' at the end of a sequence designates the center of the cluster.
                 An example of a small data set:

                 >Cluster 0
                 0       993nt, >seq128 template_6... *
                 >Cluster 1
                 0       1043nt, >seq235 template_10...
                 1       1000nt, >seq216 template_10... *
                 2       1015nt, >seq237 template_10...

    --delta      decides how many clusters are looked around in the final clustering stage.
                 Increasing it creates more accuracy, but takes more time. Default value is 5.

    --iterations specifies how many iterations in the final stage of merging are done until convergence.
                 Default value is 15.

    If the argument is not listed here, it is interpreted as an input (FASTA format) file.


### License

Academic use: The software is provided as-is under the GNU GPLv3.
Any restrictions to use for-profit or non-academics: License needed.
