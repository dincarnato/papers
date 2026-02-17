# diffShape

## Description
The diffShape tool enables comparing two experimental conditions, to identify bases having a reactivity fold change exceeding a user-defined threshold (__base__ analysis mode), or regions having whose reactivity profiles have a correlation lower than a user-defined threshold (__region__ analysis mode).


## Requirements
The pipeline relies upon the following components:

1. [RNA Framework](https://github.com/dincarnato/RNAFramework) (v2.9.3 or greater)
2. [BEDTools](https://bedtools.readthedocs.io/en/latest/) (v2.31.1 or greater) [optional]

Make sure that the `lib/` folder of RNA Framework is added to the `$PERL5LIB` environment variable:

```bash
$ export PERL5LIB=$PERL5LIB:/path/to/RNAFramework/lib
```

## Running the tool
The diffShape tool must be run as it follows:

```bash
$ diffShape <mode> [params] -c ctrl_1,...,ctrl_n -t treat_1,...,treat_n
```

where:

- `mode` must be either __base__ or __region__
- `-c` and `-t` are comma-separated lists of folders of RNA Framework-compliant XML files containing per-transcript normalized reactivities, respectively for the control and treatment conditions (where each folder represents a replicate experiment)

A full list of parameters is available via `-h` (or `--help`).<br/>
Two *key* parameters are:

- `--compareAll`: by default diffShape performs a paired analysis, therefore, if, for example, 2 control experiments were provided (`-c ctrl_1,ctrl_2`), the tool will then expect 2 treatment experiments (`-t treat_1,treat_2`), and it will compare them in a paired fashion (so, *ctrl_1* vs. *treat_1* and *ctrl_2* vs. *treat_2*). However, when this parameter is specified, diffShape won't expect anymore an equal number of control and treatment experiments, and will perform an all vs. all comparison (in this case: *ctrl_1* vs. *treat_1*, *ctrl_1* vs. *treat_2*, *ctrl_2* vs. *treat_1* and *ctrl_2* vs. *treat_2*)
- `--minComparisons`: specifies the minimum fraction of comparisons that must pass the user-defined thresholds for a base/region to be reported. For example, if 3 control and 3 treatment experiments are passed, along with the `--compareAll` parameter, and `--minComparisons` is set to __0.8__, the tool will require 8 out of 9 comparisons to pass the thresholds.

<br/>

### Base analysis mode
In base analysis mode, the tool selects bases that are covered across all experiments and calculates the reactivity fold change between the control and treatment conditions.

For a base to be considered differentially reactive between the two conditions, it must have a log<sub>2</sub> fold change &ge; `--upDiffThresh` and &le; `--downDiffThresh`, and the reactivity difference between the two conditions must be &ge; `--minDiff`. Therefore, if, for example, `--downDiffThresh` is set to __-1__, `--minDiff` is set to __0.1__, and the reactivity in the control is __0.1__ and in the treatment is __0.05__, the __log<sub>2</sub>(0.05 / 0.1)__ will exceed the `--downDiffThresh` threshold, but the base will still not be reported because the reactivity difference between the two conditions is &lt; `--minDiff`.

Similarly, for a base to be considered unchanged between the two conditions, it must have a log<sub>2</sub> fold change &lt; `--upSimThresh` and &gt; `--downSimThresh`, and the reactivity difference between the two conditions must be &lt; `--minDiff`.

<br/>

### Region analysis mode
In region analysis mode, the tool selects regions of sizes specified via `--winLengths`, slid in `--winOffset` increments, having at least `--minBasesInWin` informative bases across all experiments. 

For a region to be considered differentially structured between the two conditions, the correlation of their reactivity profiles must be &ge; `--diffThresh`, while if it is &le; `--simThresh` the region is considered to be unchanged. 

If `--simThresh` or `--diffThresh` are not specified, the tool will estimate them automatically. Briefly, the distributions of inter-replicate (within the same condition) correlations will be calculated across all `--winLengths` values, and for each window length the `--simThresh` and `--diffThresh` will be respectively set to the median and 25<sup>th</sup> percentile of the distribution.

This automatic threshold estimation requires that at least two replicates were provided for at least one experimental condition. While this procedure is a bit more computationally intensive, it allows estimating window length-specific thresholds.

After differential (or unchanged) regions have been identified across all window lengths, they are merged into a final non-redundant set.