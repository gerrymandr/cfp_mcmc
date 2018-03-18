chain v1.0
==========

This program runs the sqrt(ep) test for redistriciting Markov Chains on the Congressional Districtings of the Commonwealth of Pennsylvania.

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-generate-toc again -->
**Table of Contents**

- [chain v1.0](#chain-v10)
    - [Setup](#setup)
    - [Getting Started](#getting-started)
        - [Running a sample district](#running-a-sample-district)
        - [Polsby-Popper Compactness](#polsby-popper-compactness)
    - [Pesentation](#pesentation)
        - [Generating SVG maps](#generating-svg-maps)
        - [Generating histogram plots](#generating-histogram-plots)

<!-- markdown-toc end -->


## Setup

To install, type
make

then
sudo make install


## Getting Started

The program requires as input the included file `CurrentRep.txt`, which includes the graph and statistical data for the 9059 precincts.  Party affiliation in this file is taken from the 2010 Senate race between Sestak and Toomey.


### Running a sample district

A simple run of the program is then
    chain -f CurrentRep.txt -n 22 --variance --median_mean

This runs the chain for 2^22 steps, with no constraints other than that districts are contiguous and simply connected, and with the default population difference threshold (2%).  It outputs outlier and significance with respect to the variance and median_mean metrics.

Further options can be explored by running
    chain --help


### Polsby-Popper Compactness

For example, to constrain the average inverse Polsby-Popper compactness at 160 and preserve counties not divided by the current districting, one would run:
    chain -f CurrentRep.txt -n 22 --variance --median_mean --counties --L1-compactness 160


Note that not all choices allowed by the program are equally reasonable.  For example, metrics based on the efficiency gap and seat count are insensitive to small changes when the number of districts is small (say, <50).


## Pesentation

### Generating SVG maps

The included file input.svg is used by the program when outputting svg files.


### Generating histogram plots

Start with running the program with the `--histogram` 
    chain -f CurrentRep.txt -n 22 --efficiency_gap --histogram --L1-compactness 160

And record the initial efficiency gap as computed by the program. Then run use the plotting script provided in the `scripts` folder as follows:
    python scripts/plot.py eg_hist.txt [initial efficiency_gap] eg_plot.png


This program accompanies the paper
"Assessing significance in a Markov chain without mixing"
by
Maria Chikina, mchikina@pitt.edu
Alan Frieze, alan@random.math.cmu.edu
Wesley Pegden, wes@math.cmu.edu
