#!/usr/bin/env python
r"""Algorithm to find balanced equivalence relations and lattices.

Copyright 2010-2012 by Hiroko Kamei & Peter Cock.  All rights reserved.
Revisions copyright 2012-2020 Peter Cock. All rights reserved.

Version 1.0.0 of this script was provided as a supplementary file for the
manuscript:

    Hiroko Kamei and Peter Cock (2013) "Computation of Balanced Equivalence
    Relations and their Lattice for a Coupled Cell Network", SIAM Journal
    on Applied Dynamical Systems (SIADS) 12(1), pp. 352-382.
    http://dx.doi.org/10.1137/100819795 http://arxiv.org/abs/1211.6334

Version 1.1.0 restructured the code into a conventional Python library
(which can be used via ``import couple_cell_networks`` or similar),
retaining ``graphs.py`` as a runable example as in the original release.

You are welcome to use and modify this code provided this copyright notice
is retained, however we request you cite this manuscript in any scientific
publications using it.

This code and any public updates to it are available on GitHub, at:
https://github.com/peterjc/ccn


History & Requirements
======================

See the README file. This script requires Python, NumPy, and additional
libraries for drawing the graphs and lattices as images.


Usage
=====

At the command prompt (terminal window), assuming this file has been saved
as graphs.py in the current directory, on Mac or Linux you would run this
code by typing:

    python graphs.py

On Windows, assuming you have installed Python 3.7 in the default location,
you would use the following in the command prompt (also called a "DOS box"
or terminal window):

    C:\Python37\python graphs.py

The script should then run, printing output to the screen, and (assuming the
GraphViz requirements are installed) also generate a number of image files.

The structure of this script is as follows:

    1. This introduction
    2. Import the relevant libraries
    3. Example usage

As provided this computes the balanced equivalence relations and lattices
for the graphs in the manuscript. Please edit a copy of this script to look
at other graphs of interest.


Introduction
============

The file includes a number of embedded examples (which double as self-tests
using the Python doctest facility). These start with >>> which represents the
Python prompt, and ... for continued lines. If you are not already familiar
with Python and the interactive Python prompt for simplicity we recommend
simply modifying this file and running it as shown above.

As a simple example, consider this regular network with just one edge type
(Graph #7 in the manuscript)::

    *---------------------------*
    |                           |
    |  (5) <--- (1) <--- (4)    |
    |            |       /|\    |
    |            |        |     |
    |           \|/       |     |
    |           (2) ---> (3)    |
    |                           |
    *---------------------------*

To enter this topology, just one integer edge matrix is required, which can
be typed as a nested list of integers:

    >>> network7 = CoupledCellNetwork([[0,0,0,1,0],
    ...                                [1,0,0,0,0],
    ...                                [0,1,0,0,0],
    ...                                [0,0,1,0,0],
    ...                                [1,0,0,0,0]])
    >>> print(network7)
    0 0 0 1 0 node 1
    1 0 0 0 0 node 2
    0 1 0 0 0 node 3
    0 0 1 0 0 node 4
    1 0 0 0 0 node 5

The printed output shows the adjacency matrix on the left (5 by 5) with the
node captions on the right.

To handle typical usage there is a go() function provided which will print the
matrix, draw it (if GraphViz is setup), compute the lattice, print it, and
draw it. Just pass a CoupledCellNetwork object to the go function with a name
(a string), e.g. go(network7, "n7")

The idea is you can edit the examples below in a copy of this file to run the
algorithm on your particular networks of interest.
"""

import sys

import numpy as np

from coupled_cell_networks import __version__
from coupled_cell_networks import CoupledCellNetwork
from coupled_cell_networks import go

if "-v" in sys.argv or "--version" in sys.argv:
    print("Coupled Cell Networks " + __version__)
    sys.exit(0)

##########################################################
#
# As supplied, this first shows the seven graphs from the manuscript Table 2.1
#
##########################################################

# Network one has three edge types,
e1 = [[0, 0, 0, 0], [1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]]
e2 = [[0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 1], [0, 0, 0, 0]]
e3 = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
go(CoupledCellNetwork(e1, e2, e3), "n1")

##########################################################

e1 = [[0, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]]
e2 = [[0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
e3 = [[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
e4 = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 1, 0]]
go(CoupledCellNetwork(e1, e2, e3, e4), "n2")

##########################################################

e1 = [
    [0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [1, 0, 0, 0, 0],
]
e2 = [
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
]
go(CoupledCellNetwork(e1, e2), "n3")

##########################################################

e1 = [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]]
e2 = [[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
e3 = [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]
go(CoupledCellNetwork(e1, e2, e3), "n4")

##########################################################

e1 = [
    [0, 1, 0, 1, 0],
    [1, 0, 0, 1, 0],
    [0, 0, 1, 0, 1],
    [1, 1, 0, 0, 0],
    [0, 0, 1, 0, 1],
]
e2 = [
    [0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0],
]
go(CoupledCellNetwork(e1, e2), "n5")

##########################################################

e1 = [
    [0, 0, 0, 1, 1],
    [1, 0, 0, 1, 0],
    [1, 0, 0, 0, 1],
    [1, 1, 0, 0, 0],
    [1, 0, 1, 0, 0],
]
go(CoupledCellNetwork(e1), "n6")

##########################################################

e1 = [
    [0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [1, 0, 0, 0, 0],
]
go(CoupledCellNetwork(e1), "n7")

##########################################################

e1 = [
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
]
e2 = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 1],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
]
# go(CoupledCellNetwork(e1,e2), "n8")

##########################################################
# fmt: off
# This should match Figure 3(b) in Belykh and Hasler (2011) Chaos,
# Mesoscale and clusters of synchrony in networks of bursting neurons.
#                      1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3
#    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
e1 = [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],  # 1  # noqa: E231
    [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],  # 2  # noqa: E231
    [0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # 3  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0],  # 4  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],  # 5  # noqa: E231
    [0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],  # 6  # noqa: E231
    [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],  # 7  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1],  # 8  # noqa: E231
    [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],  # 9  # noqa: E231
    [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],  # 10 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],  # 11 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1],  # 12 # noqa: E231
    [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # 13 # noqa: E231
    [1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # 14 # noqa: E231
    [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0],  # 15 # noqa: E231
    [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0],  # 16 # noqa: E231
    [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0],  # 17 # noqa: E231
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],  # 18 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0],  # 19 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0],  # 20 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],  # 21 # noqa: E231
    [0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],  # 22 # noqa: E231
    [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],  # 23 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # 24 # noqa: E231
    [0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # 25 # noqa: E231
    [0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # 26 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],  # 27 # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0],  # 28 # noqa: E231
    [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],  # 29 # noqa: E231
    [0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # 30 # noqa: E231
]
# Double check the row sums match the number of edges expected,
# Node                1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3
#   1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
assert [
    2,2,2,3,2,3,2,2,2,2,2,3,1,3,3,3,3,2,2,2,2,3,2,2,2,2,1,2,2,2,  # noqa: E231
] == [sum(row) for row in e1], [sum(row) for row in e1]
# fmt: off
assert np.all(np.array(e1) == np.array(e1).T), "Not symmetric"
go(CoupledCellNetwork(e1), "n9")

# Modified version, remove edge between 4 and 17 (3 and 16 in Python
# counting), replace with edge joining 13 and 27 (12 and 26 in Python
# counting):
e1[3][16] = 0
e1[16][3] = 0
e1[12][26] = 1
e1[26][12] = 1
assert np.all(np.array(e1) == np.array(e1).T), "Not symmetric"
go(CoupledCellNetwork(e1), "n9-modified")

##################################################################

# This should match Figure 5 in in Belykh and Hasler (2011) Chaos,
# Mesoscale and clusters of synchrony in networks of bursting neurons.
e1 = [
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 1, 1, 0, 0, 0, 0],
    [1, 0, 0, 1, 1, 1, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 1, 1, 1, 1],
    [0, 1, 1, 0, 0, 0, 1, 1, 1, 1],
    [0, 1, 1, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
]
assert np.all(np.array(e1) == np.array(e1).T), "Not symmetric"

print("")
print("=" * 50)
print("Warning this next one prints out a lot... (big lattice!)")
go(CoupledCellNetwork(e1), "n10")

######################################################

# This row of six nodes is based on Zhang et al (2001) Physical Review E, vol 63
e1 = [
    [1, 1, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0],
    [0, 1, 0, 1, 0, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 0, 0, 1, 0, 1],
    [0, 0, 0, 0, 1, 1],
]
assert np.all(np.array(e1) == np.array(e1).T), "Not symmetric"
go(CoupledCellNetwork(e1), "n11")

# This row of seven nodes is based on Zhang et al (2001) Physical Review E, vol 63
e1 = [
    [1, 1, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 0, 0],
    [0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 1, 0, 1],
    [0, 0, 0, 0, 0, 1, 1],
]
assert np.all(np.array(e1) == np.array(e1).T), "Not symmetric"
go(CoupledCellNetwork(e1), "n12")

######################################################

# Casado (2003) Figure 2, repeated in the review paper
# Ibarz et al (2010) as Figure 61(a).
e1 = [
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 1, 0, 0, 1],
    [0, 1, 1, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
]
go(CoupledCellNetwork(e1), "Casado_F2_original")

# This is Fig. 2 from Rabinovich et al (2001),
# FitzHugh-Nagumo models in neuronal dynamics.
# Same as Ibaraz example, but with an edge from node 1 to 5.
e1 = [
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 1, 1, 0],
    [1, 0, 0, 1, 0, 1, 0, 0, 1],
    [0, 1, 1, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
]
go(CoupledCellNetwork(e1), "Rabinovich_F2_original")

# Modified version of Casado (2003) example to use two edge types,
# those from n3 and n8 excitatory (edge type two), the rest are
# inhibitory (edge type one).
e1 = [
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 1, 0, 0, 1],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]
e2 = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
]
go(CoupledCellNetwork(e1, e2), "Casado_F2_mod_2e_i38")

# Modified from Rabinovich et al (2000), not Ibaraz et al. (2010)
# Two edge types, outputs from n3 and n8 excitatory (edge type 2),
# rest are inhibitory (edge type 1)
# i.e. Like the above, but with an output from one to five.
e1 = [
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 1, 0, 0],
    [1, 0, 0, 1, 0, 1, 0, 0, 1],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]
e2 = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0],
]
go(CoupledCellNetwork(e1, e2), "Rabinovich_F2_mod_2e_i38")

######################################################

# To check Belykh's conjecture (2.8) in:
# Belykh et al. Hierarchy and stability of partially synchronous
# oscillations of diffusively coupled dynamical systems.
# Phys Rev E Stat Phys Plasmas Fluids Relat Interdiscip Topics
# (2000) vol. 62 (5 Pt A) pp. 6332-45
e1 = [
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
]
# i.e. make_chain(15)
# Finding the full lattice (not just the top node) took 3.5 days.
go(CoupledCellNetwork(e1), "Belykh_2000_chain15", top_only=True)

######################################################

# From: 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16   To:
e1 = [
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n1 yellow bottom
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # n2 blue bottom
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n3 red bottom-right
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n4 green bottom-righ
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # n5 yellow right
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n6 blue right
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n7 red top-right
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n8 green top-right
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # n9 yellow top
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],  # n10 blue top
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],  # n11 red top-left
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],  # n12 green top-left
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],  # n13 yellow left
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # n14 blue left
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],  # n15 red bottom-left
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],  # n16 green bottom-left
]
assert [sum(row) for row in e1] == [1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1], [
    sum(row) for row in e1
]
# fmt: off
assert [sum(e1[i][j] for i in range(16)) for j in range(16)] == [
    1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1
], [sum(e1[i][j] for j in range(16)) for i in range(16)]
# fmt: on
go(CoupledCellNetwork(e1), "Lehnert_2010_Syncline_Poster")

# Modified to include self coupling to make this a regular netork,
# From:
#    1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16   To:
e1 = [
    [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n1 yellow bottom
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # n2 blue bottom
    [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n3 red bottom-right
    [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n4 green bottom-right
    ######################################################
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # n5 yellow right
    [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n6 blue right
    [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n7 red top-right
    [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # n8 green top-right
    ######################################################
    [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],  # n9 yellow top
    [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],  # n10 blue top
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],  # n11 red top-left
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],  # n12 green top-left
    ######################################################
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],  # n13 yellow left
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # n14 blue left
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0],  # n15 red bottom-left
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0],  # n16 green bottom-left
]
for row in e1:
    assert sum(row) == 2
# fmt: off
assert [sum(e1[i][j] for i in range(16)) for j in range(16)] == [
    1, 2, 2, 3, 3, 2, 3, 2, 2, 1, 1, 2, 2, 3, 2, 1,
], [sum(e1[i][j] for j in range(16)) for i in range(16)]
# fmt: on
go(CoupledCellNetwork(e1), "Lehnert_2010_Syncline_Poster_reg", top_only=True)

###################################################################

# 30 node example from here:
# http://jasss.soc.surrey.ac.uk/12/2/3.html
# fmt: off
e1 = [
    [0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1],  # noqa: E231
    [1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],  # noqa: E231
    [1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0],  # noqa: E231
    [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,1],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0],  # noqa: E231
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0],  # noqa: E231
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1],  # noqa: E231
    [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0],  # noqa: E231
]
# fmt: on
go(CoupledCellNetwork(e1), "sw30")

###################################################################

# 24 node example
e1 = [
    [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
    [1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
    [1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1],
    [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
]
go(CoupledCellNetwork(e1), "sw24")

###################################################################

# 22 node example
e1 = [
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
]
go(CoupledCellNetwork(e1), "tangled22")

print("Done")
