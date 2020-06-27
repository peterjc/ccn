Algorithm to find balanced equivalence relations and lattices
=============================================================

Copyright 2010-2012 by Hiroko Kamei & Peter Cock.  All rights reserved.
This script was provided as a supplementary file for the manuscript:

  Hiroko Kamei and Peter Cock (2013) "Computation of Balanced Equivalence
  Relations and their Lattice for a Coupled Cell Network", SIAM Journal on
  Applied Dynamical Systems (SIADS) 12(1), pp. 352-382.
  http://dx.doi.org/10.1137/100819795
  http://arxiv.org/abs/1211.6334

You are welcome to use and modify this code provided this copyright notice
is retained (see the open source license for full details), however we
request you cite this manuscript in any scientific publications using it.

This code and any public updates to it are available on GitHub, at:
https://github.com/peterjc/ccn


History
-------

| Version | Released   | Description                                         |
|---------|------------|-----------------------------------------------------|
| v1.0.5  | 2019-08-10 | - Adopted black code style.                         |
|         |            | - Using flake8 to check style etc.                  |
| v1.0.4  | 2019-08-09 | - Testing Python 3.5 to 3.8 only                    |
|         |            | - Dropped Python 2 support.                         |
|         |            | - Fixed casting 'int64' to 'uint8' on recent numpy. |
| v1.0.3  | 2015-09-28 | - Minor coding style updates.                       |
|         |            | - Can use pydot_ng as Python 3 pydot alternative.   |
|         |            | - Fix for graphical output under Python 3.          |
| v1.0.2  | 2014-12-02 | - Made the Python script executable; fixed typo.    |
| v1.0.1  | 2014-05-10 | - Dropped support for Python 2.4 and 2.5            |
|         |            | - Supported Python 2.6, 2.7 and 3.2, 3.3, 3.4.      |
| v1.0.0  | 2012-11-04 | - Supported Python 2.4, 2.5, 2.6, 2.7 directly.     |
|         |            | - Worked on Python 3 via 2to3 conversion script.    |

Requirements
------------

This script requires Python. Python is pre-installed on Apple Mac OS X, and on
Linux distributions is easily installable via the package manager. For Windows
use one of the installers available from http://www.python.org

While it originally supported older versions of Python, this script has been
updated to work unmodified on Python 2.6, 2.7, and Python 3.2, 3.3, 3.4 (and
is expected to work on future Python 3 releases).

It requires NumPy (Numerical Python), available from http://numpy.scipy.org
This is used for array support, in particular matrix multiplication. This
should work with NumPy 1.6 or later.

In order to draw the graphs or lattices, it requires the free tool GraphViz
available from http://www.graphviz.org and a small python library to call
GraphViz called pydot from https://github.com/erocarrera/pydot which in turn
requires pyparsing from http://sourceforge.net/projects/pyparsing/

Currently pydot has not been updated to work on Python 3, instead you can
install pydot_ng https://github.com/pydot/pydot-ng


Usage
-----

At the command prompt (terminal window), assuming the Python script has been
saved as graphs.py in the current directory, on Mac or Linux you would run
this code by typing:

    python graphs.py

On Windows, assuming you have installed Python 3.7 in the default location,
you would use the following in the command prompt (also called a "DOS box"
or terminal window):

    C:\Python37\python graphs.py

The script should then run, printing output to the screen, and (assuming the
GraphViz requirements are installed) also generate a number of image files.

The structure of the script is as follows:
 1. An introduction (partly duplicated in this README.md file)
 2. A few key functions
 3. A few classes
 4. Code for self testing
 5. User editable section, as provided this computes the balanced equivalence
    relations and lattices for the graphs in the manuscript. Please copy or edit
    these entries to look at other graphs of interest.


Testing
-------

The graphs.py script includes its own basic self tests, run automatically when
you execute `python graphs.py` on your computer.

[![Build Status](https://secure.travis-ci.org/peterjc/ccn.png?branch=master)](https://travis-ci.org/peterjc/ccn)

These tests are also run automatically whenever changes are uploaded to this
repository using TravisCI, a continuous testing service.


Introduction
------------

The file includes a number of embedded examples (which double as self-tests
using the Python doctest facility). These start with >>> which represents the
Python prompt, and ... for continued lines. If you are not already familiar
with Python and the interactive Python prompt for simplicity we recommend
simply modifying this file and running it as shown above.

As a simple example, consider this regular network with just one edge type
(Graph #7 in the manuscript):

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
    >>> print network7
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

We will now briefly explain some of the details here, in case you want to
modify this.  To obtain the network as an image file, assuming GraphViz etc is
setup, use network7.plot(filename), where the filename can end with ".png",
".pdf", etc. For instance, network7.plot("n7.pdf")

The network object has a method to calculate those cell partitions which are
balanced equivalence relations (balanced colourings):

    >>> for p in network7.partitions():
    ...     print "%r or %s" % (p, cyclic_partition(p))
    [0, 0, 0, 0, 0] or (12345)
    [0, 0, 0, 0, 1] or (1234)(5)
    [0, 1, 0, 1, 1] or (13)(245)
    [0, 1, 0, 1, 2] or (13)(24)(5)
    [0, 1, 2, 3, 1] or (1)(25)(3)(4)
    [0, 1, 2, 3, 4] or (1)(2)(3)(4)(5)

Similar methods allow access to the associated quotient networks instead or as
well. There is also a method which uses the balanced equivalence relations to
build a lattice:

    >>> lattice7 = network7.lattice()
    >>> print lattice7
    0 0 0 0 0 0 (12345)
    1 0 0 0 0 0 (1234)(5)
    1 0 0 0 0 0 (13)(245)
    0 1 1 0 0 0 (13)(24)(5)
    0 0 1 0 0 0 (1)(25)(3)(4)
    0 0 0 1 1 0 (1)(2)(3)(4)(5)

The lattice is printed using a directed graph adjacency matrix (here with
six nodes, listed on the right of the output in cyclic notation) which is
shown on the left as a lower triangular matrix. Rather than being drawn as a
directed graph with arrows, it is conventional to use lattice diagrams with
undirected edges with the directionality of the partition cover relationship
implicit in the vertical placement of the nodes. To obtain the lattice diagram
as an image file, assuming GraphViz etc is setup, use lattice7.plot(filename),
e.g lattice7.plot("n7_lattice.pdf") for a PDF file. Here is a simple text
graphic of this lattice diagram:

    *--------------------------------------*
    |                                      |
    | Rank 1:          (12345)             |
    |                  /     \             |
    |                 /       \            |
    | Rank 2:    (1234)(5)  (13)(245)      |
    |                |     /   |           |
    |                |    /    |           |
    | Rank 3:  (13)(24)(5)     |           |
    |                \         |           |
    |                 \        |           |
    | Rank 4:          \  (1)(25)(3)(4)    |
    |                   \     /            |
    |                    \   /             |
    | Rank 5:       (1)(2)(3)(4)(5)        |
    |                                      |
    *--------------------------------------*

For a multiple edge type example, consider this inhomogeneous version of the
previous network, with two edge types (Graph #3 in the manuscript), where edge
type one (single thickness) arrows run from node 1 to 5, 1 to 2 and 3 to 4,
and edge type two (dotted) arrows run from node 4 to 1 and from 2 to 3.

    *---------------------------*
    |                           |
    |  [5] <--- (1) <... [4]    |
    |            |       /|\    |
    |            |        |     |
    |           \|/       |     |
    |           [2] ...> (3)    |
    |                           |
    *---------------------------*

To enter this topology, two integer edge matrices are required (for the two
edge types):

    >>> network3 = CoupledCellNetwork([[0,0,0,0,0],
    ...                                [1,0,0,0,0],
    ...                                [0,0,0,0,0],
    ...                                [0,0,1,0,0],
    ...                                [1,0,0,0,0]],
    ...                               [[0,0,0,1,0],
    ...                                [0,0,0,0,0],
    ...                                [0,1,0,0,0],
    ...                                [0,0,0,0,0],
    ...                                [0,0,0,0,0]])
    >>> print network3
    (0,0) (0,0) (0,0) (0,1) (0,0) node 1
    (1,0) (0,0) (0,0) (0,0) (0,0) node 2
    (0,0) (0,1) (0,0) (0,0) (0,0) node 3
    (0,0) (0,0) (1,0) (0,0) (0,0) node 4
    (1,0) (0,0) (0,0) (0,0) (0,0) node 5

When printed as shown above, the network adjacency matrix is represented as a
single combined matrix where each element shows the edge types, e.g. (1,0) for
one edge of the first type, no edges of the second type. Internally however
the data remains as two separate matrices. The network object works just the
same as before, for example notice how there are only 4 balanced equivalence
relations (compared to the regular network #7 used above with 6):

    >>> for p in network3.partitions():
    ...     print "%r or %s" % (p, cyclic_partition(p))
    [0, 1, 0, 1, 1] or (13)(245)
    [0, 1, 0, 1, 2] or (13)(24)(5)
    [0, 1, 2, 3, 1] or (1)(25)(3)(4)
    [0, 1, 2, 3, 4] or (1)(2)(3)(4)(5)

Taking the second of these partitions gives a three cell quotient network:

    >>> print network3.quotient([0, 1, 0, 1, 2])
    (0,0) (0,1) (0,0) node 1+3
    (1,0) (0,0) (0,0) node 2+4
    (1,0) (0,0) (0,0) node 5

And using the balanced equivalence relations to build the lattice:

    >>> lattice = network3.lattice()
    >>> print lattice
    0 0 0 0 (13)(245)
    1 0 0 0 (13)(24)(5)
    1 0 0 0 (1)(25)(3)(4)
    0 1 1 0 (1)(2)(3)(4)(5)

The lattice is an undirected graph (here with four nodes, listed on the right
of the output in cyclic notation) which is represented on the left as a lower
triangular matrix. Graphically:

    *--------------------------------------*
    |                                      |
    | Rank 1:    (no nodes at rank one)    |
    |                                      |
    |                                      |
    | Rank 2:        (13)(245)             |
    |                 /    \               |
    |                /      \              |
    | Rank 3:  (13)(24)(5)   \             |
    |                \        \            |
    |                 \        \           |
    | Rank 4:          \   (1)(25)(3)(4)   |
    |                   \      /           |
    |                    \    /            |
    | Rank 5:       (1)(2)(3)(4)(5)        |
    |                                      |
    *--------------------------------------*

The idea is you can edit the examples in last section of the graphs.py file to
run this program on particular networks of interest. In the long term if the
tool is extended, restructuring this into a typical Python library would be
sensible. For now however, a single self contained Python file was simplest.
