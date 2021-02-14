#!/usr/bin/env python
r"""Algorithm to find lattices of synchrony subspaces & their reduced lattices.

Copyright 2010-2014 by Hiroko Kamei & Peter J. A. Cock.
Revisions copyright 2014-2022 by Hiroko Kamei, Haibo Ruan, and Peter Cock.

The initial version (v1.0.0) was originally published as a supplementary file
for the manuscript:

  Hiroko Kamei and Peter Cock (2013) "Computation of Balanced Equivalence
  Relations and their Lattice for a Coupled Cell Network", SIAM Journal
  on Applied Dynamical Systems (SIADS) 12(1), pp. 352-382.
  http://dx.doi.org/10.1137/100819795 http://arxiv.org/abs/1211.6334

Public releases of version (v1.0.0) are available on GitHub, at:
https://github.com/peterjc/ccn

The current version (v1.1.0) was extended to construct a reduced lattice, and
the modified script provided as a supplementary file for the manuscript:

  Hiroko Kamei and Haibo Ruan (2021) "Reduced Lattices of Synchrony Subspaces
  and their Indices", submitted to SIAM Journal on Applied Dynamical Systems
  (SIADS). https://arxiv.org/abs/2007.07414

You are welcome to use and modify this code provided this copyright notice
retained, however we request you cite the manuscripts in any scientific
publications using it. It is released under the 3-clause BSD license.


History & Requirements
======================

See the README file. This script requires Python 3, NumPy, and additional
libraries for drawing the graphs and lattices as images.


Usage
=====

At the command prompt (terminal window), assuming this file has been saved
as graphs.py in the current directory, on Mac or Linux you would run this
code by typing::

    python graphs.py

On Windows, assuming you have installed Python 3.7 in the default location,
you would use the following in the command prompt (also called a "DOS box"
or terminal window)::

    C:\Python37\python graphs.py

The script should then run, printing output to the screen, and (assuming the
GraphViz requirements are installed) also generate a number of image files.

The structure of this script is as follows:
(1) This introduction
(2) A few key functions
(3) A few classes
(4) Code for self testing
(5) User editable section, as provided this computes the balanced equivalence
relations and lattices for the graphs in the first manuscript, and lattice
reductions from the second manuscript. Please copy or edit these entries to
look at other graphs of interest.


Introduction
============

The file includes a number of embedded examples (which double as self-tests
using the Python doctest facility). These start with >>> which represents the
Python prompt, and ... for continued lines. If you are not already familiar
with Python and the interactive Python prompt for simplicity we recommend
simply modifying this file and running it as shown above.

As a simple example, consider this regular network with just one edge type::

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

    >>> network1 = CoupledCellNetwork([[0,0,0,1,0],
    ...                                [1,0,0,0,0],
    ...                                [0,1,0,0,0],
    ...                                [0,0,1,0,0],
    ...                                [1,0,0,0,0]])
    >>> print(network1)
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
(a string), e.g. go(network1, "n1")

We will now briefly explain some of the details here, in case you want to
modify this.  To obtain the network as an image file, assuming GraphViz etc is
setup, use network1.plot(filename), where the filename can end with ".png",
".pdf", etc. For instance, network1.plot("n1.pdf")

The network object has a method to calculate those cell partitions which are
balanced equivalence relations (balanced colourings):

    >>> for p in network1.partitions():
    ...     print("%r or %s" % (p, cyclic_partition(p)))
    [0, 0, 0, 0, 0] or (12345)
    [0, 0, 0, 0, 1] or (1234)(5)
    [0, 1, 0, 1, 1] or (13)(245)
    [0, 1, 0, 1, 2] or (13)(24)(5)
    [0, 1, 2, 3, 1] or (1)(25)(3)(4)
    [0, 1, 2, 3, 4] or (1)(2)(3)(4)(5)

Similar methods allow access to the associated quotient networks instead or as
well. There is also a method which uses the balanced equivalence relations to
build a lattice:

    >>> lattice1 = network1.lattice()
    >>> print(lattice1)
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
as an image file, assuming GraphViz etc is setup, use lattice1.plot(filename),
e.g lattice1.plot("n1_lattice.pdf") for a PDF file. Here is a simple text
graphic of this lattice diagram::

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
and edge type two (dotted) arrows run from node 4 to 1 and from 2 to 3::

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
    >>> print(network3)
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
    ...     print("%r or %s" % (p, cyclic_partition(p)))
    [0, 1, 0, 1, 1] or (13)(245)
    [0, 1, 0, 1, 2] or (13)(24)(5)
    [0, 1, 2, 3, 1] or (1)(25)(3)(4)
    [0, 1, 2, 3, 4] or (1)(2)(3)(4)(5)

Taking the second of these partitions gives a three cell quotient network:

    >>> print(network3.quotient([0, 1, 0, 1, 2]))
    (0,0) (0,1) (0,0) node 1+3
    (1,0) (0,0) (0,0) node 2+4
    (1,0) (0,0) (0,0) node 5

And using the balanced equivalence relations to build the lattice:

    >>> lattice = network3.lattice()
    >>> print(lattice)
    0 0 0 0 (13)(245)
    1 0 0 0 (13)(24)(5)
    1 0 0 0 (1)(25)(3)(4)
    0 1 1 0 (1)(2)(3)(4)(5)

The lattice is an undirected graph (here with four nodes, listed on the right
of the output in cyclic notation) which is represented on the left as a lower
triangular matrix. Graphically::

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

The idea is you can edit the examples in last section of this file to run this
program on particular networks of interest.

Building a reduced lattice (Example 6.1 in the 2021 manuscript):

    >>> network = CoupledCellNetwork([[0, 0, 0, 2],
    ...                               [0, 0, 0, 2],
    ...                               [0, 1, 0, 1],
    ...                               [0, 1, 0, 1]])
    >>> print(network.lattice())
    0 0 0 0 0 0 0 0 0 (1234)
    1 0 0 0 0 0 0 0 0 (124)(3)
    1 0 0 0 0 0 0 0 0 (12)(34)
    1 0 0 0 0 0 0 0 0 (13)(24)
    1 0 0 0 0 0 0 0 0 (1)(234)
    0 1 1 0 0 0 0 0 0 (12)(3)(4)
    0 1 0 1 1 0 0 0 0 (1)(24)(3)
    0 0 1 0 1 0 0 0 0 (1)(2)(34)
    0 0 0 0 0 1 1 1 0 (1)(2)(3)(4)
    >>> for reduced in network.reduced_lattices():
    ...     print("Reduction to %i nodes:" % reduced.n)
    ...     print(reduced)  # reduced lattice
    ...
    Reduction to 7 nodes:
    0 0 0 0 0 0 0 (1234) eta=1
    1 0 0 0 0 0 0 (124)(3)+(1)(234) eta=1
    1 0 0 0 0 0 0 (12)(34) eta=1
    1 0 0 0 0 0 0 (13)(24) eta=1
    0 1 1 0 0 0 0 (12)(3)(4)+(1)(2)(34) eta=0
    0 1 0 1 0 0 0 (1)(24)(3) eta=0
    0 0 0 0 1 1 0 (1)(2)(3)(4) eta=0

As an alternative, the example function "go" can be used to do this and
additional output and graph plotting, go(network, "file_name", reduce=True),
as shown in the user-editable section below.

The reduced lattice is an undirected graph, here with seven nodes, drawn by
hand::

    *------------------------------------------*
    |                                          |
    | Rank 1:                 (1234)           |
    |                        /  |   \          |
    |                       /   |    \         |
    |                      /    |     \        |
    | Rank 2:       (124)(3) (12)(34) (13)(24) |
    |               (1)(234)  /         /      |
    |                |      \/         /       |
    |                |      /\        /        |
    |                |     /  \      /         |
    | Rank 3:       (12)(3)(4) (1)(24)(3)      |
    |               (1)(2)(34)     /           |
    |                     \       /            |
    |                      \     /             |
    |                       \   /              |
    | Rank 4:            (1)(2)(3)(4)          |
    |                                          |
    *------------------------------------------*

In general there will be zero or more candidate reductions returned.
"""

import doctest
import os
import sys
import time
from itertools import product

if "-v" in sys.argv or "--version" in sys.argv:
    print("Coupled Cell Network graphs.py v1.0.5")
    sys.exit(0)

try:
    import numpy as np
except ImportError:
    sys.exit("Please install NumPy")

try:
    import pydot
except ImportError:
    try:
        # Try this API compatible version which works on Python 3
        import pydot_ng as pydot
    except ImportError:
        pydot = None
        sys.stderr.write(
            "Please install pyparsing, pydot & graphviz if you want "
            "to draw graphs or lattice diagrams\n"
        )


MAXLAT = 42294  # Increase this if you have a big powerful computer and are patient ;)


class NotBalancedError(ValueError):
    """Not a balanced equivalence relation."""

    pass


def bell_number(n):
    """Return bell number (number of ways to partition a set of n elements).

    The size of a lattice of balanced equivalence relations for n-cell
    regular network is the bell number of a set of n elements when the
    adjacency matrix of the regular network has only two distinct eigenvalues
    lambda_1 (valency) and lambda_2, where the geometric multiplicity of lambda_2
    is the same as the algebraic multiplicity of lambda_2.
    """
    bell = np.zeros((n + 1, n + 1), np.int)
    bell[0, 0] = 1

    for i in range(n):
        bell[i + 1, 0] = bell[i, i]
        for j in range(i + 1):
            bell[i + 1, j + 1] = bell[i, j] + bell[i + 1, j]

    return bell[n, 0]


def make_quotient(adj_matrix, partition):
    """Return quotient adjacency matrix, or raises an exception.

    This is the core algorithm described in the manuscript.

    The first argument (adj_matrix) should be a square matrix with
    non-negative integer elements representing an adjacency matrix (single
    edge type), with the second argument (partition) a partition of the
    nodes.

    For an n by n adjacency matrix (i.e. a graph with n nodes), the partition
    should be given as a list of n integers between 0 and k-1, where k is the
    number of equivalence classes. The function will return a k by k quotient
    matrix where the partition is a balanced equivalence relation (balanced
    colouring), or raise an exception if the partition is not balanced.

    This partition representation was chosen to simplify this function. See
    function cyclic_partition for displaying such a partition representations
    in a human friendly way, and function possible_partitions for generating
    them.

    For example, consider this 3 by 3 matrix:

    >>> import numpy as np
    >>> a = np.array([[1,0,1],[2,0,1],[0,0,2]])
    >>> print(a)
    [[1 0 1]
     [2 0 1]
     [0 0 2]]

    Take the trivial partition (no node merging):

    >>> print(make_quotient(a, [0, 1, 2]))
    [[1 0 1]
     [2 0 1]
     [0 0 2]]

    The following partition is a balanced equivalence relation merging the
    first and last nodes (nodes 0 and 2 in python zero based counting, 1 and 3
    in normal one based counting) giving a 2x2 quotient matrix:

    >>> print(make_quotient(a, [0, 1, 0]))
    [[2 0]
     [3 0]]

    Finally, here is a partition which is not a balanced equivalence relation,
    and therefore hasn't got a quotient matrix:

    >>> print(make_quotient(a, [0, 1, 1]))
    Traceback (most recent call last):
    ...
    NotBalancedError: Not a balanced equivalence relation
    """
    n = len(adj_matrix)  # gives number of rows, but should be square matrix
    if (n, n) != adj_matrix.shape:
        raise ValueError(
            "Adjacency matrix should be square, not %s" % repr(adj_matrix.shape)
        )
    assert len(partition) == n, "Bad partition"

    # Get number of merged nodes (equivalence classes):
    k = max(partition) + 1
    # Check partition entries are 0,1,2,...,k-1:
    assert set(partition) == set(range(k)), "Bad partition %r" % partition

    if k == n:
        # Trivial case of no merging
        return adj_matrix

    # Take column sums to go from n by n matrix to n by k matrix,
    q_matrix_step_1 = np.zeros([n, k], np.uint8)
    for (old_col, new_col) in enumerate(partition):
        q_matrix_step_1[:, new_col] += adj_matrix[:, old_col].astype(np.uint8)
    # Take representative rows to go from n by k matrix to k by k matrix,
    # note that we needlessly over-write n-k of the rows (simpler to just
    # do it rather than checking as it doesn't matter).
    q_matrix_step_2 = np.zeros([k, k], np.uint8)
    for (old_row, new_row) in enumerate(partition):
        q_matrix_step_2[new_row, :] = q_matrix_step_1[old_row, :]
    # Check that all the merged rows agreed (by comparing them to the
    # representative row picked above):
    for (old_row, new_row) in enumerate(partition):
        # NumPy's == gives element wise equality, giving k True/False values
        # We use the .all() to check if they are all True, thus all equal.
        if not (q_matrix_step_2[new_row, :] == q_matrix_step_1[old_row, :]).all():
            # Not a balanced colouring...
            raise NotBalancedError("Not a balanced equivalence relation")
    return q_matrix_step_2


def make_reduced_lattice_eigenvalue_matrix(adj_matrix, partition):
    """Make a reduced lattice eigenvalue matrix.

    Used as part of reducing a lattice, returns a new matrix or an exception.

    The operations are much like that for computing a quotient matrix
    or testing if an equivalence relation is balanced - see the function
    make_quotient defined above.

    First we combine columns, but not simply using addition - more of
    a binary operation where 0+0 -> 0; a+0 -> a; 0+a -> a; a+b-> error
    (where a and b are non-zero and distinct). Then we compare rows.
    """
    n = len(adj_matrix)  # gives number of rows, but should be square matrix

    assert (n, n) == adj_matrix.shape, "Matrix not square"
    assert len(partition) == n, "Bad partition"

    k = max(partition) + 1
    # Check partition entries are 0,1,2,...,k-1:
    assert set(partition) == set(range(k)), "Bad partition %r" % partition

    if k == n:
        # Trivial case of no merging
        return adj_matrix

    # Take column 'sums' to go from n by n matrix to n by k matrix,
    q_matrix_step_1 = np.zeros([n, k], np.uint8)
    for (old_col, new_col) in enumerate(partition):
        for row in range(n):
            if adj_matrix[row, old_col] == 0:
                pass
            elif q_matrix_step_1[row, new_col] == 0:
                q_matrix_step_1[row, new_col] = adj_matrix[row, old_col]
            elif q_matrix_step_1[row, new_col] != adj_matrix[row, old_col]:
                raise ValueError("Nope, can't merge row %r" % row)
    # Take representative rows to go from n by k matrix to k by k matrix,
    # note that we needlessly over-write n-k of the rows (simpler to just
    # do it rather than checking as it doesn't matter).
    q_matrix_step_2 = np.zeros([k, k], np.uint8)
    for (old_row, new_row) in enumerate(partition):
        q_matrix_step_2[new_row, :] = q_matrix_step_1[old_row, :]
    # Check that all the merged rows agreed (by comparing them to the
    # representative row picked above):
    for (old_row, new_row) in enumerate(partition):
        if not (q_matrix_step_2[new_row, :] == q_matrix_step_1[old_row, :]).all():
            raise ValueError("Not a valid lattice reduction")
    return q_matrix_step_2


def possible_partitions(n):
    """Return all possible partitions of n nodes as a Python list.

    Builds lists of n digits representing possible partitions of n nodes.
    For example:

    >>> for partition in possible_partitions(1):
    ...     print(partition)
    [0]
    >>> for partition in possible_partitions(2):
    ...     print(partition)
    [0, 0]
    [0, 1]
    >>> for partition in possible_partitions(3):
    ...     print(partition)
    [0, 0, 0]
    [0, 0, 1]
    [0, 1, 0]
    [0, 1, 1]
    [0, 1, 2]

    In this example we have five partitions: All three nodes merged ([0,0,0]),
    two cases where two nodes are merged, and the other trivial case of no
    nodes merged ([0,1,2]).

    Note that this outputs [0, 1, 1], but not [1, 0, 0] which is a different
    representation of the same node partition. The choice is arbitrary, but
    was made since Python considers [0, 1, 1] < [1, 0, 0].

    >>> for partition in possible_partitions(4):
    ...     print(partition)
    [0, 0, 0, 0]
    [0, 0, 0, 1]
    [0, 0, 1, 0]
    [0, 0, 1, 1]
    [0, 0, 1, 2]
    [0, 1, 0, 0]
    [0, 1, 0, 1]
    [0, 1, 0, 2]
    [0, 1, 1, 0]
    [0, 1, 1, 1]
    [0, 1, 1, 2]
    [0, 1, 2, 0]
    [0, 1, 2, 1]
    [0, 1, 2, 2]
    [0, 1, 2, 3]

    In this example we have five partitions: All three nodes merged ([0,0,0]),
    two cases where two nodes are merged, and the other trivial case of no
    nodes merged ([0,1,2]).

    Note that this outputs [0, 1, 1], but not [1, 0, 0] which is a different
    representation of the same node partition. The choice is arbitrary, but
    was made since Python considers [0, 1, 1] < [1, 0, 0].
    """
    # This is a recursive function
    if n < 1:
        raise ValueError("Require n at least one, not %r" % n)
    elif n == 1:
        yield [0]
    else:
        # Get the possible first n-1 digits of the possible
        # partitions by recursion:
        for p in possible_partitions(n - 1):
            # The possible final digit for these partitions is
            # given x = 0, 1, ..., max(p), max(p)+1
            # (which by construction means it will be at most n).
            for x in range(max(p) + 2):
                yield p + [x]


def possible_partitions_of_required_size(n, min_size):
    """Specialised version of function possible_partitions.

    >>> for partition in possible_partitions_of_required_size(2, 1):
    ...     print(partition)
    [0, 0]
    [0, 1]

    >>> for partition in possible_partitions_of_required_size(2, 2):
    ...     print(partition)
    [0, 1]

    >>> for partition in possible_partitions_of_required_size(3, 2):
    ...     print(partition)
    [0, 0, 1]
    [0, 1, 0]
    [0, 1, 1]
    [0, 1, 2]

    >>> for partition in possible_partitions_of_required_size(4, 3):
    ...     print(partition)
    [0, 0, 1, 2]
    [0, 1, 0, 2]
    [0, 1, 1, 2]
    [0, 1, 2, 0]
    [0, 1, 2, 1]
    [0, 1, 2, 2]
    [0, 1, 2, 3]

    >>> for partition in possible_partitions_of_required_size(4, 4):
    ...     print(partition)
    [0, 1, 2, 3]

    >>> len(list(possible_partitions_of_required_size(3, 2)))
    4
    >>> len(list(possible_partitions_of_required_size(4, 2)))
    14

    """
    # This is a recursive function
    if n < 1:
        raise ValueError("Require n at least one, not %r" % n)
    elif min_size < 1:
        raise ValueError("Require positive size, not %r" % min_size)
    elif min_size > n:
        raise StopIteration
        # raise ValueError("Size (%r) must be at most n (%r)" % (min_size, n))
    elif n == 1:
        if min_size == 1:
            yield [0]
    else:
        # Get the possible first n-1 digits of the possible
        # partitions by recursion:
        for p in possible_partitions_of_required_size(n - 1, max(1, min_size - 1)):
            # The possible final digit for these partitions is
            # given by x = 0, 1, ..., max(p), max(p)+1
            # (which by construction means it will be at most n).
            for x in range(max(p) + 2):
                if max(max(p), x) + 1 >= min_size:
                    yield p + [x]


def possible_partition_refinements(top):
    """Given a partition of n nodes, how can it be sub-partitioned.

    The main intended usage of this function is an optimisation in
    finding all possible balanced colourings. Given the top lattice
    node via the Aldis (2008) / Belykh and Hasler (2011) algorithm
    or otherwise, all the other latticed nodes (balanced colourings)
    will be partitions which are refinements of the top node patition.

    For example, consider 4 nodes and their partitions:

    >>> for partition in possible_partitions(4):
    ...     print("%s %s" % (partition, cyclic_partition(partition)))
    [0, 0, 0, 0] (1234)
    [0, 0, 0, 1] (123)(4)
    [0, 0, 1, 0] (124)(3)
    [0, 0, 1, 1] (12)(34)
    [0, 0, 1, 2] (12)(3)(4)
    [0, 1, 0, 0] (134)(2)
    [0, 1, 0, 1] (13)(24)
    [0, 1, 0, 2] (13)(2)(4)
    [0, 1, 1, 0] (14)(23)
    [0, 1, 1, 1] (1)(234)
    [0, 1, 1, 2] (1)(23)(4)
    [0, 1, 2, 0] (14)(2)(3)
    [0, 1, 2, 1] (1)(24)(3)
    [0, 1, 2, 2] (1)(2)(34)
    [0, 1, 2, 3] (1)(2)(3)(4)

    Suppose want to know the refinements of partition (12)(34) only?
    i.e. partition [0, 0, 1, 1] in our notation. You could make a
    short list of candidates based on the number of clusters, but
    that would include some false positives which are not relevant.

    >>> top = [0, 0, 1, 1]
    >>> print(cyclic_partition(top))
    (12)(34)
    >>> for partition in possible_partition_refinements(top):
    ...     print("%s %s" % (partition, cyclic_partition(partition)))
    [0, 0, 1, 1] (12)(34)
    [0, 0, 1, 2] (12)(3)(4)
    [0, 1, 2, 2] (1)(2)(34)
    [0, 1, 2, 3] (1)(2)(3)(4)

    This is especially important with larger networks. Here
    with five nodes:

    >>> top = [0, 0, 1, 1, 2]
    >>> print(cyclic_partition(top))
    (12)(34)(5)
    >>> for partition in possible_partition_refinements(top):
    ...     print("%s %s" % (partition, cyclic_partition(partition)))
    [0, 0, 1, 1, 2] (12)(34)(5)
    [0, 0, 1, 2, 3] (12)(3)(4)(5)
    [0, 1, 2, 2, 3] (1)(2)(34)(5)
    [0, 1, 2, 3, 4] (1)(2)(3)(4)(5)

    This should be identical to the brute-force solution:

    >>> for partition in possible_partitions(len(top)):
    ...     if partition_refinement(partition, top) or partition==top:
    ...         print("%s %s" % (partition, cyclic_partition(partition)))
    [0, 0, 1, 1, 2] (12)(34)(5)
    [0, 0, 1, 2, 3] (12)(3)(4)(5)
    [0, 1, 2, 2, 3] (1)(2)(34)(5)
    [0, 1, 2, 3, 4] (1)(2)(3)(4)(5)

    Note the original top node argument itself is included in the
    returned values.

    For a more extreme example,

    >>> len(list(possible_partitions(10)))
    115975
    >>> len(list(possible_partition_refinements([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])))
    115975
    >>> len(list(possible_partition_refinements([0, 0, 0, 1, 0, 1, 0, 0, 0, 1])))
    4385
    >>> len(list(possible_partition_refinements([0, 0, 0, 0, 0, 1, 0, 0, 0, 1])))
    8280

    """
    size = max(top)
    k = size + 1
    assert top[0] == 0
    assert set(range(k)) == set(top), "k = %i, yet %r" % (k, set(top))

    # for i in range(k):
    #    print(
    #        "Class %i, count %i, sub-partitions %r"
    #        % (i, top.count(i), list(possible_partitions(top.count(i))))
    # )
    pN = [possible_partitions(top.count(i)) for i in range(k)]
    for sub_parts in product(*pN):
        # print(sub_parts, "-->")
        new_partition = []
        old_index = [0] * k
        mapping = {}
        for old in top:
            sub = sub_parts[old][old_index[old]]
            # print(sub_parts, "-->", old, sub_parts[old], old_index[old], sub)
            if (old, sub) in mapping:
                new_partition.append(mapping[old, sub])
            else:
                if new_partition:
                    v = max(new_partition) + 1
                else:
                    v = 0
                new_partition.append(v)
                mapping[old, sub] = v
            old_index[old] += 1
        # print("-"*40, ">", new_partition)
        yield new_partition


def cyclic_partition(partition, sep="", captions=None):
    """Display a partition in cyclic form.

    If the matrix has n nodes, then the partition should be given as
    a list of n integers between 0 and k-1, where k is the number of
    merged nodes (number of equivalence classes). Such a partition
    representations is returned as a string using disjoint cyclic
    notation (by default counting the nodes from one rather than zero
    as is normal in Python):

    >>> print(cyclic_partition([0,0,1,0,2]))
    (124)(3)(5)

    You can insert a separator between nodes (useful if not all single
    digits):

    >>> print(cyclic_partition([0,0,1,0,2], "+"))
    (1+2+4)(3)(5)

    You can also supply node names:

    >>> print(cyclic_partition([0,0,1,0,2], "+", "abcde"))
    (a+b+d)(c)(e)
    >>> print(cyclic_partition([0,0,1,0,2], "+",
    ...                        ["red", "green", "blue", "black", "white"]))
    (red+green+black)(blue)(white)
    """
    # Get number of merged nodes (equivalence classes):
    if not captions:
        captions = [str(i + 1) for i in range(len(partition))]
    else:
        assert len(captions) == len(partition), "%i captions for %i nodes" % (
            captions,
            partition,
        )
    k = max(partition) + 1
    # assert set(partition) == set(range(k)), "Bad partition"
    # if not sep and k > 9: sep="+"
    mapping = {}
    for old, new in enumerate(partition):
        if new in mapping:
            mapping[new] += sep + captions[old]
        else:
            mapping[new] = captions[old]
    # Join up the captions for the k nodes (0, 1, ..., k-1)
    return "".join("(%s)" % mapping[new] for new in range(k))


def partition_refinement(partition_a, partition_b):
    """Check if a refines b, returns Boolean (True/False).

    Trivial example, for a homogeneous network the bottom nodes (nothing
    merged) refines the top node (all merged), but not the other way
    round:

    >>> partition_refinement([0,1,2,3],[0,0,0,0])
    True
    >>> partition_refinement([0,0,0,0],[0,1,2,3])
    False

    Less trivially,

    >>> partition_refinement([0,1,2,2,1],[0,0,1,1,0])
    True
    >>> partition_refinement([0,0,1,1,0],[0,1,2,2,1])
    False

    Note that a partition is not considered a refinement of itself:

    >>> partition_refinement([0,1,1,2],[0,1,1,2])
    False

    A few more examples, six balanced colourings of a particular graph:

    >>> partitions = [[0,0,0,0,0],
    ...               [0,0,0,0,1],
    ...               [0,1,0,1,1],
    ...               [0,1,0,1,2],
    ...               [0,1,2,3,1],
    ...               [0,1,2,3,4]]
    >>> for a in partitions:
    ...    cyc_a = cyclic_partition(a)
    ...    for b in partitions:
    ...        cyc_b = cyclic_partition(b)
    ...        if partition_refinement(a,b):
    ...           print("%s refines %s" % (cyc_a, cyc_b))
    (1234)(5) refines (12345)
    (13)(245) refines (12345)
    (13)(24)(5) refines (12345)
    (13)(24)(5) refines (1234)(5)
    (13)(24)(5) refines (13)(245)
    (1)(25)(3)(4) refines (12345)
    (1)(25)(3)(4) refines (13)(245)
    (1)(2)(3)(4)(5) refines (12345)
    (1)(2)(3)(4)(5) refines (1234)(5)
    (1)(2)(3)(4)(5) refines (13)(245)
    (1)(2)(3)(4)(5) refines (13)(24)(5)
    (1)(2)(3)(4)(5) refines (1)(25)(3)(4)
    """
    # if partition_a == partition_b:
    #    return False
    # assert len(partition_a) == len(partition_b)
    # rank_a = max(partition_a) #works but assumes partition format
    # rank_b = max(partition_b)
    rank_a = len(set(partition_a))
    rank_b = len(set(partition_b))
    if rank_a <= rank_b:
        return False

    for i in range(rank_a):
        # assert i in partition_a, "Bad partition? %r" % partition_a
        # See where number "i" occurs in partition a,
        positions = [p for p, v in enumerate(partition_a) if v == i]
        # Make sure these all belong to the same partition in b
        if len({partition_b[p] for p in positions}) > 1:
            # Failed - b is not a refinement (sub partition) of a
            return False
    return True


# A few extra tests,
assert not partition_refinement([0, 0, 1, 1, 0], [0, 1, 2, 3, 2])
assert not partition_refinement([0, 1, 2, 3, 2], [0, 0, 1, 1, 0])
assert not partition_refinement([0, 0, 1, 1, 0], [0, 1, 2, 2, 1])
assert partition_refinement([0, 1, 2, 2, 1], [0, 0, 1, 1, 0])

assert not partition_refinement([0, 1, 0, 1, 2], [0, 1, 0, 1, 2])
assert partition_refinement([0, 1, 0, 1, 2], [0, 1, 0, 1, 1])
assert partition_refinement([0, 1, 0, 1, 2], [0, 0, 0, 0, 0])
assert partition_refinement([0, 1, 0, 1, 1], [0, 0, 0, 0, 0])
assert partition_refinement([0, 1, 2, 3, 4], [0, 0, 0, 0, 0])

assert partition_refinement([0, 1, 2, 3, 4], [0, 1, 2, 3, 2])
assert partition_refinement([0, 1, 2, 3, 4], [0, 1, 0, 0, 2])
assert partition_refinement([0, 1, 2, 3, 2], [0, 0, 1, 0, 1])
assert partition_refinement([0, 0, 1, 0, 2], [0, 0, 1, 0, 1])
assert partition_refinement([0, 0, 1, 0, 1], [0, 0, 0, 0, 0])

assert partition_refinement([0, 1, 2, 3, 4], [0, 0, 0, 0, 0])
assert not partition_refinement([0, 0, 0, 0, 0], [0, 1, 2, 3, 4])
assert not partition_refinement([0, 1, 0, 0, 1], [0, 0, 0, 1, 1])

assert partition_refinement([0, 1, 1, 2], [0, 1, 1, 0])
assert partition_refinement([0, 1, 1, 2], [0, 1, 1, 1])
assert not partition_refinement([0, 1, 1, 0], [0, 1, 1, 2])
assert not partition_refinement([0, 1, 1, 1], [0, 1, 1, 2])


def make_partition(classifiers):
    """Make a partition from list/tuple/sequence of hashable elements.

    Constructs a partition using our representation as a list of integers,
    for example:

    >>> make_partition("apple")
    [0, 1, 1, 2, 3]
    >>> make_partition(["A", "B", "B", "A"])
    [0, 1, 1, 0]

    This can also be used to 'normalise' an integer classification:

    >>> make_partition([3, 2, 2, 1])
    [0, 1, 1, 2]

    Note a list of lists is not suitable, convert this into a list of tuples
    or another hashable datatype first.
    """
    partition = []
    mapping = {}
    for x in classifiers:
        if x in mapping:
            partition.append(mapping[x])
        elif mapping:
            partition.append(max(partition) + 1)
            mapping[x] = partition[-1]
        else:
            # First element
            assert not partition
            partition = [0]
            mapping[x] = 0
    return partition


def go(a, name="", format="png", top_only=False, reduce=False):
    """Take a graph, draw it, then print and draw the lattice.

    Given a CoupledCellNetwork object (and an optional name for it), it shows
    the adjacency matrix (and draws it to a file), those partitions which are
    balanced equivalence classes and their associated quotient networks, then
    the resulting lattice (which is also drawn to a file).

    The optional name argument is used to assign image filenames (requires
    GraphViz etc), together with the optional the format argument which
    specifies the image type (e.g. png or pdf, defaulting to png).

    If the optional argument top_only is set to True, rather than computing all
    the balanced colourings and the lattice, only the top lattice node is found
    using an algorithm combining Aldis (2008) and Belykh and Hasler (2011).

    If the optional argument reduce is set to True (not compatible with using
    the top_only setting), then any candidate lattice reductions are found.
    """
    if not isinstance(a, CoupledCellNetwork):
        raise TypeError("Expected a CoupledCellNetwork object.")
    print("=" * 50)
    print("")
    print("%s network adjacency matrix:" % name)
    print(a)
    if a.n > 9:
        sep = "+"
    else:
        sep = ""
    if name and not os.path.isfile("%s.%s" % (name, format)):
        a.plot("%s.%s" % (name, format))
    # Find top node (quick)
    start = time.time()
    p, q = a.top_lattice_node()
    taken = time.time() - start
    if name and not os.path.isfile("%s_top_node.%s" % (name, format)):
        q.plot("%s_top_node.%s" % (name, format))
    if top_only:
        print("")
        print(
            "Lattice top node "
            "(balanced equlivalence relationship with least clusters) %0.1fs:" % taken
        )
        print("")
        print("Partition: %s" % cyclic_partition(p, sep))
        print("")
        print("Quotient matrix:")
        print(q)
    else:
        # Calculate the quotients with partitions and the resulting lattice
        # (all in one go to avoid duplicated computation, you could also use
        # the a.quotients_with_partitions() method or similar if you didn't
        # want the lattice).
        start = time.time()
        lattice = a.lattice()
        taken = time.time() - start
        for p in lattice.partitions:
            print("")
            print("Partition %s, quotient matrix:" % cyclic_partition(p, sep))
            q = a.quotient(p)
            if max(p) + 1 < a.n:
                # Non-trivial
                print(q)
            else:
                print("(trivial)")
        # l = a.lattice()
        print("")
        print("%s Lattice matrix:" % name)
        print(lattice)
        print("")
        if taken < 60:
            print("Lattice took %0.1fs" % taken)
        elif taken < 360:
            print("Lattice took %0.1fmins" % (taken / 60))
        else:
            print("Lattice took %0.1fhours" % (taken / 360))
        if name and not os.path.isfile("%s_lattice.%s" % (name, format)):
            lattice.plot("%s_lattice.%s" % (name, format))
        print("(%i lattice nodes)" % lattice.n)
        if reduce:
            print("")
            print("Looking for lattice reductions...")
            start = time.time()
            candidates = 0
            for reduction in a.reduced_lattices():
                candidates += 1
                print("")
                print(
                    "Reduction candidate %i, with %i nodes:" % (candidates, reduction.n)
                )
                print(reduction)
                # Could draw it too, e.g.
                # reduction.plot("%s_reduced_lattice_%i.png" % (name, candidate))")
            taken = time.time() - start
            print("")
            if taken < 60:
                print("Lattice reduction took %0.1fs" % taken)
            elif taken < 360:
                print("Lattice reduction took %0.1fmins" % (taken / 60))
            else:
                print("Lattice reduction took %0.1fhours" % (taken / 360))


class CoupledCellLattice(object):
    """Balanced equivalence relation lattice."""

    def __init__(self, *partitions):
        r"""Initialise a CoupledCellLattice object.

        You can create a lattice object "by hand" if required. However, you are
        normally expected to create a CoupledCellNetwork object and use its
        lattice method to construct the lattice for you.

        For example, take the five cell network #1, which has
        four partitions which are balanced equivalence relations:

        [0,1,0,1,1] or (13)(245)
        [0,1,0,1,2] or (13)(24)(5)
        [0,1,2,3,1] or (1)(25)(3)(4)
        [0,1,2,3,4] or (1)(2)(3)(4)(5)

        To create a four node lattice from these four partitions, using the
        algorithm described in the manuscript "Computation of Balanced Equivalence
        Relations and their Lattice for a Coupled Cell Network", SIAM Journal
        on Applied Dynamical Systems (SIADS) 12(1), pp. 352-382:

        >>> lattice = CoupledCellLattice([0,1,0,1,1],
        ...                              [0,1,0,1,2],
        ...                              [0,1,2,3,1],
        ...                              [0,1,2,3,4])
        >>> print(lattice)
        0 0 0 0 (13)(245)
        1 0 0 0 (13)(24)(5)
        1 0 0 0 (1)(25)(3)(4)
        0 1 1 0 (1)(2)(3)(4)(5)

        The output represents a four by four lower triangular adjacency matrix
        with 1 and 0 entries (with and without a connection) on the left, and
        the node captions on the right (in cyclic notation). Notice this has
        just four edges, represented graphically::

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

        The node captions and ranks are determined from the partitions (balanced
        equivalence relations, also known as balanced colourings) automatically.

        You can access a number of key properties as attributes of the object:

        >>> lattice.n  # number of nodes
        4
        >>> lattice.matrix  # edge matrix as numpy integer array
        array([[0, 0, 0, 0],
               [1, 0, 0, 0],
               [1, 0, 0, 0],
               [0, 1, 1, 0]], dtype=uint8)
        >>> lattice.partitions  # partitions for each node
        ([0, 1, 0, 1, 1], [0, 1, 0, 1, 2], [0, 1, 2, 3, 1], [0, 1, 2, 3, 4])
        >>> lattice.ranks  # rank for each node
        [2, 3, 4, 5]
        >>> lattice.colors  # automatically assigned for use in plotting
        ['white', 'white', 'white', 'grey']
        >>> lattice.captions  # automatically assigned for printing & plotting
        ['(13)(245)', '(13)(24)(5)', '(1)(25)(3)(4)', '(1)(2)(3)(4)(5)']
        """
        self.n = n = len(partitions)  # number of lattice nodes
        if n > MAXLAT:
            raise ValueError(
                "Excessive lattice size %i nodes, MAXLAT = %i" % (n, MAXLAT)
            )

        self.partitions = partitions
        if max(max(p) for p in partitions) > 9:
            sep = "+"
        else:
            sep = ""
        self.captions = [cyclic_partition(p, sep) for p in partitions]

        # trivial_partitions = [[0] * len(partitions[0]), range(len(partitions[0]))]

        colors = ["white"] * n  # default white, grey for top and bottom nodes
        # This assumes the partitions are sorted:
        if partitions[0] == [0] * len(partitions[0]):
            colors[0] = "grey"
        if partitions[-1] == list(range(len(partitions[0]))):
            colors[-1] = "grey"
        self.colors = colors

        self.ranks = ranks = [len(set(p)) for p in partitions]

        try:
            refinement = np.zeros((n, n), np.uint8)
            # Note this includes zeros (False)  on the diagonal:
            # refinement = np.array([[partition_refinement(a,b) for b in partitions]\
            #                        for a in partitions], np.bool)
        except MemoryError:
            print("Out of memory problem? Lattice of %i partitions" % n)
            for p, r in zip(partitions, ranks):
                print(p, r)
            raise MemoryError("Out of memory problem? Lattice of %i partitions" % n)

        # Making a lower triangular matrix
        ep = list(enumerate(partitions))
        for row, a in ep:
            # This is equivalent to but faster than:
            # for col, b in ep:
            #    if col < row:
            #        refinement[row,col] = partition_refinement(a,b)
            rank_a = ranks[row]
            enu_a = list(enumerate(a))
            all_positions = [[p for p, v in enu_a if v == i] for i in range(rank_a)]
            for col, b in ep:
                # Want col < row as building a lower triangular matrix
                # The ranks must change too, a quick comparison
                if col < row and rank_a > ranks[col]:
                    # refinement[row,col] = partition_refinement(a,b)
                    linked = True
                    for positions in all_positions:
                        if len({b[p] for p in positions}) > 1:
                            # Failed - b is not a refinement (sub partition) of a
                            linked = False
                            break
                    if linked:
                        refinement[row, col] = linked

        # Now remove redundant edges
        try:
            edge_matrix = refinement.copy()
            # TODO - Can we do this in-situ (so needing less RAM)?
        except MemoryError:
            print(
                "Out of memory problem removing redundant edges "
                "in lattice of %i partitions" % n
            )
            raise MemoryError("Out of memory problem? Lattice of %i partitions" % n)
        for row in range(n):
            for col in range(n):
                if col < row and refinement[row, col] and ranks[col] + 1 < ranks[row]:
                    # This edge jumps at least one rank, so could be redundant
                    # print("Checking %s <-> %s" % (cyclic_partition(partitions[row]),
                    #                               cyclic_partition(partitions[col]))
                    for mid in range(col, row):
                        if refinement[row, mid] and refinement[mid, col]:
                            # This other point provides a route => edge was redundant
                            # print(
                            #     "Found %s <-> %s <-> %s"
                            #     % (cyclic_partition(partitions[row]),
                            #        cyclic_partition(partitions[mid]),
                            #        cyclic_partition(partitions[col]))
                            # )
                            edge_matrix[row, col] = False
                            break

        self.matrix = edge_matrix

    def __str__(self):
        """Return string representation of the matrix, used by the print command."""
        answer = []
        n = self.n
        x = max(len(str(self.matrix[i, j])) for i in range(n) for j in range(n))
        for i in range(self.n):
            answer.append(
                " ".join([str(self.matrix[i, j]).ljust(x) for j in range(n)])
                + " %s" % self.captions[i]
            )
        return "\n".join(answer)

    def __repr__(self):
        """Return string representation of the matrix."""
        return "%s(%r)" % (self.__class__.__name__, self.matrix.tolist())

    def plot(self, filename):
        """Use this function to produce an image file of the graph.

        It uses pydot library to talk to GraphViz to do this.
        """
        n = self.n
        if pydot is None:
            sys.stderr.write(
                "Please install graphviz, pydot and pyparsing "
                "to draw lattice diagrams\n"
            )
            return

        # Use the OS path module's "Split Extension" function
        extension = os.path.splitext(filename)[-1]
        # Check that there was a leading "."
        assert extension[0] == os.path.extsep
        # Remove the leading "." and make lower case (to give to graphviz as format):
        extension = extension[1:].lower()

        g = pydot.Dot(graph_type="graph")

        # Divide into subgraphs, one for each lattice rank.
        # This is a GraphViz layout trick so that each rank is
        # shown on its own horizontal level in the image.
        for rank in range(1, max(self.ranks) + 1):
            sub_g = pydot.Subgraph(
                graph_type="graph", graph_name="Rank%i" % rank, rank="same"
            )  # GraphViz rank is about layout
            sub_g.add_node(pydot.Node("Rank %i" % rank, shape="none"))
            for i in range(self.n):
                # Only add this node to the sub graph if it has the right rank.
                if self.ranks[i] == rank:
                    sub_g.add_node(
                        pydot.Node(
                            self.captions[i], fillcolor=self.colors[i], style="filled"
                        )
                    )
            g.add_subgraph(sub_g)
            if rank != 1:
                g.add_edge(pydot.Edge("Rank %i" % (rank - 1), "Rank %i" % (rank)))

        # Now that we have added all the nodes, we can do the edges:
        for EdgeTo in range(n):
            for EdgeFrom in range(n):
                if self.matrix[EdgeTo, EdgeFrom]:
                    g.add_edge(
                        pydot.Edge(
                            self.captions[EdgeFrom],
                            self.captions[EdgeTo],
                            style="solid",
                        )
                    )

        # Now get the graph object to write itself to an image file.
        try:
            g.write(filename, prog="dot", format=extension)
        except pydot.InvocationException:
            sys.stderr.write("Please check graphviz is installed\n")


class CoupledCellNetwork(object):
    """Class to represent a Couple Cell Network (CCN) with multiple edge types.

    Multiple edge types are supported, specified and represented internally
    as a separate adjacency matrix for each edge type.
    """

    def __init__(self, *edge_matrices):
        """Initialise a CoupledCellNetwork object.

        edge_matrices - One or more adjacency matrices to define the
                        connections. Use one matrix for each edge type.
                        Each matrix should either be a square numpy array
                        object, square numpy matrix, or a list of lists of
                        integers. The matrix elements should be non-negative
                        integers.

        e.g.

        >>> ccn1 = CoupledCellNetwork([[1,0,1],[1,1,0],[0,0,1]])
        >>> print(ccn1)
        1 0 1 node 1
        1 1 0 node 2
        0 0 1 node 3

        >>> ccn2 = CoupledCellNetwork([[1,0,1],[1,1,0],[0,0,1]],
        ...                           [[0,2,1],[0,0,0],[2,0,0]])
        >>> print(ccn2)
        (1,0) (0,2) (1,1) node 1
        (1,0) (1,0) (0,0) node 2
        (0,2) (0,0) (1,0) node 3

        You can access two key properties as attributes of the object:

        >>> ccn2.n  # number of nodes
        3
        >>> len(ccn2.matrices)  # one numpy integer array per edge type
        2
        """
        assert isinstance(edge_matrices, tuple), edge_matrices
        assert edge_matrices
        n = None
        self.matrices = []
        for matrix in edge_matrices:
            if n is None:
                n = len(matrix)
                self.n = n
            # Turn it into a numpy array (in case it wasn't already)
            matrix = np.array(matrix, np.uint8)
            if (n, n) != matrix.shape:
                raise ValueError(
                    "All edge matrixes must be (%i, %i), not %r" % (n, n, matrix.shape)
                )
            self.matrices.append(matrix)

    def __str__(self):
        """Return string representation of the matrix, used by the print command."""
        answer = []
        n = self.n
        try:
            captions = self.captions
        except AttributeError:
            captions = [str(i + 1) for i in range(n)]
        if len(self.matrices) == 1:
            m = self.matrices[0]
            x = max(len(str(m[i, j])) for i in range(n) for j in range(n))
            for i in range(self.n):
                answer.append(
                    " ".join([str(m[i, j]).ljust(x) for j in range(n)])
                    + " node "
                    + captions[i]
                )
        else:
            x = max(
                len(str(m[i, j]))
                for i in range(n)
                for j in range(n)
                for m in self.matrices
            )
            for i in range(self.n):
                answer.append(
                    " ".join(
                        [
                            "(%s)"
                            % ",".join([str(m[i, j]).ljust(x) for m in self.matrices])
                            for j in range(n)
                        ]
                    )
                    + " node "
                    + captions[i]
                )
        return "\n".join(answer)

    def __repr__(self):
        """Return string representation of the matrix."""
        answer = ", ".join(repr(m.tolist()) for m in self.matrices)
        return "%s(%s)" % (self.__class__.__name__, answer)

    def input_driven_refinement(self, partition):
        """Input driven partition refinement based on Belykh and Hasler (2011)'s method.

        Part of the Belykh and Hasler (2011) algorithm computes what they called the
        input driven refinement of a partition.  Their description only considers
        single edge types. This implementation generalises their method to consider
        multiple edge types as well, giving a simplified version of Aldis (2008).

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

        The trivial partition where all 5 nodes are in the same cluster is
        represented as the list [0, 0, 0, 0, 0], thus:

        >>> network3.input_driven_refinement([0, 0, 0, 0, 0])
        [0, 1, 0, 1, 1]
        >>> network3.input_driven_refinement([0, 1, 0, 1, 1])
        [0, 1, 0, 1, 1]
        >>> print(cyclic_partition([0, 1, 0, 1, 1]))
        (13)(245)

        See also the method top_lattice_node which calls this method repeatedly
        starting with the trivial partition into one cluster, until it reaches
        the top lattice node (aka balanced colouring with minimal number of colours).

        Here is a larger example using repeated edges (but only one edge type):

        >>> big = CoupledCellNetwork([
        ...     [0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0],
        ...     [0, 2, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 2, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 1],
        ...     [0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 2, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 1, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 2, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 2, 0]])
        >>> print(big.input_driven_refinement([0]*16))
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        """
        colors = max(partition) + 1

        # Basic testing of partition validity
        assert len(partition) == self.n
        assert 0 == min(partition)
        assert max(partition) < self.n
        assert set(partition) == set(range(colors))

        assert isinstance(self.matrices, list), self.matrices

        # This could be written far more elegantly, but it isn't a computational
        # bottleneck so I have not attempted to polish or speed it up at all.

        def input_color_count(node, matrix, n, color, partition):
            # For each other potential input node i,
            # is it connected and is it the requested colour?
            # Note the NUMBER of input edged from each node is important
            return sum(matrix[node, i] for i in range(n) if partition[i] == color)

        def single_edge_color_counts(node, matrix, n, colors, partition):
            assert node < n
            assert matrix.shape == (n, n), "%r vs %r" % (matrix.shape, n)
            return tuple(
                input_color_count(node, matrix, n, c, partition) for c in range(colors)
            )

        def color_counts(node, matrix_list, n, colors, partition):
            assert node < n
            assert isinstance(matrix_list, list), matrix_list
            assert matrix_list[0].shape == (n, n)
            return tuple(
                single_edge_color_counts(node, m, n, colors, partition)
                for m in matrix_list
            )

        # Make list of tuples, one entry for each node containing:
        # -current color (partition class, integer), c
        # -input color counts for edge type 1
        # -input color counts for edge type 2
        # -...
        # -input color counts for edge type M
        answer = []
        assigned_color = {}
        for i, c in enumerate(partition):
            assert isinstance(self.matrices, list), self.matrices
            input_colors = color_counts(i, self.matrices, self.n, colors, partition)
            key = (c, input_colors)
            if key in assigned_color:
                answer.append(assigned_color[key])
            else:
                if answer:
                    new_c = max(answer) + 1
                else:
                    new_c = 0
                answer.append(new_c)
                assigned_color[key] = new_c
        return answer

    def top_lattice_node(self):
        """Return the top/minimal partition and its quotient matrix.

        Uses an algorithm combining Aldis (2008) and Belykh and Hasler (2011)
        to find the lattice top node, which is the balanced equivalence relation
        using the fewest partition classes (aka balanced coloring with minimal
        number of colors). It starts with the trivial partition where all nodes
        are in the same class (all nodes given the same colour), and then
        calculates the input driven refinement. If this returns the same
        partition, then it is balanced.

        This case be regarded as a generalisation of Belykh and Hasler (2011)
        where we consider multiple edge types, and could the inputs to each
        node from each edge type and node class combination.

        Alternatively, this can equally be regarded as a simplification of the
        Aldis (2008) algorithm where the edge types are not tracked explicitly,
        skipping 'phase one' of his method. Instead all possible edge classes
        are considered (each arrow type and tail node class), including those
        not present which explicit tracking would skip.

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
        >>> p, q = network3.top_lattice_node()
        >>> print(p)
        [0, 1, 0, 1, 1]

        This more complicated example should match Figure 5 in Belykh and Hasler (2011),

        >>> figure5 = CoupledCellNetwork([[0,1,1,0,0,0,0,0,0,0],
        ...                               [1,0,0,1,1,1,0,0,0,0],
        ...                               [1,0,0,1,1,1,0,0,0,0],
        ...                               [0,1,1,0,0,0,1,1,1,1],
        ...                               [0,1,1,0,0,0,1,1,1,1],
        ...                               [0,1,1,0,0,0,1,1,1,1],
        ...                               [0,0,0,1,1,1,0,0,0,0],
        ...                               [0,0,0,1,1,1,0,0,0,0],
        ...                               [0,0,0,1,1,1,0,0,0,0],
        ...                               [0,0,0,1,1,1,0,0,0,0]])
        >>> figure5.input_driven_refinement([0]*10)
        [0, 1, 1, 2, 2, 2, 3, 3, 3, 3]

        In this case we only needed one call to reach a balanced equivalence
        relationship, thus:

        >>> p, q = figure5.top_lattice_node()
        >>> p
        [0, 1, 1, 2, 2, 2, 3, 3, 3, 3]

        Now a 16 cell example, which is interesting for using repeated edges (although
        there is only one edge type here):

        >>> big = CoupledCellNetwork([
        ...     [0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0],
        ...     [0, 2, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 2, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 1],
        ...     [0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 2, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  2, 0, 0, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 1, 0,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0],
        ...     [0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  1, 0, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 2, 0, 0],
        ...     [0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 2, 0],
        ... ])
        >>> print(big.input_driven_refinement([0]*16))
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        >>> p, q = big.top_lattice_node()
        >>> print(cyclic_partition(p, "+"))
        (1+2+3+4+5+6+7+8+9+10+11+12+13+14+15+16)

        """
        p = [0] * self.n
        while True:
            # Note that rather than calculating each input driven refinement
            # in isolation like this, many of the values could be recycled at
            # each iteration, reducing the computation (at the cost of more
            # complexity). Currently this is not a performance bottleneck.
            new = self.input_driven_refinement(p)
            if new == p:
                break
            p = new
        try:
            q = self.quotient(p)
            return p, q
        except NotBalancedError:
            raise RuntimeError("Top lattice node %r was not balanced" % p)

    def quotient(self, partition):
        """Return another quotient network adjacency matrix, or raise exception.

        If the matrix has n nodes, then the partition should be given as
        a list of n integers between 0 and k-1, where k is the number of
        merged nodes (equivalence classes).

        Consider this graph with two edge types:

        >>> e1 = [[0,1,0,1,0],
        ...       [1,0,0,1,0],
        ...       [0,0,1,0,1],
        ...       [1,1,0,0,0],
        ...       [0,0,1,0,1]]
        >>> e2 = [[0,1,0,0,0],
        ...       [0,0,0,1,0],
        ...       [1,0,0,0,0],
        ...       [1,0,0,0,0],
        ...       [1,0,0,0,0]]
        >>> a = CoupledCellNetwork(e1,e2)
        >>> print(a)
        (0,0) (1,1) (0,0) (1,0) (0,0) node 1
        (1,0) (0,0) (0,0) (1,1) (0,0) node 2
        (0,1) (0,0) (1,0) (0,0) (1,0) node 3
        (1,1) (1,0) (0,0) (0,0) (0,0) node 4
        (0,1) (0,0) (1,0) (0,0) (1,0) node 5

        For this graph there are many valid balanced equivalence relations
        (balanced colourings), including the trivial partition which merges
        all the nodes:

        >>> print(a.quotient([0,0,0,0,0]))
        (2,1) node 1+2+3+4+5

        This more interesting partition gives a 3 by 3 quotient network:

        >>> print(a.quotient([0,0,1,0,2]))
        (2,1) (0,0) (0,0) node 1+2+4
        (0,1) (1,0) (1,0) node 3
        (0,1) (1,0) (1,0) node 5

        """
        # If the partition is not a balanced colouring, then at least one
        # of the edge type specific matrices will raise a ValueError:
        q = CoupledCellNetwork(*[make_quotient(m, partition) for m in self.matrices])
        q.captions = [
            "+".join(str(j + 1) for j, p in enumerate(partition) if p == i)
            for i in range(max(partition) + 1)
        ]
        return q

    def quotients_with_partitions(self, resume_file=None):
        """Return quotients (smaller matrices) and partitions (generator function).

        This allows you to iterate over all valid quotient networks:

        >>> a = CoupledCellNetwork([[0,1,0,1],
        ...                         [0,0,1,0],
        ...                         [1,1,0,0],
        ...                         [0,0,1,0]])
        >>> for q, p in a.quotients_with_partitions():
        ...     print("Quotient with %i nodes: %s" % (q.n, cyclic_partition(p)))
        ...     print(q)
        Quotient with 3 nodes: (1)(24)(3)
        0 2 0 node 1
        0 0 1 node 2+4
        1 1 0 node 3
        Quotient with 4 nodes: (1)(2)(3)(4)
        0 1 0 1 node 1
        0 0 1 0 node 2
        1 1 0 0 node 3
        0 0 1 0 node 4

        In this example (with a single edge type) there are just two.
        """
        top = self.top_lattice_node()[0]
        if top == [0] * self.n:
            # Trival top node, this is faster:
            partitions = possible_partitions(self.n)
        else:
            partitions = possible_partition_refinements(top)
        if not resume_file:
            # Handle trivial top node as special case?
            for partition in partitions:
                try:
                    yield self.quotient(partition), partition
                except ValueError:
                    pass
        else:
            # Note - we now exploit the Aldis (2008) / Belylh & Hasler (2011)
            # based shortcut
            done = False
            count = 0
            last_good_partition = None
            if os.path.isfile(resume_file):
                # print("Resuming from %s" % resume_file)
                if resume_file.endswith(".gz"):
                    import gzip

                    h = gzip.open(resume_file)
                else:
                    h = open(resume_file, "rU")
                a = h.readline()
                assert a.strip() == repr(self)
                for line in h:
                    if line[0] == "[" and line[-2:] == "]\n":
                        last_good_partition = [
                            int(i.strip()) for i in line[1:-2].split(",")
                        ]
                        assert len(last_good_partition) == self.n
                        yield self.quotient(last_good_partition), last_good_partition
                        count += 1
                    elif line == "DONE\n":
                        done = True
                        break
                    elif line.startswith("DONE\t"):
                        assert line == "DONE\t%i\n" % count
                        done = True
                        break
                h.close()
                if not done:
                    h = open(resume_file, "a")
            else:
                # print("Checkpointing in %s" % resume_file)
                h = open(resume_file, "w")
                h.write("%r\n" % self)
            if not done:
                # TODO - Jump the known partitions from the file?
                # i.e. avoid rechecking they are balanced
                for partition in partitions:
                    if not last_good_partition:
                        try:
                            yield self.quotient(partition), partition
                            count += 1
                            h.write("%r\n" % partition)
                            h.flush()
                        except ValueError:
                            pass
                    elif partition == last_good_partition:
                        last_good_partition = None
                h.write("DONE\t%i\n" % count)
                h.close()
            # All done

    def quotients(self, resume_file=None):
        """Return quotients as smaller matrices (generator function).

        This allows you to iterate over all valid quotient networks:

        >>> a = CoupledCellNetwork([[0,1,0,1],
        ...                         [0,0,1,0],
        ...                         [1,1,0,0],
        ...                         [0,0,1,0]])
        >>> for q in a.quotients():
        ...     print("Quotient with %i nodes:" % q.n)
        ...     print(q)
        Quotient with 3 nodes:
        0 2 0 node 1
        0 0 1 node 2+4
        1 1 0 node 3
        Quotient with 4 nodes:
        0 1 0 1 node 1
        0 0 1 0 node 2
        1 1 0 0 node 3
        0 0 1 0 node 4

        In this example (with a single edge type) there are just two.
        """
        return (q for q, p in self.quotients_with_partitions(resume_file))

    def partitions(self, resume_file=None):
        """Return those partitions which are balanced (generator function).

        This allows you to iterate over all valid balanced equivalence relations:

        >>> a = CoupledCellNetwork([[0,1,0,1],
        ...                         [0,0,1,0],
        ...                         [1,1,0,0],
        ...                         [0,0,1,0]])
        >>> for partition in a.partitions():
        ...     print("%s %s" % (partition, cyclic_partition(partition)))
        [0, 1, 2, 1] (1)(24)(3)
        [0, 1, 2, 3] (1)(2)(3)(4)

        In this example (with a single edge type) there are just two.
        """
        # Note considerable code duplication with quotients_with_partitions
        # BUT this avoids needlessly rechecking previously identified
        # partitions are balanced.
        if resume_file:
            done = False
            count = 0
            last_good_partition = None
            if os.path.isfile(resume_file):
                # print("Resuming from %s" % resume_file)
                if resume_file.endswith(".gz"):
                    import gzip

                    h = gzip.open(resume_file)
                else:
                    h = open(resume_file, "rU")
                a = h.readline()
                assert a.strip() == repr(self)
                for line in h:
                    if line[0] == "[" and line[-2:] == "]\n":
                        last_good_partition = [
                            int(i.strip()) for i in line[1:-2].split(",")
                        ]
                        assert len(last_good_partition) == self.n
                        yield last_good_partition
                        count += 1
                    elif line == "DONE\n":
                        done = True
                        break
                    elif line.startswith("DONE\t"):
                        assert line == "DONE\t%i\n" % count
                        done = True
                        break
                h.close()
                if not done:
                    h = open(resume_file, "a")
            else:
                # print("Checkpointing in %s" % resume_file)
                h = open(resume_file, "w")
                h.write("%r\n" % self)
            if not done:
                # TODO - Jump the known partitions from the file?
                top = self.top_lattice_node()[0]
                if top == [0] * self.n:
                    # Trival top node, this is faster:
                    partitions = possible_partitions(self.n)
                else:
                    partitions = possible_partition_refinements(top)
                for partition in partitions:
                    if not last_good_partition:
                        # Is it balanced?
                        try:
                            q = self.quotient(partition)
                        except ValueError:
                            q = None
                        if q:
                            h.write("%r\n" % partition)
                            h.flush()
                            yield partition
                            count += 1
                    elif partition == last_good_partition:
                        last_good_partition = None
                h.write("DONE\t%i\n" % count)
                h.close()
            # All done
        else:
            for _q, p in self.quotients_with_partitions():
                yield p

    def lattice(self, caption_sep="+", resume_file=None):
        r"""Find balanced equivalence relations and build lattice.

        This method implements an alogrithm described in Kamei & Cock (2013).
        Consider graph #5 from that paper, which has two edge types:

        >>> e1 = [[0,1,0,1,0],
        ...       [1,0,0,1,0],
        ...       [0,0,1,0,1],
        ...       [1,1,0,0,0],
        ...       [0,0,1,0,1]]
        >>> e2 = [[0,1,0,0,0],
        ...       [0,0,0,1,0],
        ...       [1,0,0,0,0],
        ...       [1,0,0,0,0],
        ...       [1,0,0,0,0]]
        >>> a = CoupledCellNetwork(e1,e2)
        >>> print(a)
        (0,0) (1,1) (0,0) (1,0) (0,0) node 1
        (1,0) (0,0) (0,0) (1,1) (0,0) node 2
        (0,1) (0,0) (1,0) (0,0) (1,0) node 3
        (1,1) (1,0) (0,0) (0,0) (0,0) node 4
        (0,1) (0,0) (1,0) (0,0) (1,0) node 5
        >>> q = a.lattice()
        >>> print(q)
        0 0 0 0 0 (12345)
        1 0 0 0 0 (124)(35)
        0 1 0 0 0 (124)(3)(5)
        0 1 0 0 0 (1)(2)(35)(4)
        0 0 1 1 0 (1)(2)(3)(4)(5)

        This 5 node graph has a 5 node lattice which includes both trivial
        node (12345) from merging all graph nodes (rank 1), and also node
        (1)(2)(3)(4)(5) for merging no graph nodes (rank 5). Graphically::

            *--------------------------------------*
            |                                      |
            | Rank 1:         (12345)              |
            |                    |                 |
            |                    |                 |
            | Rank 2:        (124)(35)             |
            |                 /    \               |
            |                /      \              |
            | Rank 3:  (124)(3)(5)   \             |
            |                \        \            |
            |                 \        \           |
            | Rank 4:          \   (1)(2)(35)(4)   |
            |                   \      /           |
            |                    \    /            |
            | Rank 5:       (1)(2)(3)(4)(5)        |
            |                                      |
            *--------------------------------------*

        The entries in the adjacency matrix for the lattice show which nodes
        (partitions) cover other nodes (partitions): The first node (12345)
        covers the second node, (124)(35), which in turn covers the next two
        nodes, (124)(3)(5) and (1)(2)(35)(4), and they cover the final node
        (1)(2)(3)(4)(5), which covers no nodes.
        """
        partitions = list(self.partitions(resume_file))
        # Sort by rank (increasing), then by partition lexically.
        # rank(p) = max(p) + 1, thus equivalent to sort on max(p)
        partitions.sort(key=lambda p: (max(p), p))
        return CoupledCellLattice(*partitions)

    def reduced_lattices(self, caption_sep="+", resume_file=None):
        """Return zero or more reduced lattices.

        This method implements an alogrithm described in Kamei & Ruan (2021),
        and returns zero or more candidate lattice reductions. i.e. Smaller
        lattices created by merging nodes in the lattice constructed as per
        the alogrithm described in Kamei & Cock (2013).

        In an ideal case there will be a single unique reduction. In this
        example a seven node lattice has a unique six node reduction:

        >>> network = CoupledCellNetwork([[0, 0, 0, 2],
        ...                               [0, 0, 0, 2],
        ...                               [0, 1, 1, 0],
        ...                               [0, 0, 0, 2]])
        >>> lattice = network.lattice()
        >>> print(lattice)
        0 0 0 0 0 0 0 (1234)
        1 0 0 0 0 0 0 (124)(3)
        1 0 0 0 0 0 0 (1)(234)
        0 1 0 0 0 0 0 (12)(3)(4)
        0 1 0 0 0 0 0 (14)(2)(3)
        0 1 1 0 0 0 0 (1)(24)(3)
        0 0 0 1 1 1 0 (1)(2)(3)(4)
        >>> for reduced in network.reduced_lattices():
        ...     print("Reduction to %i nodes:" % reduced.n)
        ...     print(reduced)
        Reduction to 6 nodes:
        0 0 0 0 0 0 (1234) eta=1
        1 0 0 0 0 0 (124)(3) eta=1
        1 0 0 0 0 0 (1)(234) eta=1
        0 1 0 0 0 0 (12)(3)(4)+(14)(2)(3) eta=1
        0 1 1 0 0 0 (1)(24)(3) eta=0
        0 0 0 1 1 0 (1)(2)(3)(4) eta=0

        In general there will be multiple candidate reductions returned.
        Sometimes a network has a lattice with no reductions:

        >>> network = CoupledCellNetwork([[0, 1, 0, 1, 0],
        ...                               [1, 0, 0, 0, 1],
        ...                               [1, 0, 0, 1, 0],
        ...                               [1, 0, 1, 0, 0],
        ...                               [1, 1, 0, 0, 0]])
        >>> lattice = network.lattice()
        >>> print(lattice)
        0 0 0 0 0 0 0 0 0 0 (12345)
        1 0 0 0 0 0 0 0 0 0 (123)(45)
        1 0 0 0 0 0 0 0 0 0 (145)(23)
        1 0 0 0 0 0 0 0 0 0 (1)(2345)
        0 1 1 1 0 0 0 0 0 0 (1)(23)(45)
        0 0 0 1 0 0 0 0 0 0 (1)(24)(35)
        0 0 0 1 0 0 0 0 0 0 (1)(25)(34)
        0 0 0 0 0 0 1 0 0 0 (1)(2)(34)(5)
        0 0 0 0 0 0 1 0 0 0 (1)(25)(3)(4)
        0 0 0 0 1 1 0 1 1 0 (1)(2)(3)(4)(5)
        >>> for reduced in network.reduced_lattices():
        ...     print("Reduction to %i nodes:" % reduced.n)
        ...     print(reduced)

        i.e. No output from calling reduced_lattices here.
        """
        if len(self.matrices) > 1:
            raise NotImplementedError("Only one edge type implemented so far!")
        if len({sum(row) for row in self.matrices[0]}) != 1:
            # In regular network of valency r, all nodes have r inputs.
            raise NotImplementedError("Only regular networks can be reduced so far!")

        # If the lattice (self) is the partition lattice of n elements,
        # then the reduced lattice is given by L(1,...,1).
        # If the network size is more than 4, we return L(1,...,1) as an output without
        # searching candidates of reduced lattice.

        # print("Network size is %i" % self.n)

        lattice = self.lattice()
        # print("Lattice size is %i" % lattice.n)

        b_number = bell_number(self.n)
        # print("Bell number is %i" % b_number)

        if self.n > 4 and lattice.n == b_number:
            print("Reduced lattice is given by L(1,...,1)")
            # TODO - Can we build the appropriate lattice matrix?
            return

        q_and_p = list(self.quotients_with_partitions(resume_file))
        # Sort by rank (increasing), then by partition lexically.
        # rank(p) = max(p) + 1, thus equivalent to sort on max(p)
        q_and_p.sort(key=lambda q_p: (max(q_p[1]), q_p[1]))
        # n = len(q_and_p)

        PLACES = 6  # Used for fuzzy matching of eigenvalues
        all_eigen = sorted(
            round(e, PLACES) for e in np.linalg.eigvals(self.matrices[0])
        )
        unique_eigen = sorted(set(all_eigen), key=lambda c: (c.real, c.imag))

        partitions = [q_p[1] for q_p in q_and_p]
        assert partitions == sorted(partitions, key=lambda p: (max(p), p))
        quotients = [q_p[0] for q_p in q_and_p]
        for q in quotients:
            assert len(q.matrices) == 1, "Quotient has more edge types than parent!"
        assert partitions[-1] == list(range(self.n)), (
            "Bottom lattice node is %r" % partitions[-1]
        )

        eigenvalues = []
        for q in quotients:
            evals = sorted(round(e, PLACES) for e in np.linalg.eigvals(q.matrices[0]))
            for e in evals:
                assert e in all_eigen
                assert evals.count(e) <= all_eigen.count(e)
            eigenvalues.append(tuple(evals))  # Tuple is hashable, see make_partition

        # for p, q, e in zip(partitions, quotients, eigenvalues):
        #    print("Partition", p, "eigenvalues", e)
        #    print(q.matrices[0])
        #    print("")

        fat_lattice = CoupledCellLattice(*partitions)
        # Expand the captions with more details for exploration/debug use:
        fat_lattice.captions = [
            "#%i\n%s\n%r\n%r" % (i, caption, partitions[i], eigenvalues[i])
            for (i, caption) in enumerate(fat_lattice.captions)
        ]
        # fat_lattice.plot("reduction_debug.png")
        sym_fat = fat_lattice.matrix + fat_lattice.matrix.T
        # print(sym_fat)
        # print("")
        eigen_matrices = [
            sym_fat * np.array([evals.count(e) for evals in eigenvalues])
            for e in unique_eigen
        ]
        # print("")
        # for e, tmp in zip(unique_eigen, eigen_matrices):
        #    print("E for eigenvalue", e)
        #    print(tmp)
        #    print("")

        # Now do column combination (using a sort of 'or' operation)
        # and row comparison much like that used to identify if a
        # partition is balanced.
        # Only need to consider partitions which merge within the ranks
        # i.e. sub-partitions of the rank induced partition,
        # In fact, can go further in reducing the scope, and only
        # need try merging nodes with same set of eigenvalues.

        # print("Possible reductions:")
        for p in possible_partition_refinements(make_partition(eigenvalues)):
            # print("Possible reduction", p)
            ok = True
            for m in eigen_matrices:
                # for e, m in zip(sorted(set(all_eigen)), eigen_matrices):
                # print(p, e)
                try:
                    q = make_reduced_lattice_eigenvalue_matrix(m, p)
                except ValueError:
                    ok = False
                    break
            if not ok:
                continue
            # Prepare reduced lattice captions
            new_n = max(p) + 1
            q = np.zeros((new_n, new_n), np.uint8)
            assert q.shape == (new_n, new_n)
            reduced_nodes = [[] for i in range(new_n)]
            reduced_eigens = [None] * new_n
            reduced_ranks = [None] * new_n
            for i, node_class in enumerate(p):
                reduced_nodes[node_class].append(partitions[i])
                reduced_eigens[node_class] = eigenvalues[i]
                reduced_ranks[node_class] = max(partitions[i])  # !!!
                # print("reducing", i, node_class, partitions[i])
                # Simple merge of lattice connectivity...
                for j, tmp in enumerate(p):
                    if fat_lattice.matrix[i, j]:
                        q[node_class, tmp] = 1
            # What are the ((great)grand)parents of each lattice node?
            # Considering multiple path connections up the lattice.
            parents = np.identity(new_n, np.uint8)
            for _ in range(new_n):
                parents += np.dot(q, parents)
            # Assign eta
            eta = [0] * new_n
            for i in range(new_n):
                # start with eta = node's rank (i.e. partition size)
                r = reduced_ranks[i]
                e = r + 1
                # print("Reduced lattice node %s, start eta %r" % (i, e))
                for j in range(i):
                    assert reduced_ranks[j] <= r
                    if reduced_ranks[j] < r and parents[i, j]:
                        # print("Considering parent %s, with eta %r" % (j, eta[j]))
                        e -= eta[j]
                        # print("Reduced lattice node %s, revised eta %r" % (i, e))
                eta[i] = e
            if min(eta) < 0:
                continue
            # print("")
            # for i in range(new_n):
            #    print("+".join(cyclic_partition(tmp) for tmp in reduced_nodes[i]))
            #    print("evals=", reduced_eigens[i], "eta=", eta[i])
            # print(eta, "<-- eta")
            if min(eta) < 0:
                # print(p, "<-- Nice except eta negative")
                continue
            assert sum(eta) == self.n, "eta=%r, sum=%i, but n=%i" % (
                eta,
                sum(eta),
                self.n,
            )
            # print(p, "<-- Good, %i nodes post reduction" % (max(p) + 1))

            # Pick first of each merged node set as representative...
            reduced_lattice = CoupledCellLattice(
                *[reduced_nodes[i][0] for i in range(new_n)]
            )
            assert reduced_lattice.n == new_n, reduced_lattice
            reduced_lattice.captions = [
                "%s eta=%s"
                % ("+".join(cyclic_partition(tmp) for tmp in reduced_nodes[i]), eta[i],)
                for i in range(new_n)
            ]
            reduced_lattice.eta = eta  # one value per node
            yield reduced_lattice

    def plot(self, filename):
        """Use this function to produce an image file of the graph.

        It uses pydot library to talk to GraphViz to do this. The file
        type is determined from the filename extension (e.g. PNG, PDF).
        """
        n = self.n
        if pydot is None:
            sys.stderr.write(
                "Please install graphviz, pydot and pyparsing to draw graphs\n"
            )
            return

        # Use the OS path module's "Split Extension" function,
        extension = os.path.splitext(filename)[-1]
        # Check that there was a leading "."
        assert extension[0] == os.path.extsep
        # Remove the leading "." and make lower case (to give to graphviz as format):
        extension = extension[1:].lower()

        # Following seems to work, but we lose the number of connections
        # i.e. get a single line between nodes even when have a 2, or 3 in array
        # g=pydot.graph_from_adjacency_matrix(self.matrix, directed=directed)

        g = pydot.Dot(graph_type="digraph")

        try:
            captions = self.captions
            assert len(set(captions)) == self.n, captions
        except AttributeError:
            captions = [str(i + 1) for i in range(n)]

        for i in range(n):
            g.add_node(pydot.Node(captions[i]))

        # Test colours and styles are used when drawing graphs with multiple
        # edge types
        COLORS = [
            "black",
            "red",
            "blue",
            "green",
            "magenta",
            "pink",
            "yellow",
            "navy",
            "sienna",
            "brown",
            "crimson",
            "cyan",
        ]
        STYLES = ["solid", "dashed", "dotted", "bold"]
        LINE_COLOR_STYLES = list(zip(COLORS, STYLES * 3))

        if len(self.matrices) > len(LINE_COLOR_STYLES):
            raise ValueError(
                "Too many edge types, define some more valid "
                "graphviz styles or colors!"
            )

        # Now that we have added all the nodes, we can do the edges:
        for matrix, (color, style) in zip(self.matrices, LINE_COLOR_STYLES):
            for EdgeTo in range(n):
                for EdgeFrom in range(n):
                    # We allow multiple edges of the same type
                    for _ in range(matrix[EdgeTo, EdgeFrom]):
                        g.add_edge(
                            pydot.Edge(
                                captions[EdgeFrom],
                                captions[EdgeTo],
                                style=style,
                                color=color,
                            )
                        )

        # Now get the graph object to write itself to an image file.
        try:
            g.write(filename, prog="dot", format=extension)
        except pydot.InvocationException:
            sys.stderr.write("Please check graphviz is installed\n")


print("Runing self-tests...")

tests = doctest.testmod()
if tests.failed:
    raise RuntimeError("%i/%i tests failed" % tests)

##########################################################
# User editable section below
##########################################################

##########################################################
# Reduction examples from the 2021 manuscript and some extra
##########################################################


# Example 6.1 Unique reduced lattice (in the 2021 manuscript)
network = CoupledCellNetwork(
    [
        [0, 0, 0, 2],
        [0, 0, 0, 2],
        [0, 1, 0, 1],
        [0, 1, 0, 1],
    ]
)


# Example 6.2 Multiple reduced lattices
# network = CoupledCellNetwork(
#    [
#        [0, 0, 0, 2],
#        [0, 0, 0, 2],
#        [0, 0, 0, 2],
#        [0, 0, 1, 1],
#    ]
# )

# Example 6.3 Counter example
# network = CoupledCellNetwork(
#    [
#        [0, 1, 0, 1, 0],
#        [1, 0, 0, 0, 1],
#        [1, 0, 0, 1, 0],
#        [1, 0, 1, 0, 0],
#        [1, 1, 0, 0, 0],
#    ]
# )

# 4-cell regular network - partition lattice of 4 elements
# network = CoupledCellNetwork(
#    [
#        [0, 0, 0, 1],
#        [0, 0, 0, 1],
#        [0, 0, 0, 1],
#        [0, 0, 0, 1],
#        ]
# )

# 5-cell regular network - partition lattice of 5 elements
# network = CoupledCellNetwork(
#    [
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 0, 1],
#    ]
# )


# 7 node lattice post reduction
# network = CoupledCellNetwork(
#    [
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 1, 0],
#        [0, 1, 0, 0, 0],
#        [0, 0, 0, 1, 0],
#    ]
# )


# 9 node lattice post reduction, 3 unique?
# network = CoupledCellNetwork(
#    [
#        [0, 0, 0, 0, 1],
#        [0, 0, 0, 0, 1],
#        [0, 1, 0, 0, 0],
#        [1, 0, 0, 0, 0],
#        [0, 0, 0, 0, 1],
#    ]
# )


go(network, "reduction_test", reduce=True)

print("=" * 50)
print("Done")
