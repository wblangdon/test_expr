# test_expr
Create random inputs for unix expr and answer it should give. 
Uniformly randomly samples binary trees and then converts them to random arithmetic expressions.

Binary trees of a given size can be sampled rapidly uniformly at random.

![video](reading_group_bigtree/tree_square.gif)
A binary is tree is created by random vertical and horizontal moves across square,
starting at a corner and moving towards the diagnoally opposite corner,
keeping above the diagonal [slides](reading_group_bigtree/langdon_sse_15-jan-2020.pdf)
[also](http://www.cs.ucl.ac.uk/staff/W.Langdon/gggp/bigtree/langdon_sse_15-jan-2020.pdf)
The size of the tree is given by the length of the square's sides.
Random vertical and horizontal moves ensure all possible binary trees of the chosen size
are equally likely.
Works for any size tree.
Fast: bench_rand 18,600,000 tree nodes per second.
(bench_rand based on GPquick).

bench_rand generates a random binary tree,
which is piped into gawk script test_expr.awk
which converts the tree into an arithmetic expression
which can be used to test the unix expr utility.
