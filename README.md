# multichainACF

Calculates the globally-centered autcorrelation function (ACF) plots following Agarwal and Vats (2020). Given a list of multiple Markov chains, this package offers functions for ACF and CCF plots that pool information from all chains yielding more accurate estimates of the correlation. The functions differs from the base `acf` and `ccf` by allowing the users to center the chains around the global mean from all Markov chains. A different mean argument can also be specified.

See [webpage](https://home.iitk.ac.in/~dootika/_pages/mcACF.html) for a working example.
