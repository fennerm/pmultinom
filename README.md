# pmultinom

You're probably here because you've noticed that the multinomial cumulative
distribution function is absent from the core libraries. This package fills that
gap. It's implemented in C++ for efficiency and syntax you're probably familiar
with.

## Installation

pmultinom depends on the GNU scientific library (GSL), if you have trouble with
the installation, double check that it's installed and set up correctly on your
system.

### Debian (Ubuntu/Linux Mint)

```
  sudo apt-get install gsl-bin libgsl0-dev
```

Then from within R:

```
   install.packages('devtools')
   devtools::install_github('fennerm/pmultinom')
```

### MacOS

```
   brew install gsl
```

Then from within R:

```
   install.packages('devtools')
   devtools::install_github('fennerm/pmultinom')
```

### Other GNU/Linux

Download the gsl development library for your distro, then install with devtools
as above.

## Usage

pmultinom exports a single function:

```
    pmultinom(q, n, prob, log.p = FALSE, lower.tail = TRUE)
```

##### Arguments

```
       q: Vector of quantiles; If q has length 1 then the same quantile
          is used for each bin. If q is length k, each bin is assigned
          its own quantile

       n: Total number of balls

    prob: A vector of length k, describing the probability of a ball
          being placed into bin xk

lower.tail: logical; if TRUE (default), probabilities are P[X <= x],
          otherwise, P[X > x].
```

##### Examples

For 4 balls randomly thrown into 2 bins, probability that all bins end up with 2
or fewer balls

```
pmultinom(q = 2, n = 4, prob = c(0.5, 0.5))
[1] 0.3361880568`
```

For 6 balls randomly thrown into 3 bins, probability that at least one bin ends
up with 3 or more balls

```
pmultinom(q = 2, n = 6, c(1/3, 1/3, 1/3))
[1] 0.1142714833
```

Same as above, but with different probabilities for each bin

```
pmultinom(q = 2, n = 6, c(0.2, 0.3, 0.5))
[1] 0.0749668976`
```

### Credit

- Glenn Stone (staff.scm.uws.edu.au/~glenn) for the original C implementation.
