# inferchange: Multiscale Covariance Scanning and Related Algorithms

The aim of the `inferchange` package is to make methods for inference of
changes in high-dimensional linear regression available to data analysts and
researchers in statistics.

You can track (and contribute to) the development of `inferchange` at https://github.com/tobiaskley/inferchange. If you encounter unexpected behavior
while using `inferchange`, please let us know by writing an [email](mailto:tobias.kley@uni-goettingen.de) or by filing an [issue](http://github.com/tobiaskley/inferchange/issues).

Currently, the methodology described in the following pre-print is implemented:

* Cho, H., Kley, T., and Li, H. (2024). Detection and inference of changes in
  high-dimensional linear regression with non-sparse structures
  ([arXiv](http://arxiv.org/abs/2402.06915)).


## Getting started with ``inferchange``

First, if you have not done so already, install R from http://www.r-project.org
(click on download R, select a location close to you, and download R for your
platform). Once you have the latest version of R installed and started execute
the following commands on the R shell:

 ```
 install.packages("devtools")
 devtools::install_github("tobiaskley/inferchange")
 ```

This will first install the R package ``devtools`` and then use it to install
the latest (development) version of ``inferchange`` from the GitHub repository.

Now that you have R and ``inferchange`` installed you can access all the
functions available. To load the package and access the help files:

```
library(inferchange)
help("inferchange")
```

At the bottom of the online help page to the package you will find an index to
all the help files available. The main functions are ``McScan``, ``lope``,
``clom`` and ``ci_delta``. The respective help pages can be accessed by

```
help("McScan")
help("lope")
help("clom")
help("ci_delta")
```

A "workhorse" function, named ``inferchange``, that wraps the steps of a full
change point analysis is also available. To access the help page call

```
help("inferchange")
```
