### TcrDataComparison
Master thesis (UiO) repo for TCR data comparison.

The `data` directory contains experimental mouse Igm data
obtained using sequencing. The `0_human_alpha_100seq` and
`1_human_alpha_20000seq` directories inside the `data` directory,
contain 100 and 20000 CDR3 sequences respectively that were
created (simulated) using the
[IGoR](https://github.com/qmarcou/IGoR) package.

IGoR models (from version 1.3.0) are stored in `models` and
are used to generate the CDR3 sequences.

**Note:** in order to use the Python scripts in `python`, make
sure to have installed pygor module from IGoR. This module can
be installed locally with pip by using command
`pip install ./pygor` from within the IGoR source directory.