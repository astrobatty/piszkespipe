[![Image](https://img.shields.io/badge/tutorials-%E2%9C%93-blue.svg)](https://github.com/astrobatty/piszkespipe/tree/main/examples)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/astrobatty/piszkespipe/blob/main/examples/run_piszkespipe.ipynb)
<!--- [![Image](https://img.shields.io/badge/arXiv-1909.00446-blue.svg)](https://arxiv.org/abs/1909.00446) -->

# piszkespipe - A pipeline for reducing echelle spectra obtained in Piszkesteto

This code is built upon [ceres](https://github.com/rabrahm/ceres), but the two are not identical. The code is compatible with python 3.x and has several new/improved features (see below), which are not part of `ceres`. Nonetheless, the technical background is the same, which is published in [Brahm et al., 2017,PASP,129,4002](https://ui.adsabs.harvard.edu/abs/2017PASP..129c4002B/abstract).

This code is mainly written in python, but massive calculations are in C and fortran. Running ``piszkespipe`` does not requires the knowledge of python, as it is used from the command line.

## Installation dependencies

This code depends on the GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/), which can be installed on linux via
```bash
sudo apt-get install libgsl-dev
```
or on OSX via
```bash
sudo port install gsl
```

## Installation 1: in a normal python environment (for advanced python users)

Before installation, make sure GSL is properly installed, e.g. running the followings, which should return the location of the library and the installed version of GSL:
```bash
pkg-config --cflags --libs gsl && gsl-config --version
```

piszkespipe works with python 3.x, to install it ``numpy`` must be preinstalled:
```bash
pip install numpy
```

Then the package can be installed:
```bash
pip install git+https://github.com/astrobatty/piszkespipe.git
```

## Installation 2: in a separate conda environment (for people afraid of python)

The easiest way to separate ``piszkespipe`` from other python packages is to install it in a [conda](https://docs.conda.io/en/latest/miniconda.html) environment:
```bash
conda create -n piszkespipe python=3.7
```

To activate the environment:
```bash
conda activate piszkespipe
```

Then the package can be installed separately:
```bash
pip install git+https://github.com/astrobatty/piszkespipe.git
```

If the code is not used, the environment can be deactivated:
```bash
conda deactivate
```

## (No) Installation 3: running on Google Colab (for people who hates python)

You can find a Google Colab friendly tutorial [in the examples](https://github.com/astrobatty/piszkespipe/tree/master/examples/run_piszkespipe.ipynb), which can be used without installing ``piszkespipe`` and python on a local machine.


## Running the pipeline

After installation, ``piszkespipe`` will be available from the command line:

 - ``piszkespipe <path to directory> [options]``

 Listed below are the usage instructions:

```bash
$ piszkespipe --help

usage: piszkespipe [-h] [-dirout DIROUT] [-do_class] [-just_extract]
                   [-o2do O2DO] [-reffile REFFILE] [-npools NPOOLS] [-marsch]
                   [-nocosmic] [-avoid_plot]
                   directorio

piszkespipe reduces echelle spectrum obtained in Piszkesteto.

positional arguments:
  directorio        Path to the raw data.

optional arguments:
  -h, --help        show this help message and exit
  -dirout DIROUT    Path to the directory where the pipeline products will be placed. The default path will be a new directory with the same name that the input directory but followed by a '_red' suffix.
  -do_class         This option will enable the estimation of atmospheric parameters.
  -just_extract     If activated, the code will not compute the CCF and atmospheric parameters.
  -o2do O2DO        If you want to process just one particular science object you have to enter this option followed by the filename of the object.
  -reffile REFFILE  Name of the auxiliary file. The default is './reffile.txt', a file located in the directory where the raw data is.
  -npools NPOOLS    Number of CPU cores to be used by the code. Default is all.
  -marsch           If enabled, Marsch optimized raw extraction algorithm will be saved.
  -nocosmic         If activated, no cosmic ray identification will be performed.
  -avoid_plot       If activated, the code will not generate a pdf file with the plot of the computed CCF.
```

## Output

## Spectral Classification Module

## Contributing
Feel free to open PR / Issue.

## Citing
If you find this code useful, please cite [Brahm et al., 2017,PASP,129,4002](https://ui.adsabs.harvard.edu/abs/2017PASP..129c4002B/abstract) and my name to the list of co-authors. Here is the BibTeX source of Brahm et al.:
```
@ARTICLE{2017PASP..129c4002B,
       author = {{Brahm}, Rafael and {Jord{\'a}n}, Andr{\'e}s and {Espinoza}, N{\'e}stor},
        title = "{CERES: A Set of Automated Routines for Echelle Spectra}",
      journal = {\pasp},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2017,
        month = mar,
       volume = {129},
       number = {973},
        pages = {034002},
          doi = {10.1088/1538-3873/aa5455},
archivePrefix = {arXiv},
       eprint = {1609.02279},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2017PASP..129c4002B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Acknowledgements
This project was made possible by the funding provided by the Lendület Program of the Hungarian Academy of Sciences, project No. LP2018-7.
