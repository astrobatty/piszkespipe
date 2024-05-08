[![Image](https://img.shields.io/badge/tutorials-%E2%9C%93-blue.svg)](https://github.com/astrobatty/piszkespipe/tree/main/examples)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/astrobatty/piszkespipe/blob/main/examples/run_piszkespipe.ipynb)
<!--- [![Image](https://img.shields.io/badge/arXiv-1909.00446-blue.svg)](https://arxiv.org/abs/1909.00446) -->

# piszkespipe - A pipeline for reducing echelle spectra obtained in Piszkéstető

The purpose of `piszkespipe` is to automatically extract and (optionally) pre-analyze Echelle spectra obtained with the ACE spectrograph mounted on the 1m RCC telescope in Piszkéstető. The pre-analyze step includes the calculation of CCF assuming both single and double line objects (i.e. single and binary stars) and of atmospheric parameters.

This code is built upon [ceres](https://github.com/rabrahm/ceres), but the two are not identical. `piszkespipe` is compatible with python 3.x, has new/improved features (see below) and optimized for observations taken in Piszkéstető. Nonetheless, the technical background is the same, which is published in [Brahm et al., 2017,PASP,129,4002](https://ui.adsabs.harvard.edu/abs/2017PASP..129c4002B/abstract).

The pipeline is mainly written in python, but massive calculations are in C and fortran. Running ``piszkespipe`` does not require the knowledge of python; after installation it can called from the command line.

# Installation and dependencies

The pipeline can be installed via three ways, depending on your desire (see below), but each depends on the GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/), which can be installed on linux via
```bash
sudo apt-get install libgsl-dev
```
or on OSX via
```bash
sudo port install gsl
```

Before installation, make sure GSL is properly installed, e.g. running the followings, which should return the location of the library and the installed version of GSL:
```bash
pkg-config --cflags --libs gsl && gsl-config --version
```

## Installation 1: in a normal python environment (for regular python users)

``piszkespipe`` can be installed from source:
```bash
pip install piszkespipe
```

## Installation 2: in a separate conda environment (for people afraid of python)

The easiest way to separate ``piszkespipe`` from other python packages (and everything else) is to install it in a [conda (see this link for minimal installiation)](https://docs.conda.io/en/latest/miniconda.html) environment:
```bash
conda create -n piszkespipe python=3.7
```

To activate the environment:
```bash
conda activate piszkespipe
```

Then the package can be installed:
```bash
pip install piszkespipe
```

If the code is not used, the environment can be deactivated:
```bash
conda deactivate
```

## (No) Installation 3: running on Google Colab (for people who hate python)

``piszkespipe`` can be used/tested without installing python on a local machine. You can use the Google Colab friendly tutorial, located [in the examples](https://github.com/astrobatty/piszkespipe/tree/main/examples/run_piszkespipe.ipynb).

# Running the pipeline

After installation, ``piszkespipe`` will be available from the command line:

```bash
piszkespipe <path to raw data> [options]
```

 Listed below are the usage instructions:

```bash
$ piszkespipe --help

usage: piszkespipe [-h] [-dirout DIROUT] [-do_class] [-just_extract]
                   [-o2do O2DO] [-reffile REFFILE] [-npools NPOOLS] [-marsch]
                   [-nocosmic] [-keepoutliers] [-avoid_plot]
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
  -keepoutliers     If activated, strong upper (emission line like) outliers won't be removed.
  -avoid_plot       If activated, the code will not generate a pdf file with the plot of the computed CCF.
```

By deafult, to run the pipeline you have to collect your nightly/weekly raw measurements in a directory, which path must be passed to ``piszkespipe``. __Raw measurements__ are:

- bias
- dark
- fiber flat
- ThAr
- science

images. If you forgot to take flat frames, the code will use the median of all science frames to trace Echelle orders and the blaze corrected spectra will be the same as the continuum normalized one.

## Available options

Some comments on the non self-evident options.

- `marsch` optimized extraction, which details can be found in [Marsh,1989,PASP,101,1032](https://ui.adsabs.harvard.edu/abs/1989PASP..101.1032M/abstract), sometimes yields ambiguous results, thus the result of the simple extraction will be saved by default.
- `nocosmic`, Cosmic ray identification is done using the algorithm of [L.A.Cosmic](http://www.astro.yale.edu/dokkum/lacosmic/), but if you feel like e.g. some sharp emisson lines are removed, this option will deactivate the identification.
- `reffile` is an auxiliary file corresponds to a plain text file that can be used to give data for particular targets in order to improve the precision in the radial velocity estimation. The file should contain 8 comma-separated colummns containing the following information:
- `keepoutliers`, If a target has very strong emission lines, then those might be excluded. To avoid this effect, use this flag.

```
1-  filename of the target.
2-  right ascension of the target (J2000) with format hh:mm:ss.
3-  declination of the target (J2000) with format dd:mm:ss.
4-  proper motion in RA [mas/yr].
5-  proper motion in DEC [mas/yr].
6-  integer (0 or 1). If 0, the code uses the cordinates given in the image header.
    If 1, the code uses the coordinates given in this file.
7-  mask that will be used to compute the CCF. Allowed entries are G2, K5 and M2.
8-  velocity width in km/s that is used to broaden the lines of the binary mask.
    It should be similar to the standard deviation of the Gaussian that is fitted to the CCF.
```
Here is an example of a reffile:

```
betaCVn-180s-20181226-53445_0001_f1,12:33:44.5448195,+41:21:26.924857,-704.75,292.74,1,G2,4.0
```

If the reffile is not provided, the code uses the coordinates found in the image header and the proper motions obtained from Gaia catalog to compute the barycentric correction and uses the G2 mask to compute the CCF. __Note that the CCF is way more dependent on the broadening of mask lines than on the spectral type of the applied mask. The velocity width can be estimated automatically from the v*sin(i) of the best matching synthetic spectra, but to do so `-do_class` must be enabled (see below).__

### Auxiliary files

Additionally, there are other two auxiliary files that can be placed in the same directory of the raw images for improving the reductions.

A file called ``bad_files.txt`` can be used to list all the images of the raw directory that should be ignored by the pipeline. Each line of the bad_files file must have the filename of the image without the complete path.

Another file named ``moon_corr.txt`` can be used to specify the images for which the CCF computation will include a double gaussian fit in order to correct the radial velocity for the effects of scattered moonlight contamination. _In case of the double Gaussian fit mainly used for binary stars, the code does not include a third gaussian._

## Understanding the Output

While the directory specified in the `dirout` option will contain several intermediate reduction files (which are useful if you re-run the code e.g. to apply moon_corr), the final results will be placed in the subdirectory `/proc`. This directory can contain up to five types of files:

- `_final.fits` files with the extracted and wavelength calibrated spectra.
- `_combined_barycorrected.fits` files with the 1 dimensional combined and barycentric corrected spectra.
- `.pdf` files that show the CCF plots and fits.
- `_XCs_*.fits` files with the raw and (simple linear regression) corrected CCF.
- `results.csv` that contains a summary of the results of the reduction including the radial velocity measurements and the atmospheric parameters.

The `_final.fits` files are data cubes with dimensions (11,nords,npix), where _nords_ are the total number of processed echelle orders, and _npix_ is the number of pixels in the dispersion direction. The eleven entries in the first dimension correspond to:

```
0-  Wavelength
1-  Extracted Flux (Simple or Marsch if enabled)
2-  Measurement of the error in the extracted flux [1./sqrt(Var)]
3-  Blaze corrected Flux
4-  Measurement of the error in the blaze corrected flux
5-  Continuum normalized flux
6-  Measurement of the error in the continuum normalized flux
7-  Estimated continuum
8-  Signal-to-noise ratio
9-  Continumm normalized flux multiplied by the derivative of the wavelength with respect to the pixels
10- Corresponding error of the 9th entrance
```

The colummns in the `results.csv` file are the following:

```
0-  Object name
1-  Input filename
2-  BJD
3-  RV
4-  Error in RV
5-  Bisector span
6-  Error in bisector span
7-  Barycentric velocity
8-  Barycentric corrected RV
9-  Error in barycentric corrected RV
10- Systematic error in RV
11- First RV in two Gaussian fit (RV1)
11- Barycentric corrected RV1
11- Error in RV1
11- Second RV in two Gaussian fit (RV2)
11- Barycentric corrected RV2
11- Error in RV2
12- Instrument
7-  Pipeline
8-  Resolving power
9-  Effective temperture
9-  Formal error in effective temperture
10- log(g)
10- Formal error in log(g)
11- [Fe/H]
11- Formal error in [Fe/H]
12- v*sin(i)
13- Value of the continuum normalized CCF at its lowest point
14- Standard deviation of the gaussian fitted to the CCF
15- Exposure time
16- Signal to noise ratio at ~5150 \AA
17- Path to the CCF plot file
```

_Note:_ the systematic error in RV can be elliminated by measuring the RV of a standard star and comparing its current value to the literature.

## Spectral Classification Module

In order to use the automated spectral classification routines (`-do_class` option), a set of synthetic spectra is required. In particular, the original CERES uses a modified version of the [Coelho et al 2005 models](https://ui.adsabs.harvard.edu/abs/2005A&A...443..735C/abstract), where the models have been reduced in wavelength coverage, degraded in resolution and v*sin(i), and reduced in sampling. This grid of models can be downloaded from [this link](https://cloud.konkoly.hu/s/taT2qiSjpBwWWDr/download/coelho_05_red4_R40.tar.gz) (2.85 GB; or [from the original source](http://www.astro.puc.cl/~rbrahm/coelho_05_red4_R40.tar.gz)) and the $COELHO path to models must be defined:

```bash
mkdir ~/COELHO_MODELS
cd ~/COELHO_MODELS
wget https://cloud.konkoly.hu/s/taT2qiSjpBwWWDr/download/coelho_05_red4_R40.tar.gz
# Or from the orignal source:
# wget http://www.astro.puc.cl/~rbrahm/coelho_05_red4_R40.tar.gz
tar -xf coelho_05_red4_R40.tar.gz
```

Then add the location of R_40000b dir to $COELHO path, e.g. by
```bash
export COELHO=~/COELHO_MODELS
```

Optionally add this line to your `.bash_rc` file.

## Contributing
If you encounter any issues do not hesitate to open an Issue or contact me directly. Any kind of improvements are most welcome as PRs.

## Citing
If you find this code useful, please cite [Brahm et al., 2017,PASP,129,4002](https://ui.adsabs.harvard.edu/abs/2017PASP..129c4002B/abstract) and add my name to the list of co-authors. Here is the BibTeX source of Brahm et al.:
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
