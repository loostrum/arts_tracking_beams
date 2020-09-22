# ARTS tracking beams
[![DOI](https://zenodo.org/badge/291081079.svg)](https://zenodo.org/badge/latestdoi/291081079)
[![PyPI version](https://badge.fury.io/py/arts-tracking-beams.svg)](https://badge.fury.io/py/arts-tracking-beams)
[![Build Status](https://travis-ci.com/loostrum/arts_tracking_beams.svg?branch=master)](https://travis-ci.com/loostrum/arts_tracking_beams)
[![codecov](https://codecov.io/gh/loostrum/arts_tracking_beams/branch/master/graph/badge.svg)](https://codecov.io/gh/loostrum/arts_tracking_beams)

The Apertif Radio Transient System (ARTS) archive contains tied-array beam (TAB) data. The TABs have a time-dependent and
frequency-dependent pointing. This tool is able to convert the TAB data to a tracking beam (TB), which tracks a fixed point
on the sky over the course of an observation. 

## Dependencies
* python >= 3.6
* numpy >= 1.17
* astropy
* tqdm

## Installation
To install the latest release:

`pip install arts_tracking_beams`
  
To install the latest master branch:
 
`pip install git+https://github.com/loostrum/arts_tracking_beams`

## Usage

### Input data
First download the data set of interest from the Apertif Long-Term Archive (ALTA). Tools to find which pulsars are in the 
field-of-view of a given Apertif pointing and to download the data are available as a separate
[python package](https://github.com/loostrum/arts_tools).

A data file from the archive is identified by three parameters: the task ID, compound beam (CB) index, and TAB index.
The file `ARTS200102003_CB00_TAB00.fits` would be the observation identified by task ID 200102003
(that is, the third observation on January 2nd, 2020), CB zero, TAB zero. A TB is created from the TABs of a single CB.

### Creating a tracking beam
The TB is created from the TAB data with `arts_create_tracking_beam`. 

The simplest use case is to create a tracking beam
from a folder which contains only one data set (i.e. the TABs of one CB of one observation), for a source with known 
coordinates. For example, to create a tracking beam towards the Crab pulsar:

`arts_create_tracking_beam --input_folder /path/to/data/ --source 'PSR B0531+21'`

If there are multiple data sets in the input data folder, specify the task ID and/or CB index. Instead of the source name,
it is also possible to provide a RA and Dec. The name of the output FITS file is determined automatically from the input 
source name or RA/Dec, but can also be specified manually. Using all of these options, an example command is:

`arts_create_tracking_beam --input_folder /path/to/data/ --taskid 200102003 --cb 0 --ra 05:34:32 --dec 22:00:52 --output tracking_beam.fits`

The TB creation consists of two steps:
1. Calculate the required TABs at each frequency and time
2. Reorder the data from the input TAB FITS files and create a new FITS file containing the TB.

The results of step 1 can be saved to disk with `--save_tab_indices`. To only calculate the TAB indices and 
disable step 2 completely, use `--no_fits_output`. 
To generate the FITS output from a TAB indices file on disk, use`--load_tab_indices /path/to/tab/index/file.txt`.
The script then loads the TAB indices and immediately goes to step 2.

There are a few more settings that can be customized. Run `arts_create_tracking_beam -h` for an overview of all options.

### Creating a synthesised beam
A synthesised beam (SB) is a type of beam that reorders the TABs as function of frequency, but *not* as function of time.
A single CB is covered by 71 SBs. Each SB is always made out of the same TABs. The SBs are used in the real-time 
transient search that ARTS runs. The brightest transients may also be detectable in the archival data, so we here include
a tool to create the synthesised beams as well. 

The synthesised beam tool, `arts_create_synthesised_beam`, works
in a very similar fashion as the tracking beam tool. An example command:
 
`arts_create_synthesised_beam --input_folder /path/to/data --sb 35`

Run `arts_create_synthesised_beam -h` for more options.


