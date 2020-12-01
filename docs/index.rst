.. ARTS tracking beams documentation master file, created by
   sphinx-quickstart on Tue Dec  1 10:12:35 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ARTS tracking beams
===============================================
The Apertif Radio Transient System (ARTS) archive contains tied-array beam (TAB) data. The TABs have a time-dependent and
frequency-dependent pointing. This tool is able to convert the TAB data to a tracking beam (TB), which tracks a fixed point
on the sky over the course of an observation. Additionally, it can convert TAB data to Synthesised Beams (SBs),
which are suitable for transient searches.

Installation and basic usage instructions are available at the project's
`Github page <https://github.com/loostrum/arts_tracking_beams>`_.

This page contains a tutorial for both the tracking beams and synthesised beams. The required data can be downloaded
from the `Apertif Long Term Archive (ALTA) <https://alta.astron.nl>`_. 
The `arts_tools <https://github.com/loostrum/arts_tools>`_ package provides an automated way to do so, and will be used in the tutorials.

The tracking beams track a single RA/Dec point, and are suitable for known pulsar or periodicity searches.
The synthesised beams are broad-band and have a larger field-of-view than the tracking beams, but they have a time-dependent pointing.
Hence, they are suitable for single-pulse searches or very short periodicity searches only.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   tutorials/TB
   tutorials/SB
