---
title: "cran-comments.md"
output: html_document
---
## Resubmission
This is a resubmission. In this version I have:

* Added Makevars files in order to run on win-builder successfully

* Changed some of my C++ code to fix a couple of warnings from win-builder

* Changed license to GPL-3, combined with license from code my work is derived from.

## Test environments
* local OS X install, R 3.3.3
* ubuntu 12.04 (on travis-ci), R 3.3.2
* win-builder

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs

## Downstream dependencies
This is the first version of bcs.
