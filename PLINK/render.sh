#!/usr/bin/env bash

R -e "rmarkdown::render('runPlink.Rmd', output_format='all')"
