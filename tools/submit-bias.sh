#!/bin/bash

Rscript -e "targets::tar_make_future(script='_tfcebias.R', store='_bias', workers=10)"
# Rscript -e "targets::tar_make(tfce, script='_tfcebias.R', store='_bias')"
