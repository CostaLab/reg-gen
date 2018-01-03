#!/usr/bin/env bash

################################################################
# Test for RGT tools
################################################################

export RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

# make any error stop the script
set -e

SCRIPT_DIR=$(dirname "$0")

# THOR
bash ${SCRIPT_DIR}/thor.sh

# TDF
bash ${SCRIPT_DIR}/tdf.sh

# Viz
bash ${SCRIPT_DIR}/viz.sh

# Motif Analysis
bash ${SCRIPT_DIR}/motif.sh

# HINT
bash ${SCRIPT_DIR}/hint.sh

