#!/bin/bash

# Script to download test data for Subgraph Analyser pipeline
wget --content-disposition https://osf.io/s75y4/download
tar -xzf test_data.tar.gz
rm test_data.tar.gz