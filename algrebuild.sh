#!/bin/bash
# Copyright (c) General Electric Company, 2017.  All rights reserved.

echo ""
echo "Removing algorithm image."
#docker rmi thriveitcr/multi-omics-heterogeneity-analysis

echo ""
echo "Building algorithm image."
docker build -t thriveitcr/multi-omics-heterogeneity-analysis --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy=$no_proxy .
echo ""
docker images
