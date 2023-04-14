#!/usr/bin/env bash
docker run -p 8888:8888 -it --rm -v ./src:/tf/swe  tensorflow/tensorflow:latest-gpu-jupyter
   
