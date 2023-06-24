#!/usr/bin/env bash
docker run --gpus all -p 8888:8888 -it --rm -v ./src:/tf/swe  tensorflow/tensorflow:latest-gpu-jupyter
   
