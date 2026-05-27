#!/bin/bash
make clean
make -j 8
./GlobalFit jsons/info.json
/opt/local/bin/python3.14 python/PlotGlobalFit.py