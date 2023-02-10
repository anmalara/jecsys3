#!/bin/bash
make clean
make -j 8
./GlobalFit jsons/info.json
