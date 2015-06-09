#!/bin/bash

echo "compiling externals..."

chmod +x compile.ext.mar15p1.sh
./compile.ext.mar15p1.sh

echo "compiling pandaroot...."

chmod +x compile.pandaroot.sh mar15p1 mar15
./compile.pandaroot.sh mar15p1 mar15

echo "DONE"
