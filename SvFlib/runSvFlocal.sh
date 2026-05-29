#!/bin/bash
python --version
pyomo --version
echo $1
echo  "python $HOME/codingIttp/MySvF/SvF-kfe/SvFlib/_START.py"
python $HOME/codingIttp/MySvF/SvF-kfe/SvFlib/_START.py $1
