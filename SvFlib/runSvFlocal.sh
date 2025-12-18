#!/bin/bash
python --version
pyomo --version
echo $1
echo  "python $HOME/git_work/SvFlib/_START.py"
python $HOME/git_work/SvF/SvFlib/_START.py $1
