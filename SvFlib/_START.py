# -*- SvF.ing: UTF-8 -*-
import sys
import platform
import os
# import numpy as np
import COMMON as SvF


if platform.system() == 'Windows':   SvF.platform = 'Win'       # 2022.05

from ReadMng import ReadMng

while (1) :
    SvF.Compile = True

    # Read instructions from .mng file for 1 task (till EoTask) and generate pyomo driven
    # StartModel.py file for this task.
    EoF = ReadMng ( )

    print ( '#########  Start StartModel.py in ' + os.getcwd() + ' #########\n' )
    import subprocess

    # Run SvF throught execution of StartModel.py in a separate subprocess
    result = subprocess.run(
        [sys.executable, os.getcwd()+"/StartModel.py"],
        check=True
    )

    # Stop reading of .mng file if it ended or EoF flag was met
    if EoF == 'EOF' :
        break
    else:
        print('\n\n\n *********  END OF TASK! **************')


print ('END OF FILE!')
exit(0)

