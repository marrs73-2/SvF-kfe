# -*- SvF.ing: UTF-8 -*-
import sys
import platform
import os
import numpy as np
import COMMON as SvF

prog_name = sys.argv[0]
SvF.path_SvF_Lib            = prog_name[: prog_name.rfind('/')]             #   /home/sokol/D/SvF/SvFlib
SvF.path_SvF                = SvF.path_SvF_Lib[: SvF.path_SvF_Lib.rfind('/')+1]     #   /home/sokol/D/SvF/

SvF.tmpFileDir      = SvF.path_SvF + 'TMP/'
SvF.token           = "/home/marrs73/codingIttp/MySvF/python-api/.token30d" #SvF.path_SvF + "pyomo-everest/python-api" +'/.token'
SvF.startDir        = os.getcwd()
if platform.system() == 'Windows':   SvF.platform = 'Win'       # 2022.05

print ( SvF.startDir )
sys.path.append( SvF.startDir )

# if SvF.startDir.find('/home/vladimirv/mc2/agent') == 0 :       #  опции для       svf-remote   ************** #kfe_changed
#     SvF.DrawMode         = 'File'
#     SvF.LocalSolverName  = '/opt/scipopt911/bin/ipopt'
#     SvF.SolverName       = '/opt/scipopt911/bin/ipopt'  # 3.14.09

from ReadMng import ReadMng
from Compare import compare_results

# 
while (1) :
    SvF.Compile = True

    # Read instructions from .mng file for 1 task (till EoTask) and generate pyomo driven
    # StartModel.py file for this task.
    Task = ReadMng ( )
    print ('#######################################################################')
    SvF.Compile = False
    print ('\n\nCWD', os.getcwd(),  'RUN   StartModel.py *****************' )
    sys.path.append( os.getcwd() )

    # Run SvF throught execution of StartModel.py
    exec(open("StartModel.py").read())

    # Write results of SvF in .res file
    with open(SvF.resF, 'a') as f:
        f.write('addStrToRes: ' + SvF.addStrToRes)

    # Reset values of task dependent configurable variables
    if SvF.EofTask:
        print('\n\n\n *********  END OF TASK! **************')
        SvF.SModelFile = None
        SvF.ModelBuf = None
        SvF.resF = ''
        SvF.optEstim = sys.float_info.max
        SvF.currentTab = None
        SvF.useNaN = True    #      26.02.01  False
        ValidationSets = []
        notTrainingSets = []  # notTrainingSets содержит кого выбрасываем
        TrainingSets = []  # TrainingSets содержит точки обучени
        SvF.CV_NoRs = []
        SvF.numCV = -1
        SvF.OptNames = [] #kfe_added
        SvF.Penalty = [] # kfe_added

    # Stop reading of .mng file if it ended or EoF flag was met
    else :  
        # If COMPARE: section was used in .mng start launch results comparison of calculated tasks
        print(f"SvF.ResFilesToCompare={SvF.ResFilesToCompare}")
        if SvF.ResFilesToCompare != []:
            compare_results()
        break

print ('END OF FILE!')
exit(0)

