# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import TEP
from miplearn import LearningSolver
from miplearn.instance.base import Instance
from miplearn.solvers.learning import InternalSolver
from miplearn.solvers.pyomo.base import BasePyomoSolver
from miplearn.types import ConstraintName
from miplearn import PrimalSolutionComponent
from miplearn import ObjectiveValueComponent
import time
import pickle



def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    trainingInstanceFiles = ['pglib-opf-master/pglib_opf_case3_lmbd.m',
                            'pglib-opf-master/pglib_opf_case5_pjm.m',
                            'pglib-opf-master/pglib_opf_case14_ieee.m',
                            'pglib-opf-master/pglib_opf_case24_ieee_rts.m',
                            'pglib-opf-master/pglib_opf_case30_as.m',
                            'pglib-opf-master/pglib_opf_case30_ieee.m',
                            'pglib-opf-master/pglib_opf_case39_epri.m',
                            'pglib-opf-master/pglib_opf_case57_ieee.m',
                            'pglib-opf-master/pglib_opf_case73_ieee_rts.m',
                            'pglib-opf-master/pglib_opf_case89_pegase.m',
                            'pglib-opf-master/pglib_opf_case118_ieee.m',
                            'pglib-opf-master/pglib_opf_case118_ieee_anor.m',
                            'pglib-opf-master/pglib_opf_case162_ieee_dtc.m',
                            'pglib-opf-master/pglib_opf_case200_activ.m',
                            'pglib-opf-master/pglib_opf_case240_pserc.m',
                            'pglib-opf-master/pglib_opf_case300_ieee.m']
    testInstanceFiles = ['pglib-opf-master/pglib_opf_case14_ieee.m',
                         'pglib-opf-master/pglib_opf_case118_ieee.m',]
                         #'pglib-opf-master/pglib_opf_case500_goc.m',
                         #'pglib-opf-master/pglib_opf_case793_goc_tep.m',
                         #'pglib-opf-master/pglib_opf_case500_goc_tep.m',
                         #'pglib-opf-master/pglib_opf_case2000_goc_tep.m']
    trainingInstances = []
    testInstances = []
    
    
    solver = LearningSolver()


    for instance in trainingInstanceFiles:
        trainingInstances.append(TEP.TEPInstance(instance))
    for instance in testInstanceFiles:
        testInstances.append(TEP.TEPInstance(instance))

    print("Now solving the training instances")
    for instance in trainingInstances:
        solver.solve(instance)

    print("Now fitting training instances to solver")
    solver.fit(trainingInstances)

    #pickle.dump(solver, open("solver.pickle", "wb"))

    for instance in testInstances:
        t = time.time()
        print("Now solving instance ", str(instance))
        stats = solver.solve(instance)
        #print(stats['lp_log'])
        #print(stats['mip_log'])
        print(time.time() - t)
        print()


    #comp.fit(trainingInstances)
    print("Now evaluating the component?")
    #ev = comp.evaluate(testInstances)
    #ev = solver.evaluate(testInstances)
    #ev = prim.evaluate(testInstances)
    #import pandas as pd
    #pd.DataFrame(ev['Fix one']).mean(axis=1)
    #ev2 = obj.evaluate(testInstances)
    #pd.DataFrame(ev2).mean(axis=1)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
