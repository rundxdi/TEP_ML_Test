# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import TEP
from miplearn import LearningSolver
from miplearn import LazyConstraintsComponent
import time



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
                         'pglib-opf-master/pglib_opf_case118_ieee.m',
                         'pglib-opf-master/pglib_opf_case500_goc.m',
                         'pglib-opf-master/pglib_opf_case793_goc_tep.m',
                         'pglib-opf-master/pglib_opf_case500_goc_tep.m',
                         'pglib-opf-master/pglib_opf_case2000_goc_tep.m']
    trainingInstances = []
    testInstances = []
    
    
    solver = LearningSolver(solver = 'gurobi')
    comp = LazyConstraintsComponent()
    for instance in trainingInstanceFiles:
        trainingInstances.append(TEP.TEPInstance(instance))
    for instance in testInstanceFiles:
        testInstances.append(TEP.TEPInstance(instance))

    for instance in trainingInstances:
        solver.solve(instance)

    #solver.fit(trainingInstances)

    comp.fit(trainingInstances)

    for instance in testInstances:
        t = time.time()
        solver.solve(instance)
        print(time.time() - t)
        print()


    #comp.fit(trainingInstances)
    #ev = comp.evaluate(testInstances)
    #import pandas as pd
    #pd.DataFrame(ev).mean(axis=1)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
