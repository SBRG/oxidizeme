#============================================================
# file make_stressme.py
#
# Make stressME model from cobrame model and stress 
# knowledgebase
#  
# Laurence Yang, SBRG, UCSD
#
# 11 Mar 2019:  port from tests
#============================================================

#from stressme.tests.test_construct import TestConstruct
from oxidizeme.model import StressME
from qminospy.me1 import ME_NLP1
from ecolime import build_me_model
import cloudpickle
import argparse
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parameters.')
    parser.add_argument('--savepath', metavar='savepath', type=str, default='')
    args = parser.parse_args()

    savepath = args.savepath

    if savepath=='':
        savepath='OxidizeME.pickle'

    print("Will save StressME model to file: %s"%savepath)

    #--------------------------------------------------------
    # Build ME model
    print('#--------------------------------------------------------')
    print("Making basic ME model")
    me = build_me_model.return_me_model()
    #********************************************************
    # print('#--------------------------------------------------------')
    # print("Load basic ME model")
    #********************************************************

    ## Save to file
    #me_path = 'me_nostress.pickle'
    #with open(me_path, 'wb') as f:
    #    cloudpickle.dump(me, f)
    
    #--------------------------------------------------------
    # Build StressME model
    print('#--------------------------------------------------------')
    #print("Making and saving stressME")
    print("Making stressME")
    #construct = TestConstruct()
    #construct.test_construction(me_path, savepath)
    solver = ME_NLP1(me)
    stress = StressME(solver)

    # Make stressME
    stress.make_stressme()

    # Force demetallation and mismetallation via coupling
    stress.force_demetallation(csense='G')
    stress.force_mismetallation(csense='G')

    # Force Fe-S cluster damage
    stress.force_fes_damage(csense='G')

    # Need to open some bounds
    print("Opening some bounds")
    rxns_open = ['SPODM', 'SPODMpp', 'AHPRED', 'CAT']
    for rid_open in rxns_open:
        try:
            data = me.stoichiometric_data.get_by_id(rid_open)
            for rxn in data.parent_reactions:
                rxn.upper_bound = 1000.
        except:
            pass

    # Ensure detoxification proteins can be expressed
    print("Allowing detox fluxes")
    detox_prots = ['b3908', 'b1656', 'b1646']
    try:
        for bnum in detox_prots:
            rxn = me.reactions.get_by_id('translation_'+bnum)
            rxn.upper_bound = 1000.
    except:
        pass

    #--------------------------------------------------------
    # Test solve
    print("Setting the ROS and metal concentrations")
    stress.substitute_ros()
    stress.substitute_metal()

    print("Test solving stressME")
    sol = solver.bisectmu(1e-4)
