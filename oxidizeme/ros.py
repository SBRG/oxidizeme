#============================================================
# File ros.py
# 
# Class for adding spontaneous ROS generation.
#
# Laurence Yang, SBRG, UCSD
#
# 21 Feb 2018:  refactored for publication.
#
# Changes from Brynildsen et al. (2013) while porting to ME
# 1. M metabolites that are actually complexes changed to ME complex:
#    fldox, fldrd, trdox, trdrd 
# 2. Putative h2o2 and o2s rxns for FLDR should have fldox and fldr in the rxns.
#    Corrected.
# 3. ME has FLDR21 and FLDR22 (2 fldox per nadph) instead of FLDR (1 fldox/nadph).
# 4. ME has RNTRx21 and RNTRx22 (x=1 to 4). Also 2 fldox per atp/gtp/ctp/utp.
#               
#============================================================

from cobra import Model, Reaction, Metabolite
from sympy import Symbol, Basic, sympify
import sympy

import numpy.random as rand
import pandas as pd
import numpy as np
import warnings
import os

import oxidizeme

def replace_last(src, old, new):
    head, sep, tail = src.rpartition(old)
    return head + new + tail


class SpontROS(object):
    """
    Class for adding spontaneous ROS generation.
    Port of JM Monk's add_ros.py to ME models.
    """
    def __init__(self, me, ros_coeff_dict={'h2o2_c':'c', 'o2s_c':'k'}):
            
        self.me = me
        self.ros_coeff_dict = ros_coeff_dict
        # Placeholder for spont ros rxn list
        self.ros_rxns = []
        self.ros_stoichs = []   # ditto for ros stoichiometric_data
        self.stoich_coeff_dict = {}    # map rxn to c & k coeffs

    def pathto(self, filename):
        return os.path.join(oxidizeme.__path__[0], 'data', filename)

    
    def add_ros_rxns(self, rxn_file='spont_ros_rxns_me2.csv',
            col_name='Rxn Name',
            col_h2o2='Putative h2o2 Rxn', col_o2s='Putative o2s Rxn',
            met_repl_dict={'[c]':'_c', '[p]':'_p',
                '-L':'__L', '-D':'__D', '-S':'__S'},
            complex_file='protein_complexes.txt',
            verbosity=1):
        """
        add_ros_rxns(self, in_model, cvals, kvals, name='')

        Adds spontaneous ROS generation to stoichiometric_data
        stoichiometry for each of Brynildsen's ROS-generating reactions.
        """
        me = self.me
        ros_coeff_dict = self.ros_coeff_dict
        sym_h2o2 = ros_coeff_dict['h2o2_c']
        sym_o2s  = ros_coeff_dict['o2s_c']
    
        # Temp model object that holds all the ROS reactions
        mdl = Model('ros')

        # Exclude compartment from complexes
        df_cplx = pd.read_csv(self.pathto(complex_file), sep='\t', header=None)
        complexes = df_cplx[0].tolist()

        # Load Brynildsen's reactions
        df_ros = pd.read_csv(self.pathto(rxn_file)).dropna(subset=[col_h2o2, col_o2s])

        stoichs_mod = []

        for irow,row in df_ros.iterrows():
            rxn_str = row[col_h2o2]
            rxn_name = row[col_name]
            if not me.stoichiometric_data.has_id(rxn_name):
                print rxn_name, 'not in stoichiometric_data. Skipping.'
                continue

            stoich = me.stoichiometric_data.get_by_id(rxn_name)
            stoichs_mod.append(stoich)
            
            rxn_h2o2 = Reaction(stoich.id+'_h2o2')
            rxn_o2s  = Reaction(stoich.id+'_o2s')
            mdl.add_reactions([rxn_h2o2, rxn_o2s])

            # Get h2o2 and o2s rxns
            rxn_strs = row['Putative h2o2 Rxn']
            # Sometimes the reaction string has multiple rxns separated by commas
            for rxn_str in rxn_strs.split(','):
                rxn_h2o2.build_reaction_from_string(rxn_str, verbose=False)

            rxn_strs = row['Putative o2s Rxn']
            for rxn_str in rxn_strs.split(','):
                rxn_o2s.build_reaction_from_string(rxn_str, verbose=False)

            # Keep c & k stoichiometries for h2o2 and o2s as Symbols,
            # substituted at time of simulation, like mu.
            # Not too many so probably don't need to compile expressions like we do for mu
            # Note: different c & k for each rxn
            s_h2o2 = Symbol(sym_h2o2+'_'+str(irow), positive=True)
            s_o2s  = Symbol(sym_o2s+'_'+str(irow), positive=True)
            for met,s in rxn_h2o2.metabolites.iteritems():
                rxn_h2o2._metabolites[met] = s*s_h2o2

            for met,s in rxn_o2s.metabolites.iteritems():
                rxn_o2s._metabolites[met] = s*s_o2s

            # Fix compartments and other ID inconsistencies
            for met in rxn_h2o2._metabolites.keys():
                for old, new in met_repl_dict.iteritems():
                    met.id = met.id.replace(old, new)
                    # No compartment for complexes
                    id2 = replace_last(met.id, new, '')
                    id3 = id2.split('_mod')[0]  # modified complex not in complex list
                    if id3 in complexes:
                        met.id = id2

            for met in rxn_o2s._metabolites.keys():
                for old, new in met_repl_dict.iteritems():
                    met.id = met.id.replace(old, new)
                    # No compartment for complexes
                    id2 = replace_last(met.id, new, '')
                    id3 = id2.split('_mod')[0]  # modified complex not in complex list
                    if id3 in complexes:
                        met.id = id2

            #------------------------------------------------
            # Ensure that metabolite IDs in ros rxns consistent with ME
            # If not, stop now.
            for met in rxn_h2o2.metabolites:
                if not me.metabolites.has_id(met.id):
                    raise Exception(met.id+' not in me model. Check that ids are consistent.')

            for met in rxn_o2s.metabolites:
                if not me.metabolites.has_id(met.id):
                    raise Exception(met.id+' not in me model. Check that ids are consistent.')
            #------------------------------------------------

            # Add the original (pre-ROS) stoichiometry to ROS portion
            # Option 1: alter the stoichiometric_dat and push to parent_reactions
            # Option 2: add stoich to parent_reactions
            for rxn in stoich.parent_reactions:
                # Make sure that ROS only generated, not consumed spontaneously
                if rxn.reverse:
                    if verbosity>0:
                        print 'Skipping reverse rxn:', rxn.id

                else:
                    rxn.add_metabolites(rxn_h2o2.metabolites)
                    rxn.add_metabolites(rxn_o2s.metabolites)
                    self.ros_rxns.append(rxn)
                    self.stoich_coeff_dict[stoich.id] = {'ind':irow, 'h2o2_c':s_h2o2, 'o2s_c':s_o2s}
            
        # Fix compartments and other ID inconsistencies
        mdl.repair()

        self.ros_stoichs = stoichs_mod
        
        return stoichs_mod, mdl


    def get_ros_rxns(self, rxn_file='spont_ros_rxns_me2.csv',
            col_name='Rxn Name',
            col_h2o2='Putative h2o2 Rxn', col_o2s='Putative o2s Rxn',
            verbosity=1):
        """
        Get ROS rxns but don't alter stoichiometry
        """
        me = self.me
        # Temp model object that holds all the ROS reactions
        mdl = Model('ros')
        # Load Brynildsen's reactions
        df_ros = pd.read_csv(self.pathto(rxn_file)).dropna(subset=[col_h2o2, col_o2s])

        ros_rxns = []

        for irow,row in df_ros.iterrows():
            rxn_name = row[col_name]
            if not me.stoichiometric_data.has_id(rxn_name):
                continue

            stoich = me.stoichiometric_data.get_by_id(rxn_name)
            for rxn in stoich.parent_reactions:
                # Make sure that ROS only generated, not consumed spontaneously
                if rxn.reverse:
                    if verbosity>0:
                        print 'Skipping reverse rxn:', rxn.id
                else:
                    ros_rxns.append(rxn)

        self.ros_rxns = ros_rxns

        return ros_rxns


    def gross_spont_ros(self, ros_ids = ['h2o2_c', 'o2s_c'], x_dict=None):
        """
        Sum up gross ROS generated spontaneously.
        """
        me = self.me

        if len(self.ros_rxns)==0:
            #warnings.warn('No ROS rxns present!')
            self.get_ros_rxns()

        v_ros_dict = {rid:0. for rid in ros_ids}

        if x_dict is None:
            if me.solution is not None:
                x_dict = me.solution.x_dict

        if x_dict is not None:
            for ros_id in ros_ids:
                v_ros = 0.
                ros = me.metabolites.get_by_id(ros_id)
                v_ros = sum([x_dict[rxn.id] * rxn.metabolites[ros] for rxn in self.ros_rxns])
                v_ros_dict[ros_id] = v_ros
        else:
            warnings.warn('No solution in model & no x_dict provided. Solve first or provide x_dict!')

        return v_ros_dict

    def rand_coeffs_norm(self, params={
                'h2o2_c': {'mean':0.001, 'std':0.00046},
                'o2s_c': {'mean':0.00048, 'std':0.00018}}):
        """
        Randomize ROS coefficients
        """
        ros_coeff_dict = self.ros_coeff_dict
        mult_dict = {}
        nRxns = len(self.ros_rxns)

        for ros_id, par_dict in params.iteritems():
            sym = ros_coeff_dict[ros_id]
            mean = par_dict['mean']
            std  = par_dict['std']
            mults = abs(rand.normal(mean, std, nRxns))
            mult_dict[ros_id] = mults

        return mult_dict


    def scale_coeffs(self, ref_dict={'h2o2_c':0.1233, 'o2s_c':0.044}, x_dict=None,
            fix_coeff_dict={}):
        """
        Scale ROS coefficients so total basal h2o2 and o2s production
        match measured values.
        Randomize or deterministically vary relative coeffs across reactions.
        Scale all by common centered coefficient.

        Brynildsen et al. (2013):
        14uM H2O2/s = 1233x10-4 mmol H2O2 gDW-1 hr-1 using
        a cell volume of 6.8x10-16 L 4 and cell weight of 278x10-15 gDW
        """
        me = self.me
        ros_coeff_dict = self.ros_coeff_dict

        # Gross ROS production expressions
        v_ros_dict = self.gross_spont_ros(x_dict=x_dict)

        # Generate random coeffs
        mult_dict = self.rand_coeffs_norm()
        subs_dict = {}
        scaled_subs_dict = {}

        # Pre-build subs_dict for all ros since c mixed in with k
        # in some rxns due to o2s being reactant for h2o2 production
        for ros_id, mults in mult_dict.iteritems():
            sym = ros_coeff_dict[ros_id]    # i.e., c for h2o2, k for o2s
            mid = Symbol(sym+'m', positive=True)
            for i,si in enumerate(mults):
                coeff_sym = sym+'_'+str(i)
                subs_dict[coeff_sym] = mid*si

        scaled_mid = {}

        for ros_id, v_ref in ref_dict.iteritems():
            ros = me.metabolites.get_by_id(ros_id)
            # Gross v_ros
            v_ros = v_ros_dict[ros_id]
            # Collapse coeffs into scaling * mid 
            sym = ros_coeff_dict[ros_id]    # i.e., c for h2o2, k for o2s
            mid = Symbol(sym+'m', positive=True)

            # Need to use the actual symbol objects in the expr 
            # or else sympy sometimes fails to substitute...
            sym_dict = {s:subs_dict[str(s)] for s in v_ros.free_symbols}
            # Solve: s * cm = v_ref
            try:
                smid = sympy.solve(v_ros.subs(sym_dict) - v_ref)[0]
                scaled_mid[ros_id] = smid
            except IndexError:
                print 'Error with:', 'sympy.solve(v_ros.subs(sym_dict) - v_ref)[0]'
                print 'for v_ros:', v_ros
                print 'for v_ref:', v_ref
                print 'for sym_dict:', sym_dict
                print sympy.solve(v_ros.subs(sym_dict) - v_ref)

            # Substitute in the solved midpoint into subs dict
            for sym,stoich in subs_dict.iteritems():
                subs_dict[sym] = stoich.subs(mid,smid)

        # Retrieve all the actual symbols from ros stoichs since sympy is picky during subs
        for rxn in self.ros_rxns:
            for ros_id in ros_coeff_dict.keys():
                met = me.metabolites.get_by_id(ros_id)
                stoich = rxn._metabolites[met]
                if hasattr(stoich, 'free_symbols'):
                    syms = stoich.free_symbols
                    for sym in syms:
                        scaled_subs_dict[sym] = subs_dict[str(sym)]

        #----------------------------------------------------
        # If user fixed any coeffs, overwrite random ones here
        for stoich_id, coeff_dict in fix_coeff_dict.iteritems():
            stoich = me.stoichiometric_data.get_by_id(stoich_id)
            for ros_id in ros_coeff_dict.keys():
                sym = self.stoich_coeff_dict[stoich.id][ros_id]
                scaled_subs_dict[sym] = coeff_dict[ros_id]
        #----------------------------------------------------

        # Substitute the solved value back in
        for rxn in self.ros_rxns:
            for met in rxn._metabolites:
                if hasattr(rxn._metabolites[met], 'subs'):
                    rxn._metabolites[met] = rxn._metabolites[met].subs(scaled_subs_dict)

        return subs_dict, scaled_mid, scaled_subs_dict
