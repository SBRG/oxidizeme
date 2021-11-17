#============================================================
# File stresstools.py
#
# Methods for creating, simulating, and analyzing StressME
# model and data
#
# Laurence Yang, SBRG, UCSD
#
# 02 Sep 2015:  first version. Port methods from ipynb
# 05 Sep 2015:  added the (better) L1-norm that computes mass
#               fraction exactly as part of L1 problem
# 24 Sep 2015:  changed to stresstools
#============================================================

import numpy as np
from cobra.core.Solution import Solution
from cobra import DictList
from qminospy.me1 import ME_NLP1
import pandas as pd
import copy as cp
import oxidizeme
import time
import warnings
import cobra
import cobrame
import os


def pathto(filename):
    return os.path.join(oxidizeme.__path__[0], 'data', filename)


class stressTools:
    """
    Composite class of ME_NLP containing StressME methods and properties
    Note: ME_NLP is itself a composite of me
    init:
    StressME(me)
    """
    def __init__(self, me):
        self.me = me
        self.me_nlp = ME_NLP1(me)
        self.load_geneid_map()

    def __getattr__(self, attr):
        # Any attribute not in StressME is directed to ME_NLP
        return getattr(self.me_nlp, attr)

    def locus2name(self, locus):
        """Return gene name given locus number"""
        idmap = self.geneid_map
        name = idmap.name[ idmap.locus == locus ].tolist()
        return name

    def name2locus(self, name):
        """Return locus number given gne name"""
        idmap = self.geneid_map
        locus = idmap.locus[ idmap.name == name ].tolist()
        return locus

    def load_geneid_map(self):
        """Load locus - gene name mapping"""
        self.geneid_map = pd.read_csv(pathto('gene_id_map.csv'))

    def pq_to_vo2s(self, pq, v_o2s0=0.26, vmax=60, cap_max=False):
        """
        Calculate o2s flux from PQ concentration

        v_o2s = pq_to_vo2s(pq, vo2s0=0.26)

        Arguments:
        pq : Paraquat concentration (M)
        v_o2s0 : Basal o2s flux from
                 5 uM/s
                 V_cell = 4e-15 L (Volkmer et al., 2011)
                 gDW_cell = 2.78e-13 g
                 (A bit higher than Brynildsen's because cell volume larger)

        Returns:
        v_o2s : mmol/gDW/h

        PQ to fold superoxide dismutase derived in pq_o2s_calibration.ipynb
        """
        slope = 11667
        intercept = 1
        fold_sod = pq*slope + intercept
        v_o2s = v_o2s0 * fold_sod

        if cap_max:
            v_o2s = min(v_o2s, vmax)

        return v_o2s


    def pq_e_to_o2s_c(self, pq, v_dmg=0, kcat_Km_sod=1e9, cell_v=4e-15, cell_gDW=2.78e-13):
        """
        o2s = pq_e_to_o2s_c(pq)

        External PQ to intracellular o2s concentration

        Inputs:
            pq : in M

        Returns
            o2s : in M

        vo2s = vdmg + vSOD
        vSOD = kcat/Km [o2s]
        o2s  = vSOD  / (kcat/Km)

        kcat/Km [1/M/s]
        """
        v_o2s = self.pq_to_vo2s(pq) # mmol/gDW/h
        v_sod = v_o2s - v_dmg
        keff = float(kcat_Km_sod)*1e-3*3600
        o2s_gDW  = v_sod / keff     # mmol/gDW
        o2s_mM  = o2s_gDW * cell_gDW / 4e-15  # mM
        o2s  = o2s_mM * 1e-3

        return o2s


    def o2s_c_to_pq_e(self, o2s):
        """
        pq = o2s_c_to_pq_e(o2s)

        Intracellular o2s concentration to extracellular PQ
        """


    def make_L1_prob(self, df):
        """
        Make L1-norm minimization problem to fit proteome to measured 
        mass fractions in provided dataframe.
        Useful as a first step in StressME reconstruction as it 
        identifies blocked proteins and reactions that would prevent
        simulation of certain stress responses
        """
        me0 = self.me
        # Make copy of me or not...?
        me = cp.deepcopy(me0)
        # Reset obj coeffs
        for rxn in me.reactions:
            rxn.objective_coefficient = 0.
        # Add variable and constrain it to total proteome mass
        cons_p = cobra.Metabolite('p_eq_proteome_mass')

        p = cobra.Reaction('proteome_mass')
        p.lower_bound = 0.
        p.upper_bound = 1000.
        p.objective_coefficient = 0.
        p.add_metabolites({cons_p: 1.})
        me.add_reactions([p])

        for rxn in me.reactions.query('translation'):
            rxn.add_metabolites({cons_p: -rxn.translation_data.mass})

        # Add constraint & aux var for each proteine mf to fit
        for ind,row in dfi.iterrows():
            # Add constraints
            cons_l1 = cobra.Metabolite('error_'+row.locus_tag+'_1')
            cons_l2 = cobra.Metabolite('error_'+row.locus_tag+'_2')

            # Add auxiliary var
            s = cobra.Reaction('s_'+row.locus_tag)
            s.lower_bound = 0.
            s.upper_bound = 1000.
            s.add_metabolites({cons_l1: -1.,
                               cons_l2: -1.})
            s.objective_coefficient = 1.

            # Add slack vars
            y1 = cobra.Reaction('y1_'+row.locus_tag)
            y1.lower_bound = 0.
            y1.upper_bound = 1000.
            y1.objective_coefficient = 0.
            y2 = cobra.Reaction('y2_'+row.locus_tag)
            y2.lower_bound = 0.
            y2.upper_bound = 1000.
            y2.objective_coefficient = 0.
            y1.add_metabolites({cons_l1: 1.})
            y2.add_metabolites({cons_l2: 1.})

            # Translation of i.
            rxn_trsli = me.reactions.get_by_id('translation_'+row.locus_tag)
            mfi = row.mf
            mwi = row.mw
            rxn_trsli.add_metabolites({cons_l1: mwi})
            rxn_trsli.add_metabolites({cons_l2: -mwi})

            # All proteins: can just use the p variable we made earlier
            p.add_metabolites({cons_l1: -mfi,
                               cons_l2:  mfi})

            # Add remaining new columns
            me.add_reactions([s, y1, y2])

        me_nlp_l1 = ME_NLP1(me)
        self.me_nlp_l1 = me_nlp_l1
        print('Finished making L1-norm problem: self.me_nlp_l1')

        return me_nlp_l1


    def make_pacme_prob(self, sectors):
        """
        sectors: list of dict with genes and mass fraction,
                 [{'genes': [gene], 'mf': mf}]
        """
        me2 = cp.deepcopy(me)
        for sind,sector in enumerate(sectors):
            # Make the allocation constraint
            cons = cobra.Metabolite('cons_pac_'+str(sind))
            # Add the slack variable
            svar = cobra.Reaction('slack_pac_'+str(sind))
            svar.lower_bound = 0.
            svar.upper_bound = 1000.
            svar.add_metabolites({cons:1.})
            me2.add_reaction(svar)

            # LHS: all proteins in the sector
            for gene in sector['genes']:
                rxn = me2.reactions.get_by_id('translation_'+gene)
                rxn.add_metabolites({cons: -rxn.translation_data.mass})
            # RHS: all proteins in the model
            for rxn in me2.reactions.query('translation'):
                rxn.add_metabolites({cons: sector['mf'] * rxn.translation_data.mass})

        me2_nlp = qme.ME_NLP1(me2)
        me2_nlp.make_matrices()
                
        return me2_nlp



    def make_sectors_ros(self, dfi, growth_cond='D-Glucose'):

        ## Mass for total ME proteome
        m_all = sum(dfi.fpkm * dfi.length)

        def calc_mf(genes):
            df_sect = dfi[ dfi.Gene.isin(genes)]
            mi = sum(df_sect.fpkm * df_sect.length)
            mf = mi / m_all
            return mf
        
        sector_genes = [
            # Total ROS regulon sector
            get_genes_ros(),
            # Individual regulons
            get_genes_ros_regs(),
            # Core sector
            get_genes_core(),
            # Niche sector
            get_genes_niche(growth_cond)
        ]

        sectors_all = [{'genes':genes, 'mf':calc_mf(genes)} for genes in sector_genes]

        return sectors_all


class StressTools(stressTools):
    """
    Class names should be uppercase...
    """
