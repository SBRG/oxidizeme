#============================================================
# File model.py
# 
# Classes and methods for building a OxidizeME model
#
# Laurence Yang, SBRG, UCSD
#
# 21 Feb 2018:  clean up/refactor for publication
#============================================================

from collections import defaultdict
from cobrame import mu
from cobrame import StoichiometricData
from cobrame import ComplexData
from cobrame import GenericFormationReaction, MetabolicReaction
from cobrame import TranslationReaction
from cobrame import SubreactionData, GenericComponent, GenericData
from cobrame import Constraint, Complex
from cobrame.util import building
from qminospy.me2 import get_gene
from qminospy.me1 import ME_NLP1
from sympy import Basic

from six import iteritems

from cobra import Reaction, Metabolite
from sympy import Symbol
import numpy as np
import warnings
import cobrame
import cobra
import time
import re


class StressME(object):
    """
    StressME model
    """
    def __init__(self, me_nlp, observed=False):
        self.me_nlp = me_nlp
        if observed:
            # ome = ObserveME(me_nlp)
            # self.ome = ome
            warnings.warn('Unsupported: observed=True')
            self.ome = None
            me = me_nlp.me
        else:
            self.ome = None
            me = me_nlp.me

        self.me = me
        self.ros_rxns = []
        self.symbol_dict = {
                'h2o2_c':Symbol('h2o2_c'),
                'o2s_c':Symbol('o2s_c'),
                'h2o2_e':Symbol('h2o2_e'),
                'o2s_e':Symbol('o2s_e'),
                'A':Symbol('A'),
                'P_h2o2_e':Symbol('P_h2o2_e'),
                'P_o2s_e':Symbol('P_o2s_e'),
                'fe2_c':Symbol('fe2_c'),
                'mn2_c':Symbol('mn2_c'),
                'zn2_c':Symbol('zn2_c'),
                'B_fe2_c':Symbol('B_fe2_c'),
                'B_mn2_c':Symbol('B_mn2_c'),
                'B_zn2_c':Symbol('B_zn2_c'),
                'pq2_c':Symbol('pq2_c')}
        self.kcat_Km_dict={'4fe4s':{'h2o2':1e3, 'o2s':1e6},
                'fe2':{'h2o2':1e3, 'o2s':1e6}}

        self.mismet_fe2_dict = {}
        self.builder = building


    def __getattr__(self, attr):
        """Composite of ME"""
        return getattr(self.me, attr)

    def make_stressme(self, force_damage=True, extra_dilution=True,
            lower_repair=True):#, force_o2s_from_pq=True):
        """
        Add base reactions, damage, repair for stressme
        """
        me = self.me
        # Designed to work with me.update()
        self.add_alt_metal_peptide_deformylase()
        # Below may not yet work with me.update()--changes might get cleared 
        self.add_ros_scavenging()
        self.add_btn_oxidation()
        self.add_pq_reduction()
        self.add_alt_metallation()
        self.add_demetallation()
        self.add_FeS_repair_complexes()

        #----------------------------------------------------
        # Dependent
        self.add_FeS_damage_repair_and_Isc()
        #----------------------------------------------------
        self.add_stress_transporters()
        self.add_dps()
        self.add_TF_binding()

        #----------------------------------------------------
        # Mismetallation of IscU, SufA, and possibly IscA and SufU
        self.add_demetallation_FeS_assembly()

        #----------------------------------------------------
        # Need to add extra_dilution
        if extra_dilution:
            _ = self.me_nlp.make_dilution_fluxes(UB_DIL=1000)

        #----------------------------------------------------
        # Need to force damage
        if force_damage:
            self.force_fes_damage(csense='G')
            self.force_demetallation(csense='G')
            self.force_mismetallation(csense='E')

        #----------------------------------------------------
        # Lower repair rates
        if lower_repair:
            keff_fesrep = 0.004 * 3600
            rxns_fes_rep = me.reactions.query('repair_')
            met_fes_rep = me.metabolites.generic_fes_repair

            for rxn in rxns_fes_rep:
                s = rxn._metabolites[met_fes_rep]
                s2 = -mu/keff_fesrep
                rxn._metabolites[met_fes_rep] = s2

        #----------------------------------------------------
        # Allow h2o2 to cross from periplasm to cytosol
        self.add_h2o2_influx()

        #----------------------------------------------------
        # Prevent reverse MOX
        try:
            rxn = me.reactions.get_by_id('MOX_REV_CPLX_dummy')
            rxn.lower_bound = 0.
            rxn.upper_bound = 0.
        except KeyError:
            warnings.warn('MOX_REV_CPLX_dummy not found. Make sure it cannot carry flux')




    def add_FeS_damage_repair_all(self, compt='_c'):
        """
        Add Fe-S cluster enzyme damage and repair.
        Damaged 3fe4s enzymes should also be diluted.
        Especially since 90% of enzyme observed to be in 3fe4s form
        under high oxidative stress.
        """
        me = self.me
        mod = me.process_data.get_by_id('mod_4fe4s'+compt)
        cpxs = list(set([d.complex for d in mod.get_complex_data()]))
        self.add_FeS_enzyme_damage_repair(cpxs)

    # TODO: shouldn't need to do this separately
    def add_IscU_damage_repair(self):
        """
        Damage IscU
        """
        cplx = self.me.metabolites.get_by_id('IscU_mod_1:4fe4s')
        self.add_FeS_enzyme_damage_repair([cplx])


    def add_FeS_repair_complexes(self,
            complex_dict={
                'b2962':{'stoich':1, 'name':'yggX'},
                'b4209':{'stoich':1, 'name':'ytfE'}
                },
            generic_fes_repair_id='generic_fes_repair'):
        """
        Add Fe-S repair enzymes.
        yggX (b2962) and ytfE (b4209): both are up-regulated under pq treatment (SWS)
        """
        me = self.me

        for bnum, props in list(complex_dict.items()):
            cpx_name = props['name']
            stoich = props['stoich']
            try:
                data = ComplexData(cpx_name, me)
                data.stoichiometry = {'protein_'+bnum: props['stoich']}
                data.translocation = defaultdict(int)
#                data.chaperones = defaultdict(int)
                data.create_complex_formation()
            except ValueError:
                data = me.process_data.get_by_id(cpx_name)

            # Complex to generic_fes_repair
            try:
                generic_fes_repair = me.metabolites.get_by_id(generic_fes_repair_id)
            except KeyError:
                generic_fes_repair = cobrame.GenericComponent(generic_fes_repair_id)
                me.add_metabolites([generic_fes_repair])
                #self.modeler.add_metabolites([generic_fes_repair])

            try:
                rxn_name = cpx_name+'_to_'+generic_fes_repair_id
                rxn = GenericFormationReaction(rxn_name)
                me.add_reaction(rxn)

                stoich = {data.complex: -1,
                        generic_fes_repair: 1}
                rxn.add_metabolites(stoich, combine=False)
            except ValueError:
                rxn = me.reactions.get_by_id(rxn_name)


    def add_FeS_enzyme_damage_repair(self, complexes=None, repair_enzs=[], compt='_c',
            e_donor_repair='nadh_c', generic_fes_repair_id = 'generic_fes_repair', verbosity=0):
        """
        DAMAGE Fe-S enzymes by h2o2 and o2s.
        - Mechanism (Imlay ,2013):
          IscU_mod_1:4fe4s(2+) + h2o2 + 2 h(1+) --> IscU_mod_1:3fe4s(1+) + fe3(3+) + 2 h2o
          IscU_mod_1:4fe4s(2+) + o2s + 2 h(1+) --> IscU_mod_1:3fe4s(1+) + fe2(2+) + h2o2
        - 3fe4s can decompose further under prolonged stress (Djaman et al., 2004)

        REPAIR
        - Repair process can be catalyzed by any of the hypothesized Fe-S repair enzymes.
          IscU_mod_1:3fe4s(1+) + fe2(2+) + mu/keff yggX/ytfE/... --> IscU_mod_1:4fe4s(2+) + h(1+)
          Requires electron donor. Assume NADH for now (Brynildsen et al., 2013)

        - The damaged complex itself also needs to be diluted (not just yggX, etc.)
        - Especially considering at high ROS can have 90% inactive complex.

        Also add damaged complex to me.complex_data
        """
        me = self.me

        for cpx_4fe4s in complexes:
            cpx_3fe4s_id = cpx_4fe4s.id.replace('4fe4s','3fe4s')
            try:
                cpx_3fe4s = me.metabolites.get_by_id(cpx_3fe4s_id)
            except KeyError:
                cpx_3fe4s = cobrame.Complex(cpx_3fe4s_id)
                me.add_metabolites(cpx_3fe4s)

            try:
                data_dmg = ComplexData(cpx_3fe4s.id, me)
            except ValueError as e:
                if verbosity>0:
                    print((repr(e)))

            ### h2o2 damage
            rxn_id = 'damage_'+cpx_4fe4s.id+'_h2o2'
            try:
                rxn_damage = me.reactions.get_by_id(rxn_id)
            except KeyError:
                rxn_damage = Reaction(rxn_id)
                me.add_reaction(rxn_damage)

            h2o2 = me.metabolites.get_by_id('h2o2'+compt)
            h = me.metabolites.get_by_id('h'+compt)
            fe3 = me.metabolites.get_by_id('fe3'+compt)
            h2o = me.metabolites.get_by_id('h2o'+compt)
            stoich = {cpx_4fe4s: -1.,
                      h2o2: -1.,
                      h: -2.,
                      fe3: 1.,
                      h2o: 2.,
                      cpx_3fe4s: 1.}
            rxn_damage.add_metabolites(stoich, combine=False)

            ### o2s damage
            rxn_id = 'damage_'+cpx_4fe4s.id+'_o2s'
            try:
                rxn_damage = me.reactions.get_by_id(rxn_id)
            except KeyError:
                rxn_damage = Reaction(rxn_id)
                me.add_reaction(rxn_damage)

            o2s = me.metabolites.get_by_id('o2s'+compt)
            h = me.metabolites.get_by_id('h'+compt)
            fe2 = me.metabolites.get_by_id('fe2'+compt)
            stoich = {cpx_4fe4s: -1.,
                      o2s: -1.,
                      h: -2.,
                      fe2: 1.,
                      h2o2: 1.,
                      cpx_3fe4s: 1.}
            rxn_damage.add_metabolites(stoich, combine=False)

            ### Repair
            try:
                generic_fes_repair = me.metabolites.get_by_id(generic_fes_repair_id)
            except KeyError:
                generic_fes_repair = cobrame.GenericComponent(generic_fes_repair_id)
                me.add_metabolites([generic_fes_repair])

            rxn_id = 'repair_'+cpx_3fe4s_id
            try:
                rxn_repair = me.reactions.get_by_id(rxn_id)
            except KeyError:
                rxn_repair = Reaction(rxn_id)
                me.add_reaction(rxn_repair)

            repair_enz = generic_fes_repair
            keff_repair = 65    # repairing enzyme     
            keff_rep_3fe4s  = 0.004 # Gardner and Fridovich (1992) J Biol Chem

            e_donor = me.metabolites.get_by_id(e_donor_repair)

            stoich = {cpx_3fe4s: -1 - mu/keff_rep_3fe4s,
                      fe2: -1.,
                      e_donor: -1.,
                      repair_enz: -mu/keff_repair/3600.,
                      cpx_4fe4s: 1.,
                      h: 1.
                      }
            rxn_repair.add_metabolites(stoich, combine=False)

    def add_Isc_damage(self, 
            rxn_isc_4fe4s_id = 'isc_4fe4s_formation_FWD_iscUS_cyaY_pathway_complex',
            generic_fes_repair_id='generic_fes_repair',
            keff_repair=65.,
            compt='_c'):
        """
        DEPRECATED
        -- just use add_FeS_damage_repair for IscU

        Add poisoning of Isc by h2o2
        Hypothesized mechanisms:
        - **IscS and IscA are H2O2-resistant.**
        - Therefore, "...oxidants
          disrupt Isc by oxidizing clusters as they are assembled on or transferred from
          the IscU scaffold." (Jang & Imlay, 2010)

        - Imlay (2013) & Jang & Imlay (2010): h2o2 and o2s attack nascent clusters on IscU
        IscU_mod_1:4fe4s(2+) + h2o2 + 2 h(1+) --> IscU_mod_1:3fe4s(1+) + fe3(3+) + 2 h2o
        IscU_mod_1:4fe4s(2+) + o2s + 2 h(1+) --> IscU_mod_1:3fe4s(1+) + fe2(2+) + h2o2

        #----------------------------------------------------
        # Update for iLE1678
        #  IscU_mod_1:4fe4s is its own complex now--convenient
        #----------------------------------------------------
        Damage, h2o2:


          #--------------------------------------------------
        - h2o2: net rxn summing damage and repair (to be compatible current cobrame formulation): 
          4.0 cys__L_c + 4.0 fe2_c + fadh2_c + mu/keff iscUS_cyaY_pathway_complex +
            mu/keff generic_fes_repair + h2o2 + 2 h + 
            fe2_c --> 
          4.0 ala__L_c + fad_c + 4fe4s_c + 6.0 h_c + fe3 + 2 h2o + h 
          #--------------------------------------------------
        - o2s: net rxn summing damage and repair (to be compatible current cobrame formulation): 
          4.0 cys__L_c + 4.0 fe2_c + fadh2_c + mu/keff iscUS_cyaY_pathway_complex +
            mu/keff generic_fes_repair + o2s + 2 h + 
            fe2_c --> 
          4.0 ala__L_c + fad_c + 4fe4s_c + 6.0 h_c + fe2 + h2o2 + h 
          #--------------------------------------------------
        - Damage, h2o2 (net rxn):
          4.0 cys__L_c + 4.0 fe2_c + fadh2_c + mu/keff iscUS_cyaY_pathway_complex + h2o2 + 2 h --> 
              4.0 ala__L_c + fad_c + 6.0 h_c + iscUS_3fe4s_c + fe3 + 2 h2o
        - Damage, o2s (net rxn):
          4.0 cys__L_c + 4.0 fe2_c + fadh2_c + mu/keff iscUS_cyaY_pathway_complex + o2s + 2 h --> 
              4.0 ala__L_c + fad_c + 6.0 h_c + iscUS_3fe4s_c + fe2 + h2o2
          #--------------------------------------------------
        - Assume iscUS_3fe4s_c gets repaired
          (net rxn: assuming 4fe4s_c released)
          iscUS_3fe4s_c + fe2 + mu/keff yggX/ytfE/... --> iscUS_cyaY_pathway_complex + 4fe4s_c + h
        """
        me = self.me
        rxn_isc = me.reactions.get_by_id(rxn_isc_4fe4s_id)
        generic_fes_repair = me.metabolites.get_by_id(generic_fes_repair_id)

        ### h2o2
        rxn_damrep = Reaction('damage_'+'iscUS'+'_h2o2')
        h2o2 = me.metabolites.get_by_id('h2o2'+compt)
        h = me.metabolites.get_by_id('h'+compt)
        fe2 = me.metabolites.get_by_id('fe2'+compt)
        fe3 = me.metabolites.get_by_id('fe3'+compt)
        h2o = me.metabolites.get_by_id('h2o'+compt)

        stoich_damrep = {
                  generic_fes_repair: -mu/keff_repair/3600.,
                  h2o2: -1.,
                  h: -2.,
                  fe2: -1., # need one more to repair
                  fe3: 1.,  # lose as fe3, so need to reduce it
                  h2o: 2.
                  }
        stoich_normal = {k:v for k,v in list(rxn_isc._metabolites.items())} 
        rxn_damrep.add_metabolites(stoich_normal)
        rxn_damrep.add_metabolites(stoich_damrep)
        me.add_reaction(rxn_damrep)

        ### o2s
        rxn_damrep = Reaction('damage_'+'iscUS'+'_o2s')
        o2s = me.metabolites.get_by_id('o2s'+compt)
        h = me.metabolites.get_by_id('h'+compt)
        fe2 = me.metabolites.get_by_id('fe2'+compt)

        stoich_damrep = {
                  generic_fes_repair: -mu/keff_repair/3600.,
                  h2o2: -1.,
                  h: -2.,
                  h2o2: 1.
                  }
        stoich_normal = {k:v for k,v in list(rxn_isc._metabolites.items())} 
        rxn_damrep.add_metabolites(stoich_normal)
        rxn_damrep.add_metabolites(stoich_damrep)
        me.add_reaction(rxn_damrep)


    def add_alt_metal_peptide_deformylase(self,
            metal_keff_dict = {
                'fe2':1019.5963333345715,
                'ni2':1019.5963333345715*11.5/11.9,
                'cobalt2':1019.5963333345715*4.3/11.9,
                'mn2':1019.5963333345715*2.3/11.9,
                'zn2':1019.5963333345715*0.004/11.9},
            compartment='_c'):
        """
        add_alt_metal_peptide_deformylase()

        Add alternate metallated peptide deformylases (Def).

        sum_j [def_mod:j] = [generic_def]
        sum_j mu*[def_mod:j] = mu*[generic_def]
        sum_j vdil:j = vdil_generic

        v_trsl  = sum_j keff_j * [def_mod:j]
                = sum_j keff_j * vdil_j / mu
        --> (prevent mu in denominator)

        The new constraints:
        1) constraint_trsl_def:       mu * v_trsl_i - sum_j keff_j * vdil_j = 0, i\in Proteins
        2) Def_j:                vform_j + vdil_extra_j - vdil_j = 0

        The new reactions:
        vdil_j:         Def_j -->
        vdil_extra_j:   Def_j -->
        ----------------------------------------------------
        """
        me = self.me

        # Remove the Def_mono_mod_1:fe2-specific requirement for each translation reaction
        base_id = 'Def_mono_mod_1:'
        fe2_id  = base_id+'fe2'
        data_fe2 = me.process_data.get_by_id(fe2_id)
        cplx_fe2 = data_fe2.complex
        for data in me.translation_data:
            try:
                data.subreactions.pop('peptide_deformylase_processing')
            except KeyError:
                pass

            for rxn in data.parent_reactions:
                rxn.subtract_metabolites({cplx_fe2: rxn.metabolites[cplx_fe2]})

        # New constraints for each metal
        cons_id = 'constraint_translation_def'
        try:
            cons = Constraint(cons_id)
            me.add_metabolites(cons)
        except ValueError:
            cons = me.metabolites.get_by_id(cons_id)

        cons._bound = 0
        cons._constraint_sense = 'E'

        for metal,keff in list(metal_keff_dict.items()):
            # Make or retrieve the complex
            cplx_id = base_id+metal
            try:
                data = ComplexData(cplx_id, me)
                data.stoichiometry = data_fe2.stoichiometry.copy()
                data.subreactions = {'mod_'+metal+compartment: 1.}
                data.create_complex_formation()
            except ValueError:
                data = me.process_data.get_by_id(cplx_id)

            cplx = data.complex

            dil_id = 'dilution_def_'+metal
            try:
                rxn_dil = Reaction(dil_id)
                me.add_reaction(rxn_dil)
            except ValueError:
                rxn_dil = me.reactions.get_by_id(dil_id)

            rxn_dil.add_metabolites({cplx: -1})
            rxn_dil.lower_bound = 0.
            rxn_dil.upper_bound = 1000.
            rxn_dil.add_metabolites({cons: -keff*3600})

        for data in me.translation_data:
            for rxn in data.parent_reactions:
                # TODO: use subreaction data so not cleared upon me.update()
                # combine=False: Only set the mu coefficient once per translation reaction
                rxn.add_metabolites({cons:mu}, combine=False)


    def add_alt_metallation(self,
            alt_divalents = ['mn2','zn2','cobalt2'],
            keff_over_fe2={'mn2': 1.7/4.2,
                           'zn2': 0.02/4.2,
                           'cobalt2': 0.6/4.2},
            catalyzed_types=[MetabolicReaction]):
        """
        Damage of mononuclear Fe2 enzymes:
        1) Loss of Fe2 from CPLX_fe2 by h2o2 or o2s
        2) Mismetallation with Zn2 during repair of 1)
        Assume that only Fe2-enzymes attacked.
        Response: Alternate metallation of Fe2 enzymes with Mn2.
        keff_zn2 / keff_fe2 = 0.02 / 4.2 = 0.0048
        keff_mn2 / keff_fe2 = 1.7 / 4.2 = 0.4
        keff_cobalt2 / keff_fe2 = 0.6 / 4.2 = 0.14

        Mismetallation rate is a function of h2o2 and o2s and is
        not a constant fraction of the normal modification rate.
        """
        me = self.me
        # Get all fe2 mod complex formation rxns
        # Specifically, MONONUCLEAR Fe(II) enzymes
        datas = [dat for dat in me.complex_data if \
                'mod_fe2_c' in dat.subreactions and \
                dat.subreactions['mod_fe2_c']==1 and \
                len(list(dat.subreactions.keys()))==1]
        # Model mn2-version of complex
        for alt_divalent in alt_divalents:
            mod_comp = alt_divalent+'_c'
            mod_id = 'mod_'+alt_divalent+'_c'
            try:
                mod = cobrame.SubreactionData(mod_id, me)
                mod.stoichiometry = {mod_comp: -1}
            except ValueError:
                print((mod_id, 'already in model.'))

            for data in datas:
                cpx_id = data.id.replace('fe2',alt_divalent)
                try:
                    cpx_data = me.process_data.get_by_id(cpx_id)
                    for rxn in cpx_data.parent_reactions:
                        if isinstance(rxn, cobrame.MetabolicReaction):
                            keff_fe2 = rxn.keff
                            keff_alt = keff_fe2 * keff_over_fe2[alt_divalent]
                            rxn.keff = keff_alt
                            rxn.update()

                except KeyError:
                    cpx_data = cobrame.ComplexData(cpx_id, me)
                    cpx_data.stoichiometry = data.stoichiometry
                    cpx_data.subreactions = {mod_id: data.subreactions['mod_fe2_c']}
#                    cpx_data.chaperones = data.chaperones
                    if hasattr(data,'translocation'):
                        cpx_data.translocation = data.translocation
                    cpx_data.create_complex_formation()

                    # Add metabolic rxns catalyzed by the alternate metallated enzyme with
                    # adjusted keff
                    for rxn in data.parent_reactions:
                        if isinstance(rxn, cobrame.MetabolicReaction):
                            keff_fe2 = rxn.keff
                            keff_alt = keff_fe2 * keff_over_fe2[alt_divalent]
                            # Associate alternate enzyme with the same metabolic rxn
                            # If rxn already exists, just update to alt keff
                            if rxn.reverse:
                                rxn_dir = 'Reverse'
                            else:
                                rxn_dir = 'Forward'
                            self.builder.add_metabolic_reaction_to_model(
                                me, rxn.stoichiometric_data.id,
                                rxn_dir, complex_id=cpx_id, spontaneous=False,
                                update=True, keff=keff_alt)

    def add_demetallation(self, metals=['fe2'], protein_biomass_id='protein_biomass'):
        """
        (rxns_added, mets_added) = add_demetallation(metals)

        Demetallation of mononuclear (iron) metalloenzymes.

        24 Oct 2016:
        If the demetallated product is synthesized protein,
        with explicit protein_biomass component, need to ensure
        that is produced again when demetallated, since otherwise
        get free protein.
        """
        me = self.me
        # Original metallation:
        # fe2_c + 4.0 protein_b1241 --> ADHE-CPLX_mod_fe2
        # Demetallation of fe2:
        # h2o2: 2 h + h2o2 + ADHE-CPLX_mod_fe2 --> 4.0 protein_b1241 + fe2_c + 2 h2o + x protein_biomass
        # o2s : 4 h + o2s + ADHE-CPLX_mod_fe2  --> 4.0 protein_b1241 + fe3_c + 2 h2o + x protein_biomass
        # Remetallation: assume simple mechanism
        # (e.g., no need for sulphenic reduction first)
        # zn2_c + 4.0 protein_b1241 --> ADHE-CPLX_mod_zn2
        # Thus, we can just use add_alt_metallation for re-/mis-metallation

        # The protein_biomass component
        protein_biomass = None

        if protein_biomass_id is not None:
            protein_biomass = me.metabolites.get_by_id(protein_biomass_id)

        compt = '_c'
        mets_added = []
        rxns_added = []

        for metal in metals:
            mod0 = 'mod_'+metal+compt
            datas = [dat for dat in me.complex_data if \
                    mod0 in dat.subreactions and \
                    dat.subreactions[mod0]==1 and \
                    len(list(dat.subreactions.keys()))==1]

            for data in datas:
                ### Demetallation by h2o2: 
                rxn_de_h2o2 = Reaction('demetallation_h2o2_'+data.id)
                stoich = {}
                if protein_biomass is not None:
                    stoich[protein_biomass] = 0.
                # Reverse of metallation stoich for proteins
                for k,v in list(data.stoichiometry.items()):
                    prot = me.metabolites.get_by_id(k)
                    stoich[prot] = abs(v)
                    # The protein_biomass component
                    if protein_biomass is not None:
                        # Mass in kDa
                        stoich[protein_biomass] += abs(v)*prot.formula_weight/1000.

                stoich[data.complex] = -1.
                # Demetallation
                metal_stoich0 = data.subreactions[mod0]
                met_metal = me.metabolites.get_by_id(metal+compt)
                stoich[met_metal] = abs(metal_stoich0)
                # h2o2
                h2o2 = me.metabolites.get_by_id('h2o2'+compt)
                h = me.metabolites.get_by_id('h'+compt)
                h2o = me.metabolites.get_by_id('h2o'+compt)
                # Currently limiting to 1 fe2, so stoich on h2o2 & o2s the same
                stoich[h2o2] = -1.
                stoich[h] = -2.
                stoich[h2o] = 2.
                rxn_de_h2o2.add_metabolites(stoich)

                ### Demetallation by o2s: 
                rxn_de_o2s = Reaction('demetallation_o2s_'+data.id)
                stoich = {}
                if protein_biomass is not None:
                    stoich[protein_biomass] = 0.
                # Reverse of metallation stoich for proteins
                for k,v in list(data.stoichiometry.items()):
                    met = me.metabolites.get_by_id(k)
                    stoich[met] = abs(v)
                    # The protein_biomass component
                    if protein_biomass is not None:
                        # Mass in kDa
                        stoich[protein_biomass] += abs(v)*prot.formula_weight/1000.

                stoich[data.complex] = -1.
                # Demetallation
                metal_stoich0 = data.subreactions[mod0]
                # o2s 
                o2s = me.metabolites.get_by_id('o2s'+compt)
                h = me.metabolites.get_by_id('h'+compt)
                h2o = me.metabolites.get_by_id('h2o'+compt)
                fe3 = me.metabolites.get_by_id('fe3'+compt)
                # Currently limiting to 1 fe2, so stoich on h2o2 & o2s the same
                stoich[o2s] = -1.
                stoich[h] = -4.
                stoich[h2o] = 2.
                stoich[fe3] = abs(metal_stoich0)
                rxn_de_o2s.add_metabolites(stoich)

                me.add_reactions([rxn_de_h2o2, rxn_de_o2s])
                rxns_added = rxns_added + [rxn_de_h2o2, rxn_de_o2s]

        return (rxns_added, mets_added)

    def add_demetallation_FeS_assembly(self):
        """
        Need to treat these complexes differently since
        modified by MetabolicReaction rxns
        """

    def add_btn_oxidation(self, keff=65.):
        me = self.me
        try:
            mstoich = StoichiometricData('BTNOX', me)
            mstoich.lower_bound = 0.
            mstoich.upper_bound = 1000.
            mstoich._stoichiometry = {
                    me.metabolites.h2o2_c.id: -1.,
                    me.metabolites.h2o_c.id: 1.,
                    me.metabolites.btn_c.id: -1.,
                    me.metabolites.btnso_c.id: 1.
            }
            # Spontaneous reaction
            self.builder.add_metabolic_reaction_to_model(
                me, mstoich.id, 'Forward', complex_id=None, spontaneous=True,  update=True, keff=keff)
        except ValueError as e:
            print((repr(e)))


    def add_pq_reduction(self, cpx_id='FLAVONADPREDUCT-MONOMER_mod_fad',
            keff_pq = 65., compt='_c'):
        """
        Add mechanstic process of paraquat reduction
        1. pq2 + nadph --> nadp + pq1 (via FLAVONADPREDUCT-MONOMER_mod_fad)
        2. pq1 + o2 --> o2s + pq2 (Spontaneous)
        """
        me = self.me
        try:
            pq2 = me.metabolites.get_by_id('pq2'+compt)
        except KeyError:
            pq2 = Metabolite('pq2'+compt)
            pq2.name = 'Paraquat (2+)'
            me.add_metabolites([pq2])

        try:
            pq1 = me.metabolites.get_by_id('pq1'+compt)
        except KeyError:
            pq1 = Metabolite('pq1'+compt)
            pq1.name = 'Paraquat monocation radical (1+)'
            me.add_metabolites([pq1])

        ### Rxn step 1:
        try:
            stoich = StoichiometricData('PQ2RED', me)
            stoich.lower_bound = 0.
            stoich.upper_bound = 1000.
            nadph = me.metabolites.get_by_id('nadph'+compt)
            nadp = me.metabolites.get_by_id('nadp'+compt)
            stoich._stoichiometry = {
                    pq2.id: -1.,
                    nadph.id: -1.,
                    pq1.id: 1.,
                    nadp.id: 1.
                    }
            self.builder.add_metabolic_reaction_to_model(
                me, stoich.id, 'Forward', complex_id=cpx_id, spontaneous=False, update=True, keff=keff_pq)
        except ValueError as e:
            print((repr(e)))

        ### Rxn step 2:
        try:
            stoich = StoichiometricData('PQ1OX', me)
            stoich.lower_bound = 0.
            stoich.upper_bound = 1000.
            o2s = me.metabolites.get_by_id('o2s'+compt)
            o2 = me.metabolites.get_by_id('o2'+compt)
            stoich._stoichiometry = {
                    pq1.id: -1.,
                    o2.id: -1.,
                    o2s.id: 1.,
                    pq2.id: 1.
                    }
            # Spontaneous rxn
            self.builder.add_metabolic_reaction_to_model(
                me, stoich.id, 'Forward', complex_id=None, spontaneous=True,  update=True, keff=65.)
        except ValueError as e:
            print((repr(e)))

    def add_stress_transporters(self, keff=65.):
        """
        Add stress-related transporters that were not part of
        the basal ME model, e.g., nepI.
        """
        me = self.me
        ### nepI: b3662 --> YICM-MONOMER
        # Make complex
        cpx_id = 'YICM-MONOMER'
        try:
            cpx_data = cobrame.ComplexData(cpx_id, me)
            cpx_data.stoichiometry = {'protein_b3662': 1.}
            cpx_data.create_complex_formation()
        except:
            print((cpx_id, 'already in model'))

        # Add metabolic rxns
        # gsn
        try:
            stoich = StoichiometricData('GSN_anti_pp', me)
            stoich.lower_bound = 0.
            stoich.upper_bound = 1000.
            h_c = me.metabolites.get_by_id('h_c')
            h_p = me.metabolites.get_by_id('h_p')
            gsn_c = me.metabolites.get_by_id('gsn_c')
            gsn_p = me.metabolites.get_by_id('gsn_p')
            stoich._stoichiometry = {
                    h_p.id: -1.,
                    gsn_c.id: -1.,
                    h_c.id: 1.,
                    gsn_p.id: 1.
                    }
            self.builder.add_metabolic_reaction_to_model(
                me, stoich.id, 'Forward', complex_id=cpx_id, spontaneous=False,  update=True, keff=keff)
        except ValueError as e:
            print((repr(e)))

        # ins
        try:
            stoich = StoichiometricData('INS_anti_pp', me)
            stoich.lower_bound = 0.
            stoich.upper_bound = 1000.
            h_c = me.metabolites.get_by_id('h_c')
            h_p = me.metabolites.get_by_id('h_p')
            ins_c = me.metabolites.get_by_id('ins_c')
            ins_p = me.metabolites.get_by_id('ins_p')
            stoich._stoichiometry = {
                    h_p.id: -1.,
                    ins_c.id: -1.,
                    h_c.id: 1.,
                    ins_p.id: 1.
                    }
            self.builder.add_metabolic_reaction_to_model(
                me, stoich.id, 'Forward', complex_id=cpx_id, spontaneous=False,  update=True, keff=keff)
        except ValueError as e:
            print((repr(e)))


    def add_dps(self, compt='_c',
            keff_dict={'phase12':0.36, 'phase3':0.58, 'catalase':0.12},
            n_subunits=12):
        """
        Dps complex and processes.
        Assume mineralized core does not release iron. Complex must dilute out.
        Consequently, keff should be very low. Need to sample or estimate keff.
        12 Pz sites. Each performs:
        - Phase 1 [fe2 binding]     : 2 fe2 + PZ --> 2fe2-P_Zp2 + 2 h
        - Phase 2 [fe2 oxidation]   : 2 fe2 + 2fe2-P_Zp2 + h2o2 + h2o --> 2fe3o2oh-P_Zm1 + 3 h 
            (1.8 Fe/subunit/min = 0.36 1/s/dps)
        ----------------------------------------------------
        Net reaction for Fe2 bining and oxidation:
        2 fe2 + h2o2 + h2o --> 5 h

        ----------------------------------------------------
        - Phase 3: ferric core mineralization
        - 2 fe2 + h2o2 + 2 h2o --> 2 fe3ooh_core + 4 h
            (2.9 Fe/subunit/min = 0.58 1/s/dps)
        - max 500 Fe(III) mineralized in vitro with H2O2 as oxidant

        ----------------------------------------------------
        Catalase activity (keff=0.12/s): 0.12 umol H2O2/1 umol Fe(III)24-Dps/s 
        h2o2 --> h2o + 1/2 o2

        ----------------------------------------------------
        TODO: DNA-bound shield--Dps becomes highly abundant
        """
        me = self.me
        # Net rxn at 12 PZ sites:
        try:
            stoich = StoichiometricData('HPREDdps', me)
            stoich.lower_bound = 0.
            stoich.upper_bound = 1000.
            fe2 = me.metabolites.get_by_id('fe2'+compt)
            h2o2 = me.metabolites.get_by_id('h2o2'+compt)
            h2o = me.metabolites.get_by_id('h2o'+compt)
            h = me.metabolites.get_by_id('h'+compt)
            stoich._stoichiometry = {
                    fe2.id: -2. * n_subunits,
                    h2o2.id: -1. *n_subunits,
                    h2o.id: -1. *n_subunits,
                    h.id: 5. * n_subunits}
        except ValueError as e:
            print((repr(e)))

        # Make complex
        cpx_id = 'CPLX0-1521'
        try:
            cpx_data = cobrame.ComplexData(cpx_id, me)
            cpx_data.stoichiometry = {'protein_b0812':n_subunits }
            cpx_data.create_complex_formation()
            self.builder.add_metabolic_reaction_to_model(
                me, stoich.id, 'Forward', complex_id=cpx_id, spontaneous=False, update=True,
                keff=keff_dict['phase12'])
        except ValueError as e:
            print((repr(e)))

        #----------------------------------------------------
        # Iron core mineralization 
        #----------------------------------------------------
        # Phase 3: ferric core mineralization
        # 2 fe2 + h2o2 + 2 h2o --> 2 fe3ooh_core + 4 h
        #   (2.9 Fe/subunit/min = 0.58 1/s/dps)
        try:
            stoich = StoichiometricData('FE3COREMINdps', me)
            stoich.lower_bound = 0.
            stoich.upper_bound = 1000.
            stoich._stoichiometry = {
                    fe2.id: -2. * n_subunits,
                    h2o2.id: -1. *n_subunits,
                    h2o.id: -2. *n_subunits,
                    h.id: 4. * n_subunits}
            # Associate the enzyme
            self.builder.add_metabolic_reaction_to_model(
                me, stoich.id, 'Forward', complex_id=cpx_id, spontaneous=False, update=True,
                keff=keff_dict['phase3'])
        except ValueError as e:
            print((repr(e)))

    def add_ros_scavenging(self, keff=7071):
        """
        AHPRED (AhpCF) is the primary h2o2 scavenger at low h2o2 concentrations.
        Kinetically more efficient than catalase but saturated at low (1E-5 M) of h2o2.
        """
        me = self.me

        # Check that data doesn't already exist
        stoich_id = 'AHPRED'
        cplx_id = 'CPLX0-245'

        if not me.process_data.has_id(stoich_id):
            stoichdata = StoichiometricData(stoich_id, me)
            stoichdata.lower_bound = 0.
            stoichdata.upper_bound = 1000.
            stoichdata._stoichiometry = {
                    me.metabolites.h2o2_c.id: -1.,
                    me.metabolites.nadh_c.id: -1.,
                    me.metabolites.h_c.id: -1.,
                    me.metabolites.nad_c.id: 1.,
                    me.metabolites.h2o_c.id: 2.
            }
            ahpred_id = 'CPLX0-245'
            complex_data = cobrame.ComplexData(ahpred_id, me)
            complex_data.stoichiometry = {'protein_b0605': 10, 'protein_b0606': 4}
            complex_data.create_complex_formation()
            self.builder.add_metabolic_reaction_to_model(
                me, stoichdata.id, 'FORWARD', ahpred_id, False, True, keff) 


    def add_TF_binding(self, TF_dict={
        'OxyR':{'stoich':{'protein_b3961': 1},'keff':65.},
        'SoxS':{'stoich':{'protein_b4062': 1},'keff':65.},
        'SoxR':{'stoich':{'protein_b4063': 1},'keff':65.}
        }):
        """
        Add sinks for TFs so they can be synthesized.
        Should eventually make them do something useful.
        """
        me = self.me

        for tfi, values in list(TF_dict.items()):
            cpx_stoich = values['stoich']
            keff_tf = values['keff']
            tf_id = 'TF_'+tfi
            try:
                cpx_data = cobrame.ComplexData(tf_id, me)
                cpx_data.stoichiometry = cpx_stoich
                cpx_data.create_complex_formation()

                # Create dummy utility processes for TFs
                target = cobrame.GenericComponent('targets_'+tfi)
                bound  = cobrame.GenericComponent('bound_'+tfi)
                me.add_metabolites([target, bound])

                trn_id = 'Regulation_by_'+tfi
                trn_stoich = StoichiometricData(trn_id, me)
                trn_stoich.lower_bound = 0.
                trn_stoich.upper_bound = 1000.
                trn_stoich._stoichiometry = {
                        target.id: -1.,
                        bound.id: 1.  }
                tf_cycle = cobra.Reaction('unbind_'+tfi)
                tf_cycle.add_metabolites({bound: -1., target: 1.})
                tf_cycle.lower_bound = 0.
                tf_cycle.upper_bound = 1000.

                # Rxn for TF-target binding
                self.builder.add_metabolic_reaction_to_model(
                    me, trn_stoich.id, 'Forward', complex_id=tf_id, spontaneous=False,
                    update=True, keff=keff_tf)

                # Rxn allowing unbinding
                me.add_reaction(tf_cycle)

            except ValueError as e:
                print((repr(e)))

    def make_oxidizeme(self, force_damage=True, extra_dilution=True,
            lower_repair=True):#, force_o2s_from_pq=True):
        """
        Add base reactions, damage, repair for oxidizeme
        """
        me = self.me
        # Designed to work with me.update()
        self.add_alt_metal_peptide_deformylase()
        # Below may not yet work with me.update()--changes might get cleared 
        self.add_ros_scavenging()
        self.add_btn_oxidation()
        self.add_pq_reduction()
        self.add_alt_metallation()
        self.add_demetallation()
        self.add_FeS_repair_complexes()

        #----------------------------------------------------
        # Dependent
        self.add_FeS_damage_repair_and_Isc()
        #----------------------------------------------------
        self.add_stress_transporters()
        self.add_dps()
        self.add_TF_binding()

        #----------------------------------------------------
        # Mismetallation of IscU, SufA, and possibly IscA and SufU
        self.add_demetallation_FeS_assembly()

        #----------------------------------------------------
        # Need to add extra_dilution
        if extra_dilution:
            _ = self.me_nlp.make_dilution_fluxes(UB_DIL=1000)

        #----------------------------------------------------
        # Need to force damage
        if force_damage:
            self.force_fes_damage(csense='G')
            self.force_demetallation(csense='G')
            self.force_mismetallation(csense='E')

        #----------------------------------------------------
        # Lower repair rates
        if lower_repair:
            keff_fesrep = 0.004 * 3600
            rxns_fes_rep = me.reactions.query('repair_')
            met_fes_rep = me.metabolites.generic_fes_repair

            for rxn in rxns_fes_rep:
                s = rxn._metabolites[met_fes_rep]
                s2 = -mu/keff_fesrep
                rxn._metabolites[met_fes_rep] = s2

        #----------------------------------------------------
        # Allow h2o2 to cross from periplasm to cytosol
        self.add_h2o2_influx()

        #----------------------------------------------------
        # Prevent reverse MOX
        try:
            rxn = me.reactions.get_by_id('MOX_REV_CPLX_dummy')
            rxn.lower_bound = 0.
            rxn.upper_bound = 0.
        except KeyError:
            warnings.warn('MOX_REV_CPLX_dummy not found. Make sure it cannot carry flux')


    def add_FeS_damage_repair_and_Isc(self):
        self.add_FeS_damage_repair_all()
        self.add_IscU_damage_repair()


    def make_ensemble(self, n_ens=100):
        """
        Make parameter ensemble of ROS damage models.
        May not apply as broadly in damage rxns...
        """
        me = self.me


    def substitute_ros(self, me_nlp,
            subs_dict={'h2o2_c':5e-8, 'o2s_c':2e-10,
                'h2o2_e':1.5e-6, 'o2s_e':0, 'pq2_c':0.}):
        """
        substitute_ros(self, me_nlp, subs_dict={'h2o2_c':5e-8, 'o2s_c':2e-10,
                'h2o2_e':1.5e-6, 'o2s_e':0, 'pq2_c':0.}):

        Substitute ROS concentrations into stoichiometries.
        50 nM h2o2 and 0.2 nM o2s are basal levels.
        Extracellular h2o2 and o2s also specified to explore various stress scenarios.
        OxyR activated at 200 nM h2o2.
        Growth defects at 400 nM h2o2.
        """
        for ros,v in list(subs_dict.items()):
            sym = self.symbol_dict[ros]
            me_nlp.substitution_dict[sym] = v


    def substitute_metal(self, me_nlp,
            subs_dict={'fe2_c':1e-5, 'zn2_c':1e-4, 'mn2_c':1e-5,
                'B_mn2_c':1e+10, 'B_fe2_c':1e+11, 'B_zn2_c':1e+12}):
        """
        substitute_metal(me_nlp, subs_dict)

        Set free (unincorporated) metal ion concentrations.
        """
        for met,v in list(subs_dict.items()):
            sym = self.symbol_dict[met]
            me_nlp.substitution_dict[sym] = v


    def scale_damage(self, rxn_meas_dict):
        """
        scale_damage(rxn_meas_dict)

        rxn_meas_dict:  {rxn: measured flux} dict

        Scale damage coefficients to match measurements.
        """

    def scale_damage_prob(self, df_prob, col_prob='activity',
            multimer_mode='max', bool_threshold=0.1):
        """
        scale_damage_prob(self, df_prob)

        Scale damage rates using the ssbio+Bayesian net
        Inputs
        df_prob : data frame of gene (locus) and P(activity decrease by ROS)
        multimer_mode : 'max': if multimer, use max prob. 'average': use avg prob. default='max'
        """
        solver = self.me_nlp
        me = solver.me

        sym_h2o2 = self.symbol_dict['h2o2_c']
        sym_o2s  = self.symbol_dict['o2s_c']
        ros_syms = [sym_h2o2, sym_o2s]

        dmg_cplxs = me.metabolites.query('force_(demet|damag)')
        #print '%d damaged complex constraints' % len(dmg_cplxs)
        pat = re.compile(r"force_(damage|demetallation)_(\S+)_[h2o2|o2s]")
        pat2 = re.compile(r"(\S+)_mod_\S+")
        pat_locus = re.compile(r"protein_(b\d+)(_\S)?")

        self.pert_dict = {}

        for imet,met in enumerate(dmg_cplxs):
            # Get the subunit locuses
            data = None
            did  = re.findall(pat,met.id)[0][1]
            if me.process_data.has_id(did):
                data = me.process_data.get_by_id(did)
            else:
                did2 = re.findall(pat2,did)[0]
                try:
                    data = me.process_data.get_by_id(did2)
                except KeyError:
                    print(('No data for %s nor for %s' % (did, did2)))

            if data is not None:
                loci = [re.findall(pat_locus,k)[0][0] for k in list(data.stoichiometry.keys())]
                # Get the probability of damage and deactivation
                dfi = df_prob[ df_prob.locus.isin(loci)]
                if multimer_mode=='max':
                    p_deact = dfi[col_prob].max()
                elif multimer_mode=='average':
                    p_deact = dfi[col_prob].mean()
                else:
                    warnings.warn("Defaulting to multimer_mode='average'")
                    p_deact = dfi[col_prob].mean()

                # Scale activity by deactivation probability
                self.pert_dict[met.id] = p_deact
                # Damage rate constant --> 0 as p_deact --> 0

                # If sharp threshold applied:
                if bool_threshold is not None:
                    if p_deact < bool_threshold:
                        p_deact = 0. #sharp_threshold
                    else:
                        p_deact = 1.


                if not np.isnan(p_deact):
                    for rxn in met.reactions:
                        s = rxn.metabolites[met]
                        if hasattr(s,'free_symbols') and any(
                                [sym in s.free_symbols for sym in ros_syms]):
                            s_new = rxn.metabolites[met]*p_deact
                            rxn._metabolites[met] = s_new
                            # Update compiled expression
                            update_stoich_expr(solver, met, rxn)

        return self.pert_dict

    def add_h2o2_influx(self, rxn_id='diffusion_h2o2_pc', h2o2_p='h2o2_p', h2o2_c='h2o2_c'):
        """
        Allow h2o2 to cross from periplasm to cytosol.
        """
        me = self.me

        h2o2_p = me.metabolites.get_by_id(h2o2_p)
        h2o2_c = me.metabolites.get_by_id(h2o2_c)
        try:
            rxn = me.reactions.get_by_id(rxn_id)
        except KeyError:
            rxn = Reaction(rxn_id)
            rxn.add_metabolites({h2o2_p: -1., h2o2_c: 1.})
            me.add_reaction(rxn)
        rxn.lower_bound = 0.
        rxn.upper_bound = 1000.


    def calc_ros_influx(self, Cout, Cin, P=1.6e-3, A=1.41e-7, gDW=278e-15):
        """
        calc_ros_influx(Cout, Cin, P,A,LV,gDW)

        ROS influx: default params for h2o2
        v = (Cout-Cin)*P*A*1e-3 (L/cm^3)
        """
        LV=1e-3                 # L/cm3 
        J = (Cout-Cin)*P*A*LV   # M/cell/s
        v = J/gDW * 1e3 * 3600  # mmol/gDW/h

        return v


    def force_fes_damage(self, ros_compartment='_c', damaged_complexes=None, csense='G'):
        """
        force_fes_damage()

        kcat_Km_dict: {ros: kcat/Km in 1/M 1/s}

        Add constraints to force ROS damage to Fe-S clusters.
        Can scale individual kcat/Km of damage, creating ensemble.

        v_dmg = kcat/Km * [ROS] * v_dil / mu

        v_dmg - kcat/Km*ROS/mu * v_dil = 0
        v_dil = sum_j mu/keff*v_use + v_extra_dil
        """
        me = self.me
        kcat_Km_dict = self.kcat_Km_dict['4fe4s']

        if damaged_complexes is None:
            damaged_complexes = [
                    m for rxn in me.reactions.query('damage_') for
                    m in list(rxn.metabolites.keys()) if rxn.metabolites[m] < 0 and
                    isinstance(m, Complex)]

        for ros,kcat_Km in list(kcat_Km_dict.items()):
            sym = self.symbol_dict[ros+ros_compartment]
            keff = 3600*kcat_Km*sym
            for cplx in damaged_complexes:
                cons = Constraint('force_damage_'+cplx.id+'_'+ros)
                cons._constraint_sense = csense
                cons._bound = 0
                me.add_metabolites(cons)

                for rxn in cplx.reactions:
                    if 'damage' in rxn.id and ros in rxn.id:
                        rxn.add_metabolites({cons:1}, combine=False)
                dil_dict = self.get_dilution_dict(cplx)
                # Dilution of this complex coupled to h2o2 AND o2s
                # but each dmg flux has its own constraint
                for vdil,s in list(dil_dict.items()):
                    vdil.add_metabolites({cons: -keff/mu*s}, combine=False)
                    self.ros_rxns.append(vdil)


    def force_demetallation(self, ros_compartment='_c', damaged_complexes=None, csense='G'):
        """
        force_demetallation()

        kcat_Km_dict: {ros: kcat/Km in 1/M 1/s}

        Add constraints to force demetallation.
        Can scale individual kcat/Km of damage, creating ensemble.

        v_dmg = kcat/Km * [ROS] * v_dil / mu

        v_dmg - kcat/Km*ROS/mu * v_dil = 0
        v_dil = sum_j mu/keff*v_use + v_extra_dil
        """
        me = self.me

        kcat_Km_dict = self.kcat_Km_dict['fe2']

        if damaged_complexes is None:
            damaged_complexes = [
                    m for rxn in me.reactions.query('demetallation_') for
                    m in list(rxn.metabolites.keys()) if isinstance(m, Complex) and
                    rxn.metabolites[m] < 0
                    ]

        for ros,kcat_Km in list(kcat_Km_dict.items()):
            sym = self.symbol_dict[ros+ros_compartment]
            keff = 3600*kcat_Km*sym
            for cplx in damaged_complexes:
                cons_id = 'force_demetallation_'+cplx.id+'_'+ros
                if me.metabolites.has_id(cons_id):
                    cons = me.metabolites.get_by_id(cons_id)
                else:
                    cons = Constraint(cons_id)
                cons._constraint_sense = csense
                cons._bound = 0
                for rxn in cplx.reactions:
                    if 'demetallation_'+ros in rxn.id:
                        rxn.add_metabolites({cons:1}, combine=False)
                dil_dict = self.get_dilution_dict(cplx)
                # Dilution of this complex coupled to h2o2 AND o2s
                # but each dmg flux has its own constraint
                for vdil,s in list(dil_dict.items()):
                    vdil.add_metabolites({cons: -keff/mu*s}, combine=False)
                    self.ros_rxns.append(vdil)

    def force_mismetallation(self, W=278e-15/6.8e-16,
            metal_dict={
                'zn2':{'k_demet':0.002},
                'mn2':{'k_demet':0.002}},
            complex_metal_dict=None,
            metal_compartment='_c',
            csense = 'E'
            ):
        """
        force_mismetallation()

        Upon demetallation, Fe(II) metalloenzymes are rapidly re- or mis-metallated
        with relative metallation distribution determined thermodynamically.

        vdem_i = kdem_i * [E:i]
        vdil_i = mu*[E:i]
        vdil_i >= mu*sum_j vuse/keffj
        [E:i]/Keqi/[i] + sum_j [E:j] = vsynth_i/mu * W

        Creates symbols in stoichiometries that are substituted later,
        by substitute_metal.
        """
        me = self.me
        # If complexes not explicitly specified,
        # all demetallated Fe(II) complexes are candidates for
        # mismetallation
        mismet_fe2_dict = {}
        if complex_metal_dict is None:
            complex_metal_dict = {}
            damaged_fe2_complexes = [
                    m for rxn in me.reactions.query('demetallation_') for
                    m in list(rxn.metabolites.keys()) if isinstance(m, Complex) and
                    rxn.metabolites[m] < 0
                    ]
            # Retrieve the mismetallated enzymes
            for metal in list(metal_dict.keys()):
                for fe2_cplx in damaged_fe2_complexes:
                    if 'mod_fe2' in fe2_cplx.id:
                        mod_alt = 'mod_'+metal
                        alt_id = fe2_cplx.id.replace('mod_fe2', mod_alt)
                        alt_cplx = me.metabolites.get_by_id(alt_id)
                        # Append to the mismetallated complex list
                        complex_metal_dict[alt_cplx] = metal
                        # Map mismetallated to fe2 complex
                        mismet_fe2_dict[alt_cplx] = fe2_cplx
                    elif 'mod_1:fe2' in fe2_cplx.id:
                        mod_alt = 'mod_1:'+metal
                        alt_id = fe2_cplx.id.replace('mod_1:fe2', mod_alt)
                        alt_cplx = me.metabolites.get_by_id(alt_id)
                        # Append to the mismetallated complex list
                        complex_metal_dict[alt_cplx] = metal
                        # Map mismetallated to fe2 complex
                        mismet_fe2_dict[alt_cplx] = fe2_cplx

        self.mismet_fe2_dict = mismet_fe2_dict

        print('//****************************************************')
        print((len(metal_dict), 'metals'))
        print((len(complex_metal_dict), 'additional mismetallated enzymes'))
        print('//****************************************************')

        kcat_Km_dict = self.kcat_Km_dict['fe2']

        for cplx,metal in list(complex_metal_dict.items()):
            param_dict = metal_dict[metal]
            metal_sym = self.symbol_dict[metal+metal_compartment]
            #--------------------------------------------
            # Couple mismetallation to E:Fe dilution
            # vmet_zn = B_zn/B_fe*[zn]/[fe]*(sum_rosj (kcat/Km *[ROS])+mu)/mu * vdil_E:fe
            # vmet_zn - B_zn/B_fe*[zn]/[fe]*(sum_rosj (kcat/Km *[ROS])+mu)/mu * vdil_E:fe = 0  (OR)
            # vmet_zn - B_zn/B_fe*[zn]/[fe]*(sum_rosj (kcat/Km *[ROS])+mu)/mu * vdil_E:fe >= 0
            cons_id = 'mismetallation_coupling_'+cplx.id
            try:
                cons = me.metabolites.get_by_id(cons_id)
            except KeyError:
                cons = Constraint(cons_id)
                me.add_metabolites(cons)
            cons._bound = 0.
            cons._constraint_sense = csense

            rxn_met = me.process_data.get_by_id(cplx.id).formation
            rxn_met.add_metabolites({cons: 1.})

            cplx_fe2 = mismet_fe2_dict[cplx]
            B_fe2 = self.symbol_dict['B_fe2_c']
            B_i = self.symbol_dict['B_'+metal+metal_compartment]
            sym_fe2 = self.symbol_dict['fe2'+metal_compartment]
            # Split up the complicated stoich
            # stoich: B_zn/B_fe*[zn2]/[fe2]
            stoich = B_i/B_fe2*metal_sym/sym_fe2
            for ros in list(kcat_Km_dict.keys()):
                demet_id = 'demetallation_'+ros+'_'+cplx_fe2.id
                rxn_demet = me.reactions.get_by_id(demet_id)
                rxn_demet.add_metabolites({cons: -stoich})

            dil_dict = self.get_dilution_dict(cplx_fe2)
            for vdil,s in list(dil_dict.items()):
                # s includes the mu/keff
                vdil.add_metabolites({cons: -stoich*s})

            #--------------------------------------------
            # Allow demetallation of mismetallated non-Fe(II) enzyme
            rid = 'demetallation_cys_'+cplx.id
            try:
                rxn_demet = me.reactions.get_by_id(rid)
            except KeyError:
                # Demetallation: reverse the complex formation
                rxn_mod = me.process_data.get_by_id(cplx.id).formation
                stoich = {k:-v for k,v in list(rxn_mod.metabolites.items()) if v<0}
                stoich[cplx]=-1
                rxn_demet = Reaction(rid)
                rxn_demet.add_metabolites(stoich)
                me.add_reaction(rxn_demet)

            #--------------------------------------------
            # Allow dilution of mismetallated enzymes--allow them to be not used
            dil_id = 'extra_dilution_'+cplx.id
            try:
                rxn_dil = Reaction(dil_id)
                me.add_reaction(rxn_dil)
            except :
                rxn_dil = me.reactions.get_by_id(dil_id)

            rxn_dil.add_metabolites({cplx: -1.}, combine=False)
            rxn_dil.lower_bound = 0.
            rxn_dil.upper_bound = 1000.


    def get_dilution_dict(self, cplx, extra_dil_prefix='extra_dilution_',
            excludes=['damage_','demetallation_'],
            rxn_types=[MetabolicReaction, TranslationReaction]):
        """
        get_dilution_dict
        Get total dilution for this rxn = sum_j vuse + extra_dilution
        """
        me = self.me

        # Just want the coefficient on mu (1/keff). Then, multiply mu back on.
        # I.e., don't want mu/keff + 1, etc. The +1 part does not contribute to dilution.
        # vdil = mu/keff * v
        dil_dict = {r:-r.metabolites[cplx].coeff(mu)*mu for r in cplx.reactions if
                r.metabolites[cplx]<0 and
                hasattr(r.metabolites[cplx],'subs') and
                any([isinstance(r,t) for t in rxn_types]) and
                all([s not in r.id for s in excludes])}
        rid_extra_dil = extra_dil_prefix + cplx.id

        # extra_dilution is just an extra sink for unused protein
        if me.reactions.has_id(rid_extra_dil):
            rxn = me.reactions.get_by_id(rid_extra_dil)
            dil_dict[rxn] = -rxn.metabolites[cplx]

        # Add explicit dilution rxns, too
        for rxn in cplx.reactions:
            if 'dilution_' in rxn.id and rxn.metabolites[cplx]<0:
                dil_dict[rxn] = -rxn.metabolites[cplx]

        return dil_dict


    def force_detox(self, kcat_Km_dict={
        'h2o2_c':{
            'AHPRED_FWD_CPLX0-245':4e7, #1/M 1/s,
            'CAT_FWD_HYDROPEROXIDI-CPLX_mod_pheme':9e5, #4.84e5,
            'CAT_FWD_HYDROPEROXIDII-CPLX_mod_6:hemed':1.3e6},
        'o2s_c':{
            'SPODM_FWD_SUPEROX-DISMUTMN-CPLX_mod_mn2':1e9,
            'SPODM_FWD_SUPEROX-DISMUTFE-CPLX_mod_fe3':1e9,
            #'SPODMpp_FWD_G6886-MONOMER_mod_cu2':1e9 # for o2s_p
            }
        }, csense='G'):
        """
        force_detox(self, kcat_Km_dict)

        kcat_Km_dict: {ros: {rxn: kcat/Km in 1/M 1/s}}

        Given ROS concentrations, force detoxification based on
        known kcat/Km:

        1) [ROS] = sum_j vdmg_j / (kdmg_j*E) + sum_j vdetox_j / (kdetox_j*E)
        2) vgen  = sum_j vdmg_j + sum_j vdetox_j


        sum_j vj >= sum_j kcatj/Kmj * ROS * vdil/mu
                 >= sum_j keff(ros) * vdil/mu
        Thus:
        sum_j vj - sum_j keff(ros)*vdil/mu >= 0
        """
        me = self.me

        for ros,rxn_k in list(kcat_Km_dict.items()):
            cons = Constraint('detox_'+ros)
            cons._constraint_sense = csense
            cons._bound = 0
            sym = self.symbol_dict[ros]
            for rid,kcat_Km in list(rxn_k.items()):
                rxn = me.reactions.get_by_id(rid)
                # Update keff
                keff = 3600*kcat_Km*sym
                rxn.keff = keff
                rxn.update()
                self.ros_rxns.append(rxn)
                # Add the constraint:
                # sum_j vj - sum_j keff(ros)*vdil/mu >= 0
                # For now, assume v_dil = v_formation.
                rxn.add_metabolites({cons: 1.}, combine=True)
                if hasattr(rxn, 'complex_data'):
                    #cplx = rxn.complex_data.complex
                    rxn_dil = rxn.complex_data.formation
                    rxn_dil.add_metabolites({cons: -keff/mu}, combine=True)
                    self.ros_rxns.append(rxn_dil)


    def solve(self, mu_fixed=None, basis=None, subs_ros=None, subs_metal=None):
        """
        muopt = solve(self, mu_fixed=None, basis=None)

        Should have run stress.substitute_ros and stress.substitute_metal first.
        Otherwise, subs_ros and subs_metal must be provided.

        If mu_fixed is None, solves bisection to max mu,
        else, solves LP at mu=mu_fixed
        """
        me = self.me
        me_nlp = self.me_nlp
        if subs_ros:
            self.substitute_ros(me_nlp, subs_dict=subs_ros)

        if subs_metal:
            self.substitute_metal(me_nlp, subs_dict=subs_metal)

        if basis is None:
            basis = me_nlp.lp_hs

        if mu_fixed is None:
            muopt, hs, xopt, cache = me_nlp.bisectmu(verbosity=0, check_feas0=True, basis=basis)
        else:
            x,status,hs = me_nlp.solvelp(mu_fixed, basis=basis)
            muopt = mu_fixed

        return muopt #, hs

    def get_ros_conc(self):
        """
        Back-calculate ROS concentration from:
        [ROS] = [ROS]_detox + [ROS]_damage
              = sum_j\in Detox vj*mu/(keffj*vdilj) + sum_j\in Damage vj*mu/(keffj*vdilj)

        Useful for Simulation mode B: force damage by ROS flux exceeding detoxification
        capacity
        """

    def force_o2s_from_pq(self,
            pq_o2s_ratio=1e4,
            pqred_id='PQ2RED_FWD_FLAVONADPREDUCT-MONOMER_mod_fad'):
        """
        Force o2s from PQ
        """
        from oxidizeme.stresstools import StressTools

        me = self.me
        me_nlp = self.me_nlp
        tools = StressTools(me)
        sym_o2s = self.symbol_dict['o2s_c']
        o2s_conc = self.me_nlp.substitution_dict[sym_o2s]
        pq_conc = o2s_conc * pq_o2s_ratio
        vo2s  = tools.pq_to_vo2s(pq_conc)
        rxn = me.reactions.get_by_id(pqred_id)
        rxn.lower_bound = vo2s


    def couple_metal_conc(self, metal):
        """
        Couple metal concentration to flux
        """
        if metal == 'fe2_c':
            warnings.warn('For fe2_c, using dedicated function, couple_iron_conc()')
            self.couple_iron_conc(metal)
            return

        me = self.me
        solver = self.me_nlp


    def couple_iron_conc(self, kcat_km_fenton=25000,
            fenton_constraint_sense='E',
            k_dmg_dna = 10000.):
        """
        couple_iron_conc(self, kcat_km_fenton=25000,
            fenton_constraint_sense='E',
            k_dmg_dna = 10000.)

        kcat_km_fenton: kcat/Km * [H2O2] * fe2 = v_Fenton
        Iron concentration to flux coupling
        Accounts for Fenton chemistry, DNA damage,
        Dps, total incorporated Fe(II) distribution across proteins
        """
        # x_total = sum_i E_i s_i + x_free
        # v_Fenton = kcat/Km * [H2O2] * x_free
        me = self.me
        solver = self.me_nlp

        h2o2_sym = self.symbol_dict['h2o2_c']
        h2o2_c = solver.substitution_dict[h2o2_sym]
        # Make constraint or update it
        cons_fenton_id = 'cons_fenton'
        if me.metabolites.has_id(cons_fenton_id):
            cons_fenton = me.metabolites.get_by_id(cons_fenton_id)
        else:
            cons_fenton = Constraint(cons_fenton_id)
            me.add_metabolites([cons_fenton])
        cons_fenton._constraint_sense = fenton_constraint_sense

        rxn_fenton_id = 'Fenton_FWD_SPONT'
        if me.reactions.has_id(rxn_fenton_id):
            rxn_fenton = me.reactions.get_by_id(rxn_fenton_id)
        else:
            rxn_fenton = Reaction(rxn_fenton_id)
            rxn_fenton.lower_bound = 0.
            rxn_fenton.upper_bound = 1000.
            me.add_reaction(rxn_fenton)

        # v_Fenton = kcat_km_fenton * h2o2_c * fe2_free
        # fe2_free = fe2_total - sum_i Ei si
        # v_Fenton = kcat_km_fenton * h2o2_c * (fe2_total - sum_i Ei si)
        # THEREFORE:
        # v_Fenton + kcat_km_fenton*h2o2_c*(sum_i Ei si) = kcat_km_fenton*h2o2_c*fe2_total
        # v_Fenton + kcat_km_fenton*h2o2_c*(sum_i vdili/mu si)=kcat_km_fenton*h2o2_c*fe2_tot
        fe2_total = self.symbol_dict['fe2_c']
        cons_fenton._bound = kcat_km_fenton*h2o2_sym*fe2_total

        # Force Fenton rxn
        rxn_fenton.add_metabolites({cons_fenton: 1.})
        # And define the actual Fenton chemistry
        # h2o2 + fe2 --> ho_rad + oh- + fe3
        ### SIMPLIFIED:
        h2o2 = me.metabolites.get_by_id('h2o2_c')
        fe2  = me.metabolites.get_by_id('fe2_c')
        fe3  = me.metabolites.get_by_id('fe3_c')
        dna_dmg_id = 'DNA_damaged'
        if me.metabolites.has_id(dna_dmg_id):
            dna_dmg = me.metabolites.get_by_id(dna_dmg_id)
        else:
            dna_dmg = Metabolite(dna_dmg_id)
            me.add_metabolites([dna_dmg])
        dna  = me.metabolites.get_by_id('DNA_biomass')
        # Create damaged DNA, which can go to DNA_biomass, competing with DNA_biomass_dilution,
        # or be "repaired" by Dps
        rxn_dmg_dna_id = 'DNA_damage'
        if me.reactions.has_id(rxn_dmg_dna_id):
            rxn_dmg_dna = me.reactions.get_by_id(rxn_dmg_dna_id)
        else:
            rxn_dmg_dna = Reaction(rxn_dmg_dna_id)
            rxn_dmg_dna.lower_bound = 0.
            rxn_dmg_dna.upper_bound = 1000.
            me.add_reaction(rxn_dmg_dna)
        # DNA_damaged --> DNA_biomass
        rxn_dmg_dna.add_metabolites({dna_dmg:-1, dna:1.}, combine=False)

        # Create unwanted (damaged DNA) that competes with actual DNA_replication rxn,
        # which also produces DNA_biomass
        # DNA_replication is constrained to be mu.
        # Thus, if DNA_replication rate goes down, then so must max mu
        rxn_fenton.add_metabolites({h2o2:-1, fe2:-1, dna_dmg:k_dmg_dna*3600, fe3:1})

        for data in me.process_data.mod_fe2_c.get_complex_data():
            cplx = data.complex
            dil_dict = self.get_dilution_dict(cplx)
            s_mod = data.subreactions['mod_fe2_c']
            for vdil,s_enz in iteritems(dil_dict):
                vdil.add_metabolites({
                    cons_fenton: s_enz/mu*s_mod * kcat_km_fenton*h2o2_sym},
                    combine=False)

        # Update compiled expressions (even if basis not updated)
        update_stoich_expr(solver, cons_fenton)

    def add_dna_protection(self, complex_keff_dict={
        'CPLX0-1521':65.}):
        """
        Dps protects DNA from damage
        by hydroxyl radicals?
        """
        me = self.me
        solver = self.me_nlp

        dna_dmg_id = 'DNA_damaged'
        if me.metabolites.has_id(dna_dmg_id):
            dna_dmg = me.metabolites.get_by_id(dna_dmg_id)
        else:
            dna_dmg = Metabolite(dna_dmg_id)
            me.add_metabolites([dna_dmg])

        for cplx_id,keff in iteritems(complex_keff_dict):
            cplx = me.metabolites.get_by_id(cplx_id)
            rxn_id = 'DNA_protect_FWD_'+cplx_id
            if me.reactions.has_id(rxn_id):
                rxn_protect = me.reactions.get_by_id(rxn_id)
            else:
                rxn_protect = Reaction(rxn_id)
                rxn_protect.lower_bound = 0
                rxn_protect.upper_bound = 1000.
                me.add_reaction(rxn_protect)
            # DNA_damaged --(Dps)-> 
            rxn_protect.add_metabolites({
                dna_dmg:-1,
                cplx:-mu/keff/3600.},
                combine=False)

            update_stoich_expr(solver, cplx, rxn_protect)


    def add_free_conc(self, metal):
        """
        Make metal's free concentration an explicit variable so we can constrain it
        """
        solver = self.me_nlp
        me = self.me

        #----------------------------------------------------
        # xfree + \sum_j [Ej:mi]*sij = xtotal
        # or,
        # xfree + \sum_j vdilj/mu*sij = xtotal
        #----------------------------------------------------
        xfree_id = 'Free_concentration_' + metal
        if me.reactions.has_id(xfree_id):
            xfree = me.reactions.get_by_id(xfree_id)
        else:
            xfree = Reaction(xfree_id)
            xfree.lower_bound = 0.
            xfree.upper_bound = 1000.
            me.add_reaction(xfree)

        cons_id = 'constraint_conc_free_' + metal
        if me.metabolites.has_id(cons_id):
            cons = me.metabolites.get_by_id(cons_id)
        else:
            cons = Constraint(cons_id)
            me.add_metabolites([cons])

        xtotal = self.symbol_dict[metal]
        cons._bound = xtotal

        xfree.add_metabolites({cons:1.}, combine=False)

        mod_id = 'mod_'+metal
        for data in me.process_data.get_by_id(mod_id).get_complex_data():
            cplx = data.complex
            dil_dict = self.get_dilution_dict(cplx)
            s_mod = data.subreactions[mod_id]
            for vdil,s_enz in iteritems(dil_dict):
                vdil.add_metabolites({cons: s_enz/mu*s_mod}, combine=False)

        update_stoich_expr(solver, cons)

        return xfree


    def get_metal_conc(self, metal, mufix=None):
        """
        Return metal concentration: total, incorporated, free
        """
        solver = self.me_nlp
        me = self.me
        if mufix is None:
            mufix = me.reactions.biomass_dilution.x
        sym_total = self.symbol_dict[metal]
        total_conc = solver.substitution_dict[sym_total]
        # fe2_free = fe2_total - sum_i Ei si
        incorp_conc = 0.
        mod_id = 'mod_'+metal
        subs_dict = dict(solver.substitution_dict)
        subs_dict['mu'] = mufix
        cplx_metal_dict = {}
        for data in me.process_data.get_by_id(mod_id).get_complex_data():
            cplx = data.complex
            dil_dict = self.get_dilution_dict(cplx)
            s_mod = data.subreactions[mod_id]
            conc_cplxi = 0.
            for vdil,s_enz in iteritems(dil_dict):
                conci = s_enz*vdil.x/mufix*s_mod
                if hasattr(conci,'subs'):
                    conci = conci.subs(subs_dict)
                incorp_conc += conci
                conc_cplxi += conci

            # Also track distribution of incorp metal across proteins
            cplx_metal_dict[cplx.id] = conc_cplxi

        free_conc = total_conc - incorp_conc

        return {'total':total_conc, 'free':free_conc, 'incorporated':incorp_conc,
                'cplx_metal_dict':cplx_metal_dict}


    def bisect(self, symbol, mu_fix, precision=1e-3, pmin=0., pmax=1.0, basis=None,
            check_feas0=False, pzero = 1e-3, solver_precision='quad', maxIter=100, 
            direction='max', verbosity=1):
        """
        popt, hs, xopt, cache = bisect(self, symbol, mu_fix, precision=1e-3, pmin=0., pmax=1.0, basis=None,
            check_feas0=False, pzero = 1e-3, solver_precision='quad', maxIter=100,
            direction='max', verbosity=0)

        Bisect to max any oxidizeme symbolic parameter
        """
        import copy as cp
        solver = self.me_nlp
        me = self.me
        hs = basis
        cache = {}
        solution = None
        # Check feasibility at p=zero?
        if check_feas0:
            sym = self.symbol_dict[symbol]
            solver.substitution_dict[sym] = pzero
            x0, stat0, hs0 = solver.solvelp(
                    mu_fix, verbosity=verbosity-1, precision=solver_precision, basis=hs)
            if me.solution.status is not 'optimal':
                warnings.warn('Infeasible at mu=%g and %s=%g. Returning.'% (
                    zero_mu, symbol, pzero))
                return xzero, hs0, x0, cache
            else:
                hs = hs0

        def checkp(pfix, hs):
            if pfix not in cache:
                sym = self.symbol_dict[symbol]
                solver.substitution_dict[sym] = pfix
                v_new, stat_new, hs_new = solver.solvelp(
                    mu_fix, basis=hs, verbosity=verbosity-1, precision=solver_precision)
                if me.solution.status is 'optimal':
                    hs = hs_new
                stat = me.solution.status
                sol = cp.deepcopy(me.solution)
                cache[pfix] = stat

            return cache[pfix], hs, sol, v_new

        warm = False
        converged = False
        a = pmin
        b = pmax
        muopt = a
        popt = None
        xopt = None
        iter = 0

        tic = time.time()

        if verbosity > 0:
            print(('%-10.5s%-20.17s%-20.17s%-20.17s%-20.17s%-20.17s' % (
                    'iter','popt','a','b','p1','stat1')))
        while iter < maxIter and not converged:
            # Just a sequence of feasibility checks
            p1 = (a+b)/2.
            # Retrieve evaluation from cache if it exists: golden section advantage
            stat1, hs, sol1, x1 = checkp(p1, hs)
            if direction.lower() == 'max':
                if stat1 is 'optimal':
                    a = p1
                    popt = p1
                    solution = sol1
                    xopt = x1
                else:
                    b = p1
            else:
                # Minimize
                if stat1 is 'optimal':
                    b = p1
                else:
                    a = p1
                    popt = p1
                    solution = sol1
                    xopt = x1

            converged = abs(b-a) <= precision
            warm = True
            iter = iter+1

            if verbosity > 0:
                print(('%-10.5s%-20.17s%-20.17s%-20.17s%-20.17s%-20.17s' % (
                        iter, popt, a, b, p1, stat1)))

        toc = time.time()-tic
        if verbosity > 0:
            print(('Bisection done in %g seconds'%toc))

        # Save final solution
        me.solution = solution
        # Save feasible basis
        self.feas_basis = hs

        return popt, hs, xopt, cache


    def fit_omics(self):
        """
        Fit proteome or transcriptome to measured at varying (optimal) resolution
        """


def update_stoich_expr(solver, met=None, rxn=None):
    me = solver.me
    if met is not None:
        mind = me.metabolites.index(met)
        if rxn is not None:
            s = rxn.metabolites[met]
            if isinstance(s,Basic):
                rind = me.reactions.index(rxn)
                expr = solver.compile_expr(s)
                solver.compiled_expressions[(mind,rind)] = expr
        else:
            # Update ._bound?
            if isinstance(met._bound, Basic):
                expr = solver.compile_expr(met._bound)
                solver.compiled_expressions[(mind,None)] = (expr, met._constraint_sense)
            for rxn in met.reactions:
                s = rxn.metabolites[met]
                if isinstance(s,Basic):
                    rind = me.reactions.index(rxn)
                    expr = solver.compile_expr(s)
                    solver.compiled_expressions[(mind,rind)] = expr
    elif rxn is not None:
        rind = me.reactions.index(rxn)
        for met,s in iteritems(rxn.metabolites):
            if isinstance(s,Basic):
                mind = me.metabolites.index(met)
                expr = solver.compile_expr(s)
                solver.compiled_expressions[(mind,rind)] = expr
