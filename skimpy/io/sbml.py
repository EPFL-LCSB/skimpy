# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2022 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""
from skimpy.core.kinmodel import KineticModel
from skimpy.utils.namespace import PARAMETER, VARIABLE

from skimpy.mechanisms.direct_expression import make_mechanism_from_expression
from skimpy.utils.general import make_subclasses_dict, get_stoichiometry, get_all_reactants


from libsbml import *
import sys

from skimpy.core import Item, Reactant, Parameter, Reaction, BoundaryCondition, \
    ConstantConcentration, KineticModel, ExpressionModifier
from skimpy.core.compartments import Compartment

# DEFAULT UNITS:
TIME = 'hour'
SUBSTANCE = 'micromole'
VOLUME = 'litre'


"""
Check function adapted as recommended from SBML-examples:
See: https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/create_simple_model_8py-example.html
"""
def check(value, message):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    if value is None:
        raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
    elif type(value) is int:
        if value == LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error encountered trying to ' + message + '.' \
                      + 'LibSBML returned error code ' + str(value) + ': "' \
                      + OperationReturnValue_toString(value).strip() + '"'
            raise SystemExit(err_msg)
    else:
        return



def export_sbml(kmodel, filename, time_unit=TIME, substance_unit=SUBSTANCE, volume_unit=VOLUME):

    # compile odes to make sure all rate expression are computed correctly
    kmodel.compile_ode()

    try:
        document = SBMLDocument(3, 1)
    except ValueError:
        raise SystemExit('Could not create SBMLDocument object')

    model = document.createModel()
    check(model, 'create SMBL model')
    check(model.setTimeUnits(time_unit), 'set model-wide time units')
    check(model.setExtentUnits(substance_unit), 'set model units of extent')
    check(model.setSubstanceUnits(substance_unit), 'set model substance units')

    add_units(model, time_unit=time_unit, substance_unit=substance_unit, volume_unit=volume_unit)

    has_compartments = True if kmodel.compartments else False
    if not has_compartments:
        # If model does not contain compartment definition we make default compartment
        # of size 1
        c = model.createCompartment()
        check(c, 'create compartment')
        check(c.setId('c'), 'set compartment id')
        check(c.setConstant(True), 'set compartment "constant"')
        check(c.setSize(1), 'set compartment "size"')
        check(c.setSpatialDimensions(3), 'set compartment dimensions')
        check(c.setUnits(volume_unit), 'set compartment size units')

    else:
        for compartment in kmodel.compartments:
            c = model.createCompartment()
            check(c, 'create compartment')
            check(c.setId(compartment.name), 'set compartment id')
            check(c.setConstant(True), 'set compartment "constant"')
            check(c.setSize(compartment.volume), 'set compartment "size"')
            check(c.setSpatialDimensions(3), 'set compartment dimensions')
            check(c.setUnits(volume_unit), 'set compartment size units')

    for reactant in get_all_reactants(kmodel).values():

        s = model.createSpecies()
        check(s, 'create species s1')
        check(s.setId(reactant.name), 'set species s1 id')
        if has_compartments:
            check(s.setCompartment(reactant.compartment.name), 'set species s1 compartment')
        else:
            check(s.setCompartment('c'), 'set species {} compartment'.format(reactant.name))

        # Set parameter value
        if reactant.type == VARIABLE:
            check(s.setConstant(False), 'set "constant" attribute on {}'.format(reactant.name))
            check(s.setBoundaryCondition(False), 'set "boundaryCondition" on {}'.format(reactant.name))
            try:
                check(s.setInitialAmount(kmodel.initial_conditions[reactant.name]),
                        'set initial amount for {}'.format(reactant.name))
            except KeyError:
                check(s.setInitialAmount(0),'set initial amount for {}'.format(reactant.name))
        else:
            check(s.setConstant(True), 'set "constant" attribute on {}'.format(reactant.name))
            check(s.setBoundaryCondition(True), 'set "boundaryCondition" on {}'.format(reactant.name))
            check(s.setInitialAmount(kmodel.parameters[reactant.name].value),
                  'set initial amount for {}'.format(reactant.name))

        check(s.setSubstanceUnits(substance_unit), 'set substance units for {}'.format(reactant.name))
        check(s.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on {}'.format(reactant.name))

    # Create a parameter object inside this model, set the required
    # attributes 'id' and 'constant' for a parameter in SBML Level 3, and
    # initialize the parameter with a value along with its units.
    all_reactants = get_all_reactants(kmodel)
    for parameter in kmodel.parameters.values():
        unique_parameter_name = str(parameter.symbol)
        if unique_parameter_name in all_reactants:
            # Skip boundary condition parameters
            continue

        k = model.createParameter()
        check(k, 'create parameter {}'.format(unique_parameter_name))
        check(k.setId(unique_parameter_name), 'set parameter {} id'.format(unique_parameter_name))
        check(k.setConstant(True), 'set parameter {} "constant"'.format(unique_parameter_name))

        value = parameter.value if not parameter.value is None else 0.0
        check(k.setValue(value), 'set parameter {} value'.format(unique_parameter_name))

        unit = get_unit(unique_parameter_name)
        check(k.setUnits(unit), 'set parameter {} units'.format(unique_parameter_name))

    for reaction in kmodel.reactions.values():

        r1 = model.createReaction()
        check(r1, 'create reaction')
        check(r1.setId(reaction.name), 'set reaction id')
        reversible = get_reversiblity(reaction)
        check(r1.setReversible(reversible), 'set reaction reversibility flag')
        check(r1.setFast(False), 'set reaction "fast" attribute')

        for reactant, stoich in reaction.reactant_stoichiometry.items():
            if stoich > 0:
                species_ref1 = r1.createReactant()
                check(species_ref1, 'create reactant')
                check(species_ref1.setSpecies(reactant.name), 'assign reactant species')
                check(species_ref1.setStoichiometry(abs(stoich)), 'assing stoichiometry')
                check(species_ref1.setConstant(True), 'set "constant" on species ref 1')

            elif stoich < 0:
                species_ref2 = r1.createProduct()
                check(species_ref2, 'create product')
                check(species_ref2.setStoichiometry(abs(stoich)), 'assing stoichiometry')
                check(species_ref2.setSpecies(reactant.name), 'assign product species')
                check(species_ref2.setConstant(True), 'set "constant" on species ref 2')

        # Add a kinetic rate law
        math_ast = parseL3Formula( str(reaction.mechanism.reaction_rates['v_net']) )
        check(math_ast, 'create AST for rate expression')

        kinetic_law = r1.createKineticLaw()
        check(kinetic_law, 'create kinetic law')
        check(kinetic_law.setMath(math_ast), 'set math on kinetic law')

    # Write the SBML-file
    writeSBML(document, filename)


def import_sbml(filename):
    reader = SBMLReader()
    document = reader.readSBML(filename)

    if document.getNumErrors() > 0:
        raise RuntimeError("Reading the SBML file raised several errors")

    sbml_model = document.getModel()

    new = KineticModel()

    # Rebuild the reactions
    for the_reaction in sbml_model.reactions:
        expression = the_reaction.getKineticLaw().getFormula()
        products = [(p.getSpecies(), p.getStoichiometry()) for p in the_reaction.products]
        substrates = [(p.getSpecies(), -p.getStoichiometry()) for p in the_reaction.reactants]

        reactants, stoichiometry = list(zip(*(substrates+products)))

        TheMechanism = make_mechanism_from_expression(stoichiometry, reactants, expression )
        the_reactants = TheMechanism.Reactants(**{k:v for k,v in zip(TheMechanism.reactant_list, reactants)})

        # FIXME: Needs a corresponding input in the SBML?

        the_enzyme = None

        new_reaction = Reaction(name=the_reaction.id,
                                mechanism=TheMechanism,
                                reactants=the_reactants,
                                enzyme = the_enzyme,
                                )

        new.add_reaction(new_reaction)
        
    #Resbuild compartments
    for the_comp in sbml_model.compartments:
        new_comp = Compartment(the_comp.id)
        # Set comp size parameter
        new_comp.parameters.volume.value = the_comp.size
        new_comp.parameters.cell_volume.value = 1.0
        new.add_compartment(new_comp)


    # Assing compartments
    for the_met in sbml_model.species:
        met = new.reactants[the_met.id]
        comp = new.compartments[the_met.compartment]
        met.compartment = comp


    # Populate the kinmodel.parameters TabDict
    parameter_init_dict = dict()
    for rxn_obj in new.reactions.values():
        # initalize empty param list
        parameter_init_dict[rxn_obj.name] = rxn_obj.mechanism.__class__.Parameters()
    new.parametrize_by_reaction(parameter_init_dict)

    # Boundary Conditions
    for the_met in sbml_model.species:
        if the_met.boundary_condition:
            if the_met.constant:
                TheBoundaryCondition = ConstantConcentration
            else:
                raise NotImplementedError('TODO Implement flux boundary condition')

            the_bc = TheBoundaryCondition(reactant, **the_bc_dict)
            new.add_boundary_condition(the_bc)

            # Do not forget to add the value of the BC!
            reactant.value = the_met.id

    # Fetch parameter values
    for the_param in sbml_model.parameters:
        try:
            new.parameters[the_param.id].value = the_param.value
        except KeyError:
            pass

    return new


# Utility Functions
def get_reversiblity(reaction):
    if reaction.mechanism.reaction_rates['v_bwd'] ==0:
        return False
    else:
        return True


# TODO Implement this object oriented ...
def get_unit(parameter_name, substance_unit=SUBSTANCE, ):

    if 'KM' in parameter_name\
            or 'KI' in parameter_name\
            or 'KA' in parameter_name:
        return substance_unit

    elif 'kcat' in parameter_name:
        return 'per_time'

    # THIS NEEDS TO CHANGE!!!! FIGURE OUT HOW
    elif 'vmax' in parameter_name:
        return 'per_time'

    # WE ALSO NEED TO DET THE UNITS FOR EQUILIBRIUM CONSTANTS
    # THIS SHOULD PROBABLY HANDLED BY SKIMMPY ...
    else:
        return ''



def add_units(model,time_unit=TIME, substance_unit=SUBSTANCE, volume_unit=VOLUME):

    per_time = model.createUnitDefinition()
    check(per_time, 'create unit definition')
    check(per_time.setId('per_time'), 'set unit definition id')
    unit = per_time.createUnit()
    check(unit, 'create unit on per_second')
    check(unit.setKind(UNIT_KIND_SECOND), 'set unit kind')
    check(unit.setMultiplier(1), 'set unit multiplier')
    check(unit.setExponent(-1), 'set unit exponent')
    check(unit.setScale(0), 'set unit scale')

