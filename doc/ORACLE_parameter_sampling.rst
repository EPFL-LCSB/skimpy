ORACLE parameter sampling 
===========================
We here present the first open-source  implementation of the ORACLE parametrization method [3-6]. ORACLE makes use of sampling alternative steady states and parameter backcalculation to efficiently sample large sets of locally stable parameters.


The first step in the ORACLE workflow is to integrate physiological data. Here we integrare observed concentrations and we assume that the glucose transporter operates far from equilibrium.

.. code-block:: python

  from pytfa.io.json import load_json_model, save_json_model

  path_to_tmodel = './../../models/tfa_varma.json'
  tmodel = load_json_model(path_to_tmodel)

  GLPK= 'optlang-glpk'
  CPLEX = 'optlang-cplex'
  tmodel.solver = GLPK

  tmodel.solver.configuration.tolerances.feasibility = 1e-9

  # Define the reactor medium
  #Glucose  to 10 g/L = 0.056 mol/l 1/2 Reactor concentration
  GLUCOSE_CONCENTRATION = 0.056
  PHOSPHATE = 10e-3
  CARBON_DIOXIDE = 1e-7
  OXYGEN = 8e-3*0.062 # 8 mg/L 1g = 0.062 mol

  tmodel.log_concentration.get_by_id('glc-D_e').variable.ub = np.log(GLUCOSE_CONCENTRATION*1.2)
  tmodel.log_concentration.get_by_id('glc-D_e').variable.lb = np.log(GLUCOSE_CONCENTRATION*0.8)

  tmodel.log_concentration.get_by_id('pi_e').variable.lb = np.log(PHOSPHATE*0.8)
  tmodel.log_concentration.get_by_id('pi_e').variable.ub = np.log(PHOSPHATE*1.2)

  tmodel.log_concentration.get_by_id('co2_e').variable.ub = np.log(CARBON_DIOXIDE*1.2)
  tmodel.log_concentration.get_by_id('co2_e').variable.lb = np.log(CARBON_DIOXIDE*0.8)

  tmodel.log_concentration.get_by_id('o2_e').variable.ub = np.log(OXYGEN*1.2)
  tmodel.log_concentration.get_by_id('o2_e').variable.lb = np.log(OXYGEN*0.8)

  # Enforce glucose transporter displacements
  tmodel.thermo_displacement.GLCptspp.variable.ub = -2.0

  # Test feasiblity
  print(tmodel.optimize())


In the next step we sample the flux and concentration space:

.. code-block:: python

  from pytfa.analysis.sampling import sample

  NUM_TFA_SAMPLES = 10

  samples = sample(tmodel, NUM_TFA_SAMPLES, method='achr')
  samples.to_csv('./samples_fdp1_1000.csv'.format())
  

Fluxes sampled in the TFA problem are mass fluxes in units of mmol/gDW/hr we convert these
to reaction rates in units of [concentration]/[time] therefore we require the cell density and
the wet weight to dry weight ratio. It is also useful to choose a different concentration unit than [M] since
most intracellular metabolite concentration are between 1 mM and 1 uM, which can be achived by choosing a different
concentration scaling factor: ``[C] = [M] x CONCENTRATION_SCALING``.

.. code-block:: python

  N_PARAMETER_SAMPLES = 10
  CONCENTRATION_SCALING = 1e6 # 1 umol
  TIME_SCALING = 1 # 1hr
  DENSITY = 1200 # g/L
  GDW_GWW_RATIO = 0.3 # Assumes 70% Water


Finally we call the parameter sampler for each sampled reference steady-state. The parameter sampler provides then information on the charateristic timescales this information can then later be used to discard parameter samples exhibiting unphysiological characteristics.
The parameter sampler will requires that the Jacobian is compiled ``kmodel.compile_jacobian()`` to evaluate the stablity of the samples.

.. code-block:: python

  from skimpy.io.yaml import load_yaml_model
  from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
      load_concentrations, load_equilibrium_constants
  from skimpy.core.parameters import ParameterValuePopulation, load_parameter_population
  from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
  from skimpy.utils.general import get_stoichiometry


  NCPU = 8

  path_to_kmodel = './../../models/kin_varma.yml'
  path_for_output = './paramter_pop_{}.h5'
  kmodel = load_yaml_model(path_to_kmodel)

  # Perp and compile to sample parameters
  kmodel.prepare()
  kmodel.compile_jacobian(ncpu=NCPU)

  sampler_params = SimpleParameterSampler.Parameters(n_samples=N_PARAMETER_SAMPLES)
  sampler = SimpleParameterSampler(sampler_params)

  lambda_max_all = []
  lambda_min_all = []

  S = get_stoichiometry(kmodel, kmodel.reactants).todense()

  from skimpy.utils.namespace import *
  from skimpy.analysis.ode.utils import make_flux_fun
  fluxfun= make_flux_fun(kmodel, QSSA)
  fluxes = []

  for i, sample in samples.iterrows():
      # Load fluxes and concentrations and scale them to the desired units
      fluxes = load_fluxes(sample, tmodel, kmodel,
                           density=DENSITY,
                           ratio_gdw_gww=GDW_GWW_RATIO,
                           concentration_scaling=CONCENTRATION_SCALING,
                           time_scaling=TIME_SCALING)

      concentrations = load_concentrations(sample, tmodel, kmodel,
                                           concentration_scaling=CONCENTRATION_SCALING)


      # Fetch equilibrium constants and scale them to the desired units
      load_equilibrium_constants(sample, tmodel, kmodel,
                                 concentration_scaling=CONCENTRATION_SCALING,
                                 in_place=True)


      # Generate samples and fetch slowest and fastest eigenvalues
      params, lamda_max, lamda_min = sampler.sample(kmodel, fluxes, concentrations,
                                                    only_stable=True,
                                                    min_max_eigenvalues=True)
      lambda_max_all.append(pd.DataFrame(lamda_max))
      lambda_min_all.append(pd.DataFrame(lamda_min))

      params_population = ParameterValuePopulation(params, kmodel)
      params_population.save(path_for_output.format(i))


  # Process df and save dataframe
  lambda_max_all = pd.concat(lambda_max_all, axis=1)
  lambda_min_all = pd.concat(lambda_min_all, axis=1)

