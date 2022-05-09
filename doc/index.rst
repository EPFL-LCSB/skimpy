.. skimpy documentation master file, created by
   sphinx-quickstart on Sat Jul  1 12:44:20 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the (not so) SKiMPy documentation!
=================================

SKiMPy (Symbolic Kinetic Models in Python) implements various methods and resources that allow the user to build and analyze large scale kinetic models efficiently, including i) automated draft model construction from FBA and TFA using Cobrapy [1] and pyTFA [2], respectively, ii) an extensive library of kinetic mechanisms, iii) efficient model parametrization using ORACLE [3–7], iv) local stability analysis [8], v) global stability analysis, vi) local sensitivity analysis, e.g., metabolic control analysis [9], vii) global sensitivity analysis and uncertainty propagation [10,11] viii) modal analysis [8], iix) non-linear ODE integration [12,13], ix) identification of conserved pools [14]. Therefore, the presented Python package implements an object-oriented interface to construct the symbolic expressions from a library of kinetic mechanisms using sympy [15] and then precompiling these expressions into machine code using Cython [16]. SKiMpy also integrates the SUNDIALS ODE-Solver package [13] using the interface provided by the ODES package [17].

For installtion instructions please refer to `README.rst <https://github.com/EPFL-LCSB/skimpy/blob/master/README.rst>`_

.. toctree::
   :numbered:
   :maxdepth: 3
   :caption: Contents:

   quickstart
   reading_writing_files
   build_draft_models
   ORACLE_parameter_sampling
   metabolic_control_analysis
   modal_analysis
   dynamic_simulations 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


References to methods
======================

.. [1] Ebrahim, A.; Lerman, J. A.; Palsson, B. O.; Hyduke, D. R. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology 2013, 7 (1), 74. https://doi.org/10.1186/1752-0509-7-74.
.. [2] 	Salvy, P.; Fengos, G.; Ataman, M.; Pathier, T.; Soh, K. C.; Hatzimanikatis, V. PyTFA and MatTFA: A Python Package and a Matlab Toolbox for Thermodynamics-Based Flux Analysis. Bioinformatics 2019, 35 (1), 167–169. https://doi.org/10.1093/bioinformatics/bty499.
.. [3] 	Wang, L.; Birol, I.; Hatzimanikatis, V. Metabolic Control Analysis under Uncertainty: Framework Development and Case Studies. Biophys J 2004, 87 (6), 3750–3763. https://doi.org/10.1529/biophysj.104.048090.
.. [4] 	Miskovic, L.; Hatzimanikatis, V. Production of Biofuels and Biochemicals: In Need of an ORACLE. Trends in Biotechnology 2010, 28 (8), 391–397. https://doi.org/10.1016/j.tibtech.2010.05.003.
.. [5] 	Chakrabarti, A.; Miskovic, L.; Soh, K. C.; Hatzimanikatis, V. Towards Kinetic Modeling of Genome-Scale Metabolic Networks without Sacrificing Stoichiometric, Thermodynamic and Physiological Constraints. Biotechnol J 2013, 8 (9), 1043–1057. https://doi.org/10.1002/biot.201300091.
.. [6] 	Savoglidis, G.; da Silveira Dos Santos, A. X.; Riezman, I.; Angelino, P.; Riezman, H.; Hatzimanikatis, V. A Method for Analysis and Design of Metabolism Using Metabolomics Data and Kinetic Models: Application on Lipidomics Using a Novel Kinetic Model of Sphingolipid Metabolism. Metab Eng 2016, 37, 46–62. https://doi.org/10.1016/j.ymben.2016.04.002.
.. [7] 	Tokic, M.; Hatzimanikatis, V.; Miskovic, L. Large-Scale Kinetic Metabolic Models of Pseudomonas Putida KT2440 for Consistent Design of Metabolic Engineering Strategies. Biotechnology for Biofuels 2020, 13 (1), 33. https://doi.org/10.1186/s13068-020-1665-7.
.. [8] 	Heinrich, R.; Schuster, S. The Regulation of Cellular Systems; Springer Science & Business Media, 2012.
.. [9] 	Kacser, H.; Burns, J. A. The Control of Flux. Symp. Soc. Exp. Biol. 1973, 27, 65–104.
.. [10] 	Hameri, T. E. Towards Comprehensive and Consistent Kinetic Models of Metabolism under Uncertainty. Thèse de doctorat, EPFL, Lausanne, 2019. https://doi.org/10.5075/epfl-thesis-9194.
.. [11] 	Denhardt-Eriksson, R. A. Modelling the Intracellular Environment under Uncertainty, from Post-Translational Modification to Metabolic Networks. Thèse de doctorat, EPFL, Lausanne, 2021. https://doi.org/10.5075/epfl-thesis-8764.
.. [12] 	Serban, R.; Hindmash, A. C. CVODES, the Sensitivity-Enabled ODE Solver in SUNDIALS. Proceedings of the ASME International Design Engineering Technical Conferences and Computers and Information in Engineering Conference, Vol 6, Pts A-C 2005, 257–269.
.. [13] 	Hindmarsh, A. C.; Brown, P. N.; Grant, K. E.; Lee, S. L.; Serban, R.; Shumaker, D. E.; Woodward, C. S. SUNDIALS: Suite of Nonlinear and Differential/Algebraic Equation Solvers. Acm Transactions on Mathematical Software 2005, 31 (3), 363–396. https://doi.org/Doi 10.1145/1089014.1089020.
.. [14] 	Nikolaev, E. V.; Burgard, A. P.; Maranas, C. D. Elucidation and Structural Analysis of Conserved Pools for Genome-Scale Metabolic Reconstructions. Biophysical Journal 2005, 88 (1), 37–49. https://doi.org/10.1529/biophysj.104.043489.
.. [15] 	Meurer, A.; Smith, C. P.; Paprocki, M.; Čertík, O.; Kirpichev, S. B.; Rocklin, M.; Kumar, Am.; Ivanov, S.; Moore, J. K.; Singh, S.; Rathnayake, T.; Vig, S.; Granger, B. E.; Muller, R. P.; Bonazzi, F.; Gupta, H.; Vats, S.; Johansson, F.; Pedregosa, F.; Curry, M. J.; Terrel, A. R.; Roučka, Š.; Saboo, A.; Fernando, I.; Kulal, S.; Cimrman, R.; Scopatz, A. SymPy: Symbolic Computing in Python. PeerJ Comput. Sci. 2017, 3, e103. https://doi.org/10.7717/peerj-cs.103.
.. [16] 	Behnel, S.; Bradshaw, R.; Citro, C.; Dalcin, L.; Seljebotn, D. S.; Smith, K. Cython: The Best of Both Worlds. Computing in Science Engineering 2011, 13 (2), 31–39. https://doi.org/10.1109/MCSE.2010.118.
.. [17] 	Malengier, B.; Kišon, P.; Tocknell, J.; Abert, C.; Bruckner, F.; Bisotti, M.-A. ODES: A High Level Interface to ODE and DAE Solvers. Journal of Open Source Software 2018, 3 (22), 165. https://doi.org/10.21105/joss.00165.
