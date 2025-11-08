**Summary**:

- The folders sim and real contain the full setups for the simulation studies and real data analysis presented in the manuscript.
* sim: Runs 100 repetitions comparing different methods.
* real: Contains the real data analysis.

**Quick Start**:

- To reproduce results for each method in the simulation or real data analysis, run 'main.R' in the corresponding folder.
- To change the number of observations for the simulated system, modify the variable 'obs_choice' or 'N' (default: 41 observations).
- To change the noise level, adjust the variable 'snr' or 'SNR' (default: SNR = 25).



**Simulation**: 

- The 'simulation' folder includes four subfolders, each corresponding to a specific simulation study:
* add_ode: Benchmark system (17) from Section 3 of the main text.
* large_ode: Larger version of system (17), described in Web Appendix B.1.
* abs_ode: ODE system including an absolute value function, described in Web Appendix B.2.
* interaction_ode: Lotka–Volterra equations, described in Web Appendix B.3.

- Each '*_ode' folder contains:
* Bayes_int: Our proposed Bayesian integral method (as described in the manuscript).
* other_methods: Contains implementations of:
 - JADE: Joint estimation approach for generalized sparse additive ODEs [1].
 - GRADE: Graph reconstruction via additive differential equations [2].
 The code is adapted from the original authors.

- Contents of 'Bayes_int' folder:
* main.R: Main script to run simulations. Modify 'obs_choice' and 'snr' as needed.
* prep.R: Data preparation and setup script.
* help_func.R: Helper functions for simulations and result processing.
* BayesODE.R: Core function for fitting the simulation data.
* compare.R: Script to summarize and compare results.
* seeds.txt: List of random seeds used for data generation.

- Residual boostrap for JADE (add_ode/other_methods/bootstrap.R):
* bootstrap.R: run the residual bootstrap for JADE and calculate the coverage probabilities of trajectories. 


**Real Data**: 

- This folder contains our proposed integral method for real data in the manuscript, which has the following files:
* main.R: Main script to run the analysis.
* main_prior.R: Main script to assess the sensitivity to priors (Web Appendix C).
* prep.R: Data preparation script.
* help_func.R: Helper functions for processing and analysis.
* BayesODE.R: Core function for model fitting.
* compare.R: Script to process and summarize the results.
* BayesODE_prior.R: Core function for model fitting with different prior choices (Web Appendix C).
* compare_prior.R: Script to process and summarize the results for different prior choices (Web Appendix C).

- The real dataset is sourced from the R package 'longitudinal' [3].


**Reference**:
[1] Zhang, N., Nanshan, M., and Cao, J. (2022). A joint estimation approach to sparse additive ordinary differential equations. Statistics and Computing 32, 69.
[2] Chen, S., Shojaie, A., and Witten, D. M. (2017). Network reconstruction from high-dimensional ordinary differential equations. Journal of the American Statistical Association 112, 1697–1707.
[3] Strimmer, K. (2022). longitudinal: Analysis of Multiple Time Course Data. R package version 1.1.13, <https://CRAN.R-project.org/package=longitudinal>.
