py_help.py
	Defines useful functions and constants (priors, parameter types, data names, etc.)



make_posterior.py
	Supporting function for make_likelihood.py, defines the makePost function.
	This takes a folder and iterates through all Trace_XX_YY.npz files (LOO iterates)
	compiling the ensemble of predictions given in the trace object into a posterior function.
	This is then compared to the observed predictand value to get a likelihood for each LOO run.

make_likelihood.py
	Iterates through an experiment folder, running makePost on each subfolder (choice of vars).
	The resulting likelihoods for each model are saved in the allLikelihoods.npz file
	in the experiment folder.

test_pred2018.py
	For each model, makes a trace for each time step after 1/1/2018, using all prior data.
	This creates a projection for the predictand (a posterior). 

test_bma2018.py
	These are averaged according to the weights given by allLikelihoods.npz.

