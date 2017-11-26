# BayesVarSel 1.7.1.9000

* Merged Bvs and PBvs in just one function (called Bvs and old PBvs disappears). Bvs now has two extra parameters to control parallelization: parallel and n.nodes.

* Several changes in Btest: 1) it now allows for unnamed lists of models and if unnamed lists are provided, default names are given by the function. 2) The prior probabilities argument priorprobs does not have to be named (by default the order of the models in argument models is used). 3) Deprecated argument relax.nest, replaced by the explicit definition of the null model by the user via the argument null.model
