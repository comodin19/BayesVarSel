# BayesVarSel 2.2.1
* Two versions of the intrinsic hyper-g prior are implemented.

# BayesVarSel 2.2.0
* the OBICE study is included.

# BayesVarSel 2.1.0
* intrinsic prior added to Bvs and GibbsBvs.

* Now calculations of Bayes factors for liang and ZS prior use a definition based on the exp(log()) of the expression: quite more stable

# BayesVarSel 2.0.1
* Minor updates in order to pass CRAN checks.

# BayesVarSel 2.0.0

* Added new functionalities to handle the situation with p>n and p>>n

* Added a new function to handle variable selection with factors that is independent on the parameterization (contrast) used to parameterize factors and that is able to control for multiplicity.

* Corrected a numerical error in the computation of Robust Bayes factors.

# BayesVarSel 1.9.0

* Added en estimation of posterior probabilities of models based on the normalizing constant. 
(it goes with the print function)

# BayesVarSel 1.8.2

* Added a vignette based on the R Journal article https://journal.r-project.org/archive/2018/RJ-2018-021/index.html.

# BayesVarSel 1.8.1

* Added the calculation of the normalizing constant.

# BayesVarSel 1.8.0

* Merged Bvs and PBvs in just one function (called Bvs and old PBvs disappears). Bvs now has two extra parameters to control parallelization: parallel and n.nodes.

* Several changes in Btest: 1) it now allows for unnamed lists of models and if unnamed lists are provided, default names are given by the function. 2) The prior probabilities argument priorprobs does not have to be named (by default the order of the models in argument models is used). 3) Deprecated argument relax.nest, replaced by the explicit definition of the null model by the user via the argument null.model

* In Bvs and GibbsBvs the order of arguments has slightly changed and now data is the second argument followed by formula.

* Now plotBvs is a S3 function defined as plot.Bvs

* Now predictBvs is a S3 function defined as predict.Bvs

* In Bvs the argument fixed.cov is deprecated, now replaced by the formula for the null model in the argument null.model

* Option "trace" added to plot

* print updated to show the 10 most probable models (among the visited) for method Gibbs

* removed the comment at initialization

* BMAcoeff no longer shows the graphic for all the variables, instead use histBMA to plot the posterior distribution of the coefficients. 
