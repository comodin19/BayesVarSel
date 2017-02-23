`BayesVarSel`
=============

Bayes Factors, Model Choice And Variable Selection In Linear Models
-------------------------------------------------------------------

Hypothesis testing, model selection and model averaging are important
statistical problems that have in common the explicit consideration of
the uncertainty about which is the true model. The formal Bayesian tool
to solve such problems is the Bayes factor (Kass and Raftery, 1995) that
reports the evidence in the data favoring each of the entertained
hypotheses/models and can be easily translated to posterior
probabilities.

This package has been specifically conceived to calculate Bayes factors
in linear models and then to provide a formal Bayesian answer to testing
and variable selection problems. From a theoretical side, the emphasis
in the package is placed on the prior distributions (a very delicate
issue in this context) and BayesVarSel allows using a wide range of
them: Jeffreys-Zellner-Siow (Jeffreys, 1961; Zellner and Siow,
1980,1984) Zellner (1986); Fernandez et al. (2001), Liang et al. (2008)
and Bayarri et al. (2012).

The interaction with the package is through a friendly interface that
syntactically mimics the well-known lm command of R. The resulting
objects can be easily explored providing the user very valuable
information (like marginal, joint and conditional inclusion
probabilities of potential variables; the highest posterior probability
model, HPM; the median probability model, MPM) about the structure of
the true -data generating- model. Additionally, BayesVarSel incorporates
abilities to handle problems with a large number of potential
explanatory variables through parallel and heuristic versions
(Garcia-Donato and Martinez-Beneito 2013) of the main commands.

References
----------

-   Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G. (2012).
    Criteria for Bayesian Model choice with Application to Variable
    Selection. The Annals of Statistics. 40: 1550-1577. DOI:
    [10.1214/12-aos1013](http://www.dx.doi.org/10.1214/12-aos1013)
-   Fernandez, C., Ley, E. and Steel, M.F.J. (2001). Benchmark priors
    for Bayesian model averaging. Journal of Econometrics, 100, 381-427.
    DOI:
    [10.1016/s0304-4076(00)00076-2](http://www.dx.doi.org/10.1016/s0304-4076(00)00076-2)
-   Garcia-Donato, G. and Martinez-Beneito, M.A. (2013). On sampling
    strategies in Bayesian variable selection problems with large model
    spaces. Journal of the American Statistical Association. 108:
    340-352.DOI:
    [10.1080/01621459.2012.742443](http://www.dx.doi.org/10.1080/01621459.2012.742443)
-   Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2008).
    Mixtures of g-priors for Bayesian Variable Selection. Journal of the
    American Statistical Association. 103:410-423.DOI:
    [10.1198/016214507000001337](http://www.dx.doi.org/10.1198/016214507000001337)
-   Zellner, A. and Siow, A. (1980). Posterior Odds Ratio for Selected
    Regression Hypotheses. In Bayesian Statistics 1 (J.M. Bernardo, M.
    H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603.
    Valencia: University Press. DOI:
    [10.1007/bf02888369](http://www.dx.doi.org/10.1007/bf02888369)
-   Zellner, A. and Siow, A. (1984) Basic Issues in Econometrics.
    Chicago: University of Chicago Press.
-   Zellner, A. (1986) On Assessing Prior Distributions and Bayesian
    Regression Analysis with g-prior Distributions. In Bayesian
    Inference and Decision techniques: Essays in Honor of Bruno de
    Finetti (A. Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
    DOI: [10.2307/2233941](http://www.dx.doi.org/10.2307/2233941)
