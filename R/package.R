
#' @import parallel MASS mvtnorm
#' @useDynLib BayesVarSel
#' @importFrom grDevices gray.colors gray
#' @importFrom graphics axis barplot image layout par plot hist text
#' @importFrom stats as.formula lm rbinom runif density quantile
#' @importFrom utils read.table write.table
#' @importFrom stats model.frame model.matrix terms

NULL


#' Bayes Factors, Model Choice And Variable Selection In Linear Models
#'
#' Hypothesis testing, model selection and model averaging are important
#' statistical problems that have in common the explicit consideration of the
#' uncertainty about which is the true model. The formal Bayesian tool to solve
#' such problems is the Bayes factor (Kass and Raftery, 1995) that reports the
#' evidence in the data favoring each of the entertained hypotheses/models and
#' can be easily translated to posterior probabilities.
#'
#' This package has been specifically conceived to calculate Bayes factors in
#' linear models and then to provide a formal Bayesian answer to testing and
#' variable selection problems. From a theoretical side, the emphasis in the
#' package is placed on the prior distributions (a very delicate issue in this
#' context) and BayesVarSel allows using a wide range of them:
#' Jeffreys-Zellner-Siow (Jeffreys, 1961; Zellner and Siow, 1980,1984) Zellner
#' (1986); Fernandez et al. (2001), Liang et al. (2008) and Bayarri et al.
#' (2012).
#'
#' The interaction with the package is through a friendly interface that
#' syntactically mimics the well-known lm command of R. The resulting objects
#' can be easily explored providing the user very valuable information (like
#' marginal, joint and conditional inclusion probabilities of potential
#' variables; the highest posterior probability model, HPM; the median
#' probability model, MPM) about the structure of the true -data generating-
#' model. Additionally, BayesVarSel incorporates abilities to handle problems
#' with a large number of potential explanatory variables through parallel and
#' heuristic versions (Garcia-Donato and Martinez-Beneito 2013) of the main
#' commands.
#'
#' \tabular{ll}{ Package: \tab BayesVarSel\cr Type: \tab Package\cr Version:
#' \tab 2.0.1\cr Date: \tab 2020-02-17\cr License: \tab GPL-2\cr }
#'
#' @name BayesVarSel-package
#' @aliases BayesVarSel-package BayesVarSel
#' @docType package
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: Anabel Forte \email{anabel.forte@@uv.es}
#' @seealso \code{\link[BayesVarSel]{Btest}}, \code{\link[BayesVarSel]{Bvs}},
#' \code{\link[BayesVarSel]{GibbsBvs}},
#' \code{\link[BayesVarSel]{BMAcoeff}}, \code{\link[BayesVarSel]{predict.Bvs}}
#' @references
#'
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1577
#'
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
#' Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association. 108: 340-352.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423.
#'
#' Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369>. Posterior Odds
#' Ratio for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M.
#' Bernardo, M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603.
#' Valencia: University Press.
#'
#' Zellner, A. and Siow, A. (1984) Basic Issues in Econometrics. Chicago:
#' University of Chicago Press.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#' @keywords package
#' @examples
#' demo(BayesVarSel.Hald)
#'
NULL


#' Hald data
#'
#' The following data relates to an engineering application that was interested
#' in the effect of the cement composition on heat evolved during hardening
#' (for more details, see Woods et al., 1932).
#'
#'
#' @name Hald
#' @docType data
#' @format A data frame with 13 observations on the following 5 variables.
#' \describe{ \item{y}{Heat evolved per gram of cement (in calories)}
#' \item{x1}{Amount of tricalcium aluminate} \item{x2}{Amount
#' of tricalcium silicate} \item{x3}{Amount of tetracalcium alumino
#' ferrite} \item{x4}{Amount of dicalcium silicate} }
#' @references Woods, H., Steinour, H. and Starke, H.
#' (1932)<DOI:10.1021/ie50275a002> Effect of Composition of Porland Cement on
#' Heat Evolved During Hardening. Industrial and Engineering Chemistry
#' Research, 24, 1207-1214.
#' @keywords datasets
#' @examples
#' data(Hald)
#'
"Hald"


#' Ozone35 dataset
#'
#' Polution data
#'
#' This dataset has been used by Garcia-Donato and Martinez-Beneito (2013) to
#' illustrate the potential of the Gibbs sampling method (in \code{BayesVarSel}
#' implemented in \code{\link[BayesVarSel]{GibbsBvs}}).
#'
#' This data were previously used by Casella and Moreno (2006) and Berger and
#' Molina (2005) and concern N = 178 measures of ozone concentration in the
#' atmosphere. Of the 10 main effects originally considered, we only make use
#' of those with an atmospheric meaning x4 to x10, as was done by Liang et al.
#' (2008). We then have 7 main effects which, jointly with the quadratic terms
#' and second order interactions, produce the above-mentioned p = 35 possible
#' regressors.
#'
#' @name Ozone35
#' @docType data
#' @format A data frame with 178 observations on the following 36 variables.
#' \describe{ \item{y}{Response = Daily maximum 1-hour-average ozone
#' reading (ppm) at Upland, CA} \item{x4}{500-millibar pressure height
#' (m) measured at Vandenberg AFB} \item{x5}{Wind speed (mph) at Los
#' Angeles International Airport (LAX)} \item{x6}{Humidity (percentage)
#' at LAX} \item{x7}{Temperature (Fahrenheit degrees) measured at
#' Sandburg, CA} \item{x8}{Inversion base height (feet) at LAX}
#' \item{x9}{Pressure gradient (mm Hg) from LAX to Daggett, CA}
#' \item{x10}{Visibility (miles) measured at LAX}
#' \item{x4.x4}{=x4*x4} \item{x4.x5}{=x4*x5}
#' \item{x4.x6}{=x4*x6} \item{x4.x7}{=x4*x7}
#' \item{x4.x8}{=x4*x8} \item{x4.x9}{=x4*x9}
#' \item{x4.x10}{=x4*x10} \item{x5.x5}{=x5*x5}
#' \item{x5.x6}{=x5*x6} \item{x5.x7}{=x5*x7}
#' \item{x5.x8}{=x5*x8} \item{x5.x9}{=x5*x9}
#' \item{x5.x10}{=x5*x10} \item{x6.x6}{=x6*x6}
#' \item{x6.x7}{=x6*x7} \item{x6.x8}{=x6*x8}
#' \item{x6.x9}{=x6*x9} \item{x6.x10}{=x6*x10}
#' \item{x7.x7}{=x7*x7} \item{x7.x8}{=x7*x8}
#' \item{x7.x9}{=x7*x9} \item{x7.x10}{=x7*x10}
#' \item{x8.x8}{=x8*x8} \item{x8.x9}{=x8*x9}
#' \item{x8.x10}{=x8*x10} \item{x9.x9}{=x9*x9}
#' \item{x9.x10}{=x9*x10} \item{x10.x10}{=x10*x10} }
#' @references Berger, J. and Molina, G. (2005)<DOI:j.1467-9574.2005.00275.x>
#' Posterior model probabilities via path-based pairwise priors. Statistica
#' Neerlandica, 59:3-15.
#'
#' Casella, G. and Moreno, E. (2006)<DOI:10.1198/016214505000000646> Objective
#' Bayesian variable selection. Journal of the American Statistical
#' Association, 101(473).
#'
#' Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association, 108: 340-352.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423.
#' @keywords datasets
#' @examples
#' data(Ozone35)
#'
"Ozone35"



#' SDM data
#'
#'The following data set contains 67 variables potentially related with Growth. The name of this dataset is related to its authors since it was firstly used in Sala i Martin, Doppelhofer and Miller (2004).
#'
#'
#' @name SDM
#' @docType data
#' @format A data frame with 88 observations on the following 68 variables
#' \describe{ \item{\code{y}}{Growth of GDP per capita at purchasing power parities between 1960 and 1996.}
#' \item{\code{ABSLATIT}}{Absolute latitude.}
#' \item{\code{AIRDIST}}{Logarithm of minimal distance (in km) from New York, Rotterdam, or Tokyo.}
#' \item{\code{AVELF}}{Average of five different indices of ethnolinguistic fractionalization which is the probability of two random people in a country not speaking the same language.}
#' \item{\code{BRIT}}{Dummy for former British colony after 1776.}
#' \item{\code{BUDDHA}}{Fraction of population Buddhist in 1960.}
#' \item{\code{CATH00}}{Fraction of population Catholic in 1960.}
#' \item{\code{CIV72}}{Index of civil liberties index in 1972.}
#' \item{\code{COLONY}}{Dummy for former colony.}
#' \item{\code{CONFUC}}{Fraction of population Confucian.}
#' \item{\code{DENS60}}{Population per area in 1960.}
#' \item{\code{DENS65C}}{Coastal (within 100 km of coastline) population per coastal area in 1965.}
#' \item{\code{DENS65I}}{Interior (more than 100 km from coastline) population per interior area in 1965.}
#' \item{\code{DPOP6090}}{Average growth rate of population between 1960 and 1990.}
#' \item{\code{EAST}}{Dummy for East Asian countries.}
#' \item{\code{ECORG}}{Degree Capitalism index.}
#' \item{\code{ENGFRAC}}{Fraction of population speaking English.}
#' \item{\code{EUROPE}}{Dummy for European economies.}
#' \item{\code{FERTLDC1}}{Fertility in 1960's.}
#' \item{\code{GDE1}}{Average share public expenditures on defense as fraction of GDP between 1960 and 1965.}
#' \item{\code{GDPCH60L}}{Logarithm of GDP per capita in 1960.}
#' \item{\code{GEEREC1}}{Average share public expenditures on education as fraction of GDP between 1960 and 1965.}
#' \item{\code{GGCFD3}}{Average share of expenditures on public investment as fraction of GDP between 1960 and 1965.}
#' \item{\code{GOVNOM1}}{Average share of nominal government spending to nominal GDP between 1960 and 1964.}
#' \item{\code{GOVSH61}}{Average share government spending to GDP between 1960 and 1964.}
#' \item{\code{GVR61}}{Share of expenditures on government consumption to GDP in 1961.}
#' \item{\code{H60}}{Enrollment rates in higher education.}
#' \item{\code{HERF00}}{Religion measure.}
#' \item{\code{HINDU00}}{Fraction of the population Hindu in 1960.}
#' \item{\code{IPRICE1}}{Average investment price level between 1960 and 1964 on purchasing power parity basis.}
#' \item{\code{LAAM}}{Dummy for Latin American countries.}
#' \item{\code{LANDAREA}}{Area in km.}
#' \item{\code{LANDLOCK}}{Dummy for landlocked countries.}
#' \item{\code{LHCPC}}{Log of hydrocarbon deposits in 1993.}
#' \item{\code{LIFE060}}{Life expectancy in 1960.}
#' \item{\code{LT100CR}}{Proportion of country's land area within 100 km of ocean or ocean-navigable river.}
#' \item{\code{MALFAL66}}{Index of malaria prevalence in 1966.}
#' \item{\code{MINING}}{Fraction of GDP in mining.}
#' \item{\code{MUSLIM00}}{Fraction of population Muslim in 1960.}
#' \item{\code{NEWSTATE}}{Timing of national independence measure: 0 if before 1914; 1 if between 1914 and 1945; 2 if between 1946 and 1989; and 3 if after 1989.}
#' \item{\code{OIL}}{Dummy for oil-producing country.}
#' \item{\code{OPENDEC1}}{Ratio of exports plus imports to GDP, averaged over 1965 to 1974.}
#' \item{\code{ORTH00}}{Fraction of population Orthodox in 1960.}
#' \item{\code{OTHFRAC}}{Fraction of population speaking foreign language.}
#' \item{\code{P60}}{Enrollment rate in primary education in 1960.}
#' \item{\code{PI6090}}{Average inflation rate between 1960 and 1990.}
#' \item{\code{SQPI6090}}{Square of average inflation rate between 1960 and 1990.}
#' \item{\code{PRIGHTS}}{Political rights index.}
#' \item{\code{POP1560}}{Fraction of population younger than 15 years in 1960.}
#' \item{\code{POP60}}{Population in 1960}
#' \item{\code{POP6560}}{Fraction of population older than 65 years in 1960.}
#' \item{\code{PRIEXP70}}{Fraction of primary exports in total exports in 1970.}
#' \item{\code{PROT00}}{Fraction of population Protestant in 1960.}
#' \item{\code{RERD}}{Real exchange rate distortions.}
#' \item{\code{REVCOUP}}{Number of revolutions and military coups.}
#' \item{\code{SAFRICA}}{Dummy for Sub-Saharan African countries.}
#' \item{\code{SCOUT}}{Measure of outward orientation.}
#' \item{\code{SIZE60}}{Logarithm of aggregate GDP in 1960.}
#' \item{\code{SOCIALIST}}{Dummy for countries under Socialist rule for considerable time during 1950 to 1995.}
#' \item{\code{SPAIN}}{Dummy variable for former Spanish colonies.}
#' \item{\code{TOT1DEC1}}{Growth of terms of trade in the 1960's.}
#' \item{\code{TOTIND}}{Terms of trade ranking}
#' \item{\code{TROPICAR}}{Proportion of country's land area within geographical tropics.}
#' \item{\code{TROPPOP}}{Proportion of country's population living in geographical tropics.}
#' \item{\code{WARTIME}}{Fraction of time spent in war between 1960 and 1990.}
#' \item{\code{WARTORN}}{Indicator for countries that participated in external war between 1960 and 1990.}
#' \item{\code{YRSOPEN}}{Number of years economy has been open between 1950 and 1994.}
#' \item{\code{ZTROPICS}}{Fraction tropical climate zone.}}
#' @references Sala i Martin, X., Doppelhofer, G., Miller, R.I. (2004)
#' <DOI: 10.1257/0002828042002570>.
#' Determinants of long-term growth: a Bayesian averaging of classical estimates (BACE) approach.
#'  American Economic Review 94: 813--835.
#' @keywords datasets
#' @examples
#' data(SDM)
#'
"SDM"


#' OBICE data
#'
#' Dataset corresponding to the OBICE study (Zurriaga et al 2011) where factors associated with childhood obesity
#' are studied. The data were collected in 2007 and 2008 through several questionnaries and 
#' n=1188 children were enroled in the study. It contains 155 variables.
#' This is a case and control study with 437 cases (obese) and 751 controls (not obese).
#' Purposedly the dataset is distributed without any post-processing
#' hence, many variables may contain unavailable observations coded in different way.
#'
#'
#' @name OBICE
#' @docType data
#' @format A data frame with 1188 entries and 121 variables. The more relevant are described. Contact us if you need specific information of any other.
#' \describe{
#' \item{Acostarse}{:does he/she eat before going to bed? (yes/no) }
#' \item{ActFisica}{:(physic activity) factor coded as 1 None; 2 less than monthly; 3 less than weekly; 4 less than 2/week; 5 at least 2/week }
#' \item{ActivDepor}{:weekly hours devoted to sports activity }
#' \item{Almuerzo}{: }
#' \item{Bebida}{:(main dring accompanying the main meal) 1 water tap; 2 bottle water; 3 soda; 4 natural juices; 5 bottle juices; 6 Milk (and derivatives); 7 Other }
#' \item{Caso01}{: }
#' \item{Cena}{: }
#' \item{Chuches}{:Sweets and soft drinks weekly consumption (how many times) }
#' \item{CincoComidas}{:does he/she have regularly 5 meals per day? (0 is No; 1 is Yes) }
#' \item{clSocEl}{: }
#' \item{clSocElla}{: }
#' \item{clSocXiquet}{: }
#' \item{ComedorEsc}{: }
#' \item{Comida}{: }
#' \item{Daceite}{: }
#' \item{Dcereal}{: }
#' \item{Desayuno}{: }
#' \item{Descubrimiento}{: }
#' \item{Dgalleta}{: }
#' \item{Dislipemias}{: }
#' \item{DislipeRelacion}{: }
#' \item{Dleche}{: }
#' \item{Dotros}{: }
#' \item{Dpan}{: }
#' \item{Dzumoenv}{: }
#' \item{Dzumonat}{: }
#' \item{Edad}{:years old }
#' \item{EntreHoras}{: }
#' \item{EstudiosMadre}{: }
#' \item{EstudiosMadreSinCon}{: }
#' \item{EstudiosPadre}{: }
#' \item{EstudiosPadreSinCon}{: }
#' \item{Faperitivos}{: }
#' \item{Faperitivosmp}{: }
#' \item{Farroz}{: }
#' \item{Farrozmp}{: }
#' \item{Fcarnes}{: }
#' \item{Fcarnesmp}{: }
#' \item{Fchucherias}{: }
#' \item{Fchucheriasmp}{: }
#' \item{Fdulces}{: }
#' \item{Fdulcesmp}{: }
#' \item{Ffiambres}{: }
#' \item{Ffiambresmp}{: }
#' \item{Ffritos}{: }
#' \item{Ffritosmp}{: }
#' \item{Ffruta}{: }
#' \item{Ffrutamp}{: }
#' \item{Fhuevos}{: }
#' \item{Fhuevosmp}{: }
#' \item{Flacteos}{: }
#' \item{Flacteosmp}{: }
#' \item{Flegumbres}{: }
#' \item{Flegumbresmp}{: }
#' \item{Fpan}{: }
#' \item{Fpanmp}{: }
#' \item{Fpescado}{: }
#' \item{Fpescadomp}{: }
#' \item{Fprecocina}{: }
#' \item{Fprecocinamp}{: }
#' \item{Frefrescos}{: }
#' \item{Frefrescosmp}{: }
#' \item{Fruta}{:usual consumption of fruit? (0 is No; 1 is Yes) }
#' \item{FrutaVerdura}{: }
#' \item{Fverduras}{: }
#' \item{Fverdurasmp}{: }
#' \item{HorasPantDia}{: }
#' \item{HorasPCDiaPond}{:daily hours playing videogames and/or in internet (weekends included) }
#' \item{HorasPCsem1}{: }
#' \item{HorasPCsem2}{: }
#' \item{HorasTV}{: }
#' \item{HorasTVDiaPond}{:daily hours watching TV (weekends included) }
#' \item{HorasTVsem1}{: }
#' \item{HorasTVsem2}{: }
#' \item{HoraSuenyo}{:daily hours sleeping }
#' \item{HTA}{: }
#' \item{HTARelacion}{: }
#' \item{IMC}{: }
#' \item{IndEdadComedorEscolar}{: }
#' \item{IntolGlucosa}{: }
#' \item{IntolRelacion}{: }
#' \item{LactMater}{: }
#' \item{LactMaterna}{: breast-feeding (1 is Yes; 0 is No) }
#' \item{LactMatMeses}{: }
#' \item{LactMatSemanas}{: }
#' \item{MadreObesa}{: }
#' \item{MadreObesa01}{:is the mother obese? (0 is No; 1 is Yes) }
#' \item{Merienda}{:Afternoon snack  (1 is Yes; 0 is No) }
#' \item{NumComidas}{: }
#' \item{NumContOK}{: }
#' \item{NumControles}{: }
#' \item{NumHnos}{: }
#' \item{NumHnosOb}{: }
#' \item{NumPadresEsp02}{: }
#' \item{NumPadresObesos}{: }
#' \item{OrdenadorDiario}{: }
#' \item{OrdenadorFinDe}{: }
#' \item{OsteoRelacion}{: }
#' \item{OtrosPatol}{: }
#' \item{PadreObeso}{:is the father obese? (0 is No; 1 is Yes) }
#' \item{PesoActual}{:current weight (in kilograms) }
#' \item{PesoNac}{:weight born (in grams) }
#' \item{PorcHnosObesos}{: }
#' \item{porcHnosObesosOK}{: }
#' \item{Postre}{: }
#' \item{ProbOsteo}{: }
#' \item{ProbPsico}{: }
#' \item{ProbResp}{: }
#' \item{PsicoRelacion}{: }
#' \item{ResoponF01}{: }
#' \item{RespRelacion}{: }
#' \item{Semlact}{: }
#' \item{Sexo}{:female (1); male (0) }
#' \item{TallaAct}{:current height (in meters) }
#' \item{TallaNac}{:height born (in centimeters) }
#' \item{Tipocaso}{: }
#' \item{Tipocaso.y}{: }
#' \item{TipoObeso}{: }
#' \item{TVDiario}{: }
#' \item{TVFinSemana}{: }
#' \item{Verduras}{:usual consumption of vegetables? (0 is No; 1 is Yes) }
#' }
#' @references Zurriaga, O., Perez-Panades, J., Quiles, J. , Gil, M.,
#' Anes, Y., Quiñones, C., Margolles, M., Lopez-Maside, A.,
#' Vega-Alonso, A., Miralles M. and Recent OBICE Research Group (2011)
#' Factors associated with childhood obesity in Spain. The OBICE
#' study: a case–control study based on sentinel networks. Public Health Nutrition
#' 14(6), 1105–1113.
#' @keywords datasets
#' @examples
#' data(OBICE)
#'
"OBICE"
