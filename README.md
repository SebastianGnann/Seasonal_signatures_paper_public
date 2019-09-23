# Seasonal_signatures_paper_public

Please email sg16200@bristol.ac.uk if you have any questions.

Matlab version: **Matlab R2018a**

- Repository containing code to do the analysis and results used to create (most of) the plots for the seasonal signatures paper (https://www.hydrol-earth-syst-sci-discuss.net/hess-2019-463/): sine-fitting, calculation of seasonal signatures, code to call MARRMoT, etc.
- Note that the original data is not uploaded here (see paper for data sources), but results that can be used to create (some of) the plots. This includes the results of the modelling exercise using MARRMoT and the calculated seasonal signatures. We also added examples with synthetic catchment data to show how the code can be used with any data of the correct format. See below for a description of the main scripts.
- There are a few additional open-source packages needed: the Brewerplot package (https://github.com/DrosteEffect/BrewerMap) for the colour schemes and the MARRMoT toolbox (https://zenodo.org/record/3235664#.XV-mKehKiUk) for the modelling exercise.
- For the modelling results presented in the paper we used the Latin Hypercube sampling function implemented in the SAFE toolbox (https://www.safetoolbox.info/). The example file now uses the Matlab function.

#### Additional information

> **Example_calculateSeasonalSignatures.m**: calculates seasonal signatures (amplitude ratio vs. phase shift) using example data.

> **Example_createTimeSeries.m**: creates the example time series used in Example_calculateSeasonalSignatures.m.

> **Example_makePlotsSeasonalSignatures.m**: plots seasonal signatures (amplitude ratio vs. phase shift) calculated using example data (similar to the plots in the paper).

> **makePlotsSeasonalSignatures.m**: plots seasonal signatures (amplitude ratio vs. phase shift) calculated using UK and US data. These are the results used to make the plots in the paper, yet not all catchment attributes are available here. Please go to the original data sources (described in paper) to obtain the data.

> **makePlotsTheoreticalSeasonalSignatures.m**: plots theoretical seasonal signatures (amplitude ratio vs. phase shift) for different linear reservoir configurations (theory described in paper).

> **MARRMoT_evaluateResults.m**: plots the results of the Monte Carlo runs using MARRMoT (results are stored in Data_and_results/Results_MARRMoT).

> **MARRMoT_evaluateRobustness.m**: checks how robust the results of the Monte Carlo runs with different sample sizes are.

> **MARRMoT_runExample.m**: shows how you can run MARRMoT with Latin Hypercube sampling using the example time series.
