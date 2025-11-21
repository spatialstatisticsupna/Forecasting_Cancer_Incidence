
# Multivariate Spatio-temporal Modelling for Completing Cancer Registries and Forecasting Incidence

This repository contains the R code to fit in INLA the code to replicate and reproduce the model validation results of the paper entitled "Multivariate Spatio-temporal Modelling for Completing Cancer Registries and Forecasting Incidence".

## Table of contents

1.  [Data](#Data)
2.  [R code](#Rcode)
3.  [Acknowledgements](#Acknowledgements)
4.  [References](#Ref)

# Data <a name="Data"/> {#data}

This folder contains the datasets and cartography files used to validate the models presented in the work.

-   `Lung_Male_Aggregate.Rdata`: This dataset includes complete counts for lung cancer and the corresponding population at risk.

-   `carto_england.shp`: Cartography of England.

-   `Lung_Male_Missing.Rdata`: This dataset contains counts for lung cancer with missing values, along with the complete corresponding population at risk. It is used to fit and validate the models presented in the work. 

-   `Information_missing_areas.Rdata`: This dataset includes the information for the selected missing areas, which is required to produce the figures presented in the paper.

# R code <a name="Rcode"/> {#Rcode}

This folder contains the R code to replicate and reproduce the model validation results presented in the paper.

-   `3-year-ahead_Code_model1.R`, `3-year-ahead_Code_model2.R` and `3-year-ahead_Code_model3.R`: Code to fit the models presented in the paper.

-   `Results_Missing_BestModels.R`: Code used to produce the tables and figures in section 4.2.1 of the paper.

-   `Results_Projections_BestModels.R`:  Code used to produce the figures in section 4.2.2 of the paper.

-   `Results_National_BestModels.R`: Code used to produce the tables in section 5 of the paper.


# Acknowledgements <a name="Acknowledgements"/> {#acknowledgements}

The work was supported by Project PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033, Ayudas Predoctorales Santander UPNA 2021-2022 and BIOSTATNET - PROYECTOS REDES DE INVESTIGACIÓN 2024
- RED2024-153680-T/MICIU/AEI/. ![plot](https://github.com/spatialstatisticsupna/Prior_Smoothing/blob/main/micin-aei.jpg)

# References <a name="Ref"/>

