# Polygenic-Risk-Scores-on-Mice

Background

Polygenic risk scores (PRSs) are metrics which evaluate a person’s individual risk to a trait or disease. In recent years, due to advancement of technology and data collection, PRSs have gained a lot of attention as a way to target therapies or prevention programmes to individuals. As of such statistical methods and tools for calculating PRSs are becoming increasingly important as scientists push for the application of PRSs in healthcare systems.

The basis of PRSs remains the availability of good genetic data. Unfortunately, in recent years, genetic data has been labelled by many as bias, with most genetic research being conducted on a select few populations-mainly white populations from European descent.  As the calculation of PRSs is based on genetic structures, having bias genetic datasets that exclude minority populations could be problematic and lead to the widening of health disparities.

Objectives

The primary objective of this research project was to better understand and visualize how population structures affect the calculation of genetic risk scores such as PRSs using mice as models. Secondary objectives looked at understanding how accurate PRS models were and to identify methods in the literature that demonstrate how population structures can be corrected for in genetic data.

Methods

PRS scores for a single trait in mice were calculated using two penalized regression techniques – LASSO and Elastic Net. Scores were stratified by a proxy indicator for population (litter) to understand the effect of population structure on the distribution of scores in mice. The two different PRS models developed were all trained and tested on different populations to give a measure of accurate the PRS models were.


Results

Results from the research project find that PRS scores for mice in this dataset did differ when stratified by litter. PRS models for both methods of regularization -LASSO and Elastic Net were the most accurate when they trained and tested on the same litter (population of mice) suggesting the PRS for one population have limited generalizability. 


Conclusions 

Findings from this research project on mouse models are relevant to human populations as well. This research project sheds light on the need to ensure that genetic research encompasses minority groups and that genotype information is collected for all groups before genetic risks scores are calculated. The research project also highlights the need to better evaluate the accuracy of PRS models before they are applied to make any clinical or policy decisions.
