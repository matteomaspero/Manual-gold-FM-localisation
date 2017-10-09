# Welcome to the code for analyisis of intra-prostatic gold Fiducial Marker (FM) MR-based manual localisation

This work aims at present the analysis related to the article recently submitted to Radiation Oncology by Maspero Matteo and co-workers (2017) about intra-prostatic FM localisation on sole MR images.
Note that the code will be published upon publication of the related manuscript.

## What can you find here?

The matlab code to analyse the five observer FM localisation in terms of precision and accuracy.
* **Precision** for each FM, an average position is found among the five observers and the distance to the average position of each observer is calculated. The precision is calcualted as 95% Limit of Agreement (LoA) for the single FM and for the Centre of Mass (CM). If LoA > 2 mm in one of the three plane, the FM is considered as not precisely located.
* **Accuracy** against Computed Tomography-based localisation, which is considered the gold standard. The Inter-marker Distance (ID) among corresponding fiducial markers is calculated and the absolute distance of ID on CT - ID on MR for corresponding markers is considered as accuracy. The distribution of the ID is characterised in terms of descriptive statistics (mean, standard deviation, median, minimum and maximum). The imprecisely locate FM have been excluded from the statistics.

## How to...

Run the code? Simply clone/download the repository and run within Matlab LoA_LoACM.m to analyse the precision and Ibter_Marker_Distaces.m to analyse accuracy. Note that the code has been developed and benchmarked in Matlab R2015a and we do not ensure/support workability for different versions of Matlab.

## Disclaimer

The code is released under a Commopn Development and Distribution License, under the Attribution-ShareAlike 4.0 International clausole. In practice, we ask you to cite/refer the original paper and to share your own work released under a similar license.
The code is not intended for medical use, therefore we do not take any responsability for the misuse of the code.

## Support or Contact

In case of questions or issues, please contact the owner of the repository & we’ll help you sort it out. Thank yo very much in advance for your questions!

## Authors and Contributors

Maspero Matteo ([@matteomaspero](https://github.com/matteomaspero)), PhD candidate at University Medical Center Utrecht, The Netherlands. In case you are using any of the material, please refer to the article. The link to the article will be here presented as soon as published.

_Are you interested in related code from the same author?_
Check out the following github project: [https://matteomaspero.github.io/pseudo-CT_generation/](https://matteomaspero.github.io/pseudo-CT_generation/)
