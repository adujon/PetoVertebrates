# PetoVertebrates

This repository contains the codes and raw datasets used to perform the analyses in the following manuscript: 

Towards more rigorous testing of Peto’s paradox: a Bayesian re-analysis and discussion of the effect of body mass and longevity on cancer risk across vertebrates.
Antoine M. Dujon1,2*, Peter A. Biro1, Beata Ujvari1, Frédéric Thomas2.

1Deakin University, Geelong, School of Life and Environmental Sciences, Centre for Integrative Ecology, Waurn Ponds, Vic 3216, Australia. 

2CREEC/(CREES), MIVEGEC, Unité Mixte de Recherches, IRD 224–CNRS 5290–Université de Montpellier, Montpellier, France

*Corresponding author: antoine.dujon@yahoo.fr 


-------------------------------------------------------------------------------------------------------
The rawdata folder contains the cancer prevalence raw, body mass and maximum longevities datasets included in the manuscript

The phylogenies folder contains the phylogenetic trees used in the Bayesian phylogenetic models

The Regressions folder contains the R scripts used to fit the model for each taxa and study (code are fully annotated)

The Shuffling Simulations folder contains the R scripts used to run the simulation investigating the effect of variation in the number of necropsies on effect size (here slope between the log10 BodyMass/Longevity and prevalence)

-------------------------------------------------------------------------------------------------------
Raw cancer prevalence data was sourced from the following publications

Kapsetaki, S. E., Z. T. Compton, J. Dolan, V. Κ. Harris, W. Mellon, S. M. Rupp, E. G. Duke, T. M. Harrison, S. Aksoy, M. Giraudeau, O. Vincze, K. J. McGraw, A. Aktipis, M. Tollis, A. Μ. Boddy, and C. C. Maley. 2024. Life history traits and cancer prevalence in birds. Evol Med Public Health, doi: 10.1093/emph/eoae011

Boddy, A. M., L. M. Abegglen, A. P. Pessier, A. Aktipis, J. D. Schiffman, C. C. Maley, and C. Witte. 2020. Lifetime cancer prevalence and life history traits in mammals. Evol Med Public Health 2020:187–195. https://doi.org/10.1093/emph/eoaa015

Vincze, O., F. Colchero, J.-F. Lemaître, D. A. Conde, S. Pavard, M. Bieuville, A. O. Urrutia, B. Ujvari, A. M. Boddy, C. C. Maley, F. Thomas, and M. Giraudeau. 2022. Cancer risk across mammals. Nature 601:263–267 https://doi.org/10.1038/s41586-021-04224-5

Compton, Z. T., W. Mellon, V. K. Harris, S. Rupp, D. Mallo, S. E. Kapsetaki, M. Wilmot, R. Kennington, K. Noble, C. Baciu, L. N. Ramirez, A. Peraza, B. Martins, S. Sudhakar, S. Aksoy, G. Furukawa, O. Vincze, M. Giraudeau, E. G. Duke, S. Spiro, E. Flach, H. Davidson, C. I. Li, A. Zehnder, T. A. Graham, B. V. Troan, T. M. Harrison, M. Tollis, J. D. Schiffman, C. A. Aktipis, L. M. Abegglen, C. C. Maley, and A. M. Boddy. 2024. Cancer Prevalence across Vertebrates. Cancer Discovery, https://doi.org/10.1158/2159-8290.CD-24-0573

----------------------------------------------------------------------------------------------------
Phylogenies were obtained from the following publications

Upham, N. S., J. A. Esselstyn, and W. Jetz. 2019. Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS Biol 17:e3000494.

Tonini, J. F. R., K. H. Beard, R. B. Ferreira, W. Jetz, and R. A. Pyron. 2016. Fully-sampled phylogenies of squamates reveal evolutionary patterns in threat status. Biol Conserv 204:23–31. Elsevier Ltd.

----------------------------------------------------------------------------------------------------
Body size and longevity datasets were obtained from the following publications / databases and from Kaptsetaki, Compton, Vinzce and Boddy's papers

https://genomics.senescence.info

Tacutu, R., T. Craig, A. Budovsky, D. Wuttke, G. Lehmann, D. Taranukha, J. Costa, V. E. Fraifeld, and J. P. De Magalhães. 2013. Human Ageing Genomic Resources: Integrated databases and tools for the biology and genetics of ageing. Nucleic Acids Res 41.


