# A benchmarking study of distance-based clustering methods for mixed-type data

In this repository, you can find the R code used for the project "A benchmarking study of distance-based clustering methods for mixed-type data". The project has been submitted at Imperial College London, in partial fulfilment of the requirements for the MSc in Statistics.

* Author: Efthymios Costa
* Supervisors:
  * Dr. Ioanna Papatsouma (Imperial College London)
  * Dr. Angelos Markos (Democritus University of Thrace)

## Description of folders, files & R scripts

* `clustermasking.R`: R code used in Section 4.1, illustrating the issue of Cluster Masking.
* __Folder KDDCup199__:
  * `kddcup99.RData`: Original KDD Cup 1999 data set.
  * `kddcup99sample.RData`: 10% sample of the KDD Cup 1999 used in the simulation study.
  * _Subfolder SimulationStudy_:
    * `KDDCup1999SimulationStudy.R`: R code used in Chapter 5 for the purpose of the simulation study on the real world data set.
    * `kddcup99aris.RData`: ARI values obtained from the simulation study.
  * _Subfolder Visualisations_:
    * `kmedres_kddcup9.RData`: Random cluster allocation obtained using Mixed K-Means.
    * `msRes_kddcup9.Rdata`: Random cluster allocation obtained using Modha-Spangler K-Means.
    * `nlpcakddcup9kmeans.RData`: Random cluster allocation obtained using NLPCA & K-Means.
    * `outfmix_kddcup9.RData`: Random cluster allocation obtained using Mixed Factorial K-Means.
    * `outk_kddcup9.RData`: Random cluster allocation obtained using K-Prototypes.
    * `outkm_kddcup9.RData`: Random cluster allocation obtained using FAMD & K-Means.
    * `outmix_kddcup9.RData`: Random cluster allocation obtained using Mixed Reduced K-Means.
    * `pam_fit_kddcup9.RData`: Random cluster allocation obtained using Gower's dissimilarity & PAM.
    * `tSNEVisualisationKDDCup1999.R`: R code used in Section 5.4 for visualisation using t-SNE. Code for 2D and 3D plots of the data set and cluster allocations obtained using each of the 8 clutering algorithms considered is included.
* __Folder ArtificialData__:
  * `discretise.R`: R function that has been used for discretising continuous attributes.
  * _Subfolder CDRSimulationStudy_:
    * `CDRSimStudyCatContRatio.R`: R code used for investigating the effect of categorical to continuous attributes ratio and of the kappa coefficient in Mixed CDR.
    * `CDRSimStudyCatLevels.R`: R code used for investigating the effect of the number of categorical levels used in discretisation and of the kappa coefficient in Mixed CDR.
    * `CDRSimStudyNClust.R`: R code used for investigating the effect of the number of clusters and of the kappa coefficient in Mixed CDR.
    * `CDRSimStudyOverlap.R`: R code used for investigating the effect of the level of cluster overlap and of the kappa coefficient in Mixed CDR.
    * `kappasmat.RData`: Simulation results of `CDRSimStudyCatContRatio` for 1 data set.
    * `kappasmat2.RData`: Simulation results of `CDRSimStudyCatContRatio` for 4 data sets.
    * `kappasmatlevels.RData`: Simulation results of `CDRSimStudyCatLevels` for 1 data set.
    * `kappasmatlevels2.RData`: Simulation results of `CDRSimStudyCatLevels` for 4 data sets.
    * `kappasmatnclust.RData`: Simulation results of `CDRSimStudyNClust` for 1 data set.
    * `kappasmatnclust2.RData`: Simulation results of `CDRSimStudyNClust` for 4 data sets.
    * `kappasmatoverlap.RData`: Simulation results of `CDRSimStudyOverlap` for 1 data set.
    * `kappasmatoverlap2.RData`: Simulation results of `CDRSimStudyOverlap` for 4 data sets.
  * _Subfolder TechniquesSimulationStudy_:
    * `SimStudyCatContRatio.R`: R code used for investigating the effect of categorical to continuous attributes ratio in the performance of the clustering algorithms considered.
    * `SimStudyCatLevels.R`: R code used for investigating the effect of the number of categorical levels used in discretisation in the performance of the clustering algorithms considered.
    * `SimStudyNClust.R`: R code used for investigating the effect of the number of clusters in the performance of the clustering algorithms considered.
    * `SimStudyOverlap.R`: R code used for investigating the effect of the level of cluster overlap in the performance of the clustering algorithms considered.
    * `mixsimresmat.RData`: Simulation results of `SimStudyCatContRatio` for 1 data set.
    * `mixsimresmat2.RData`: Simulation results of `SimStudyCatContRatio` for 4 data sets.
    * `mixsimresmatlevels.RData`: Simulation results of `SimStudyCatLevels` for 1 data set.
    * `mixsimresmatlevels2.RData`: Simulation results of `SimStudyCatLevels` for 4 data sets.
    * `mixsimresmatnclust.RData`: Simulation results of `SimStudyNClust` for 1 data set.
    * `mixsimresmatnclust2.RData`: Simulation results of `SimStudyNClust` for 4 data sets.
    * `mixsimresmatoverlap.RData`: Simulation results of `SimStudyOverlap` for 1 data set.
    * `mixsimresmatoverlap2.RData`: Simulation results of `SimStudyOverlap` for 4 data sets.
