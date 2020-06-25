Code to run analyses of MSC data from Gordon et al., 2020, PNAS

Make sure code from https://github.com/MidnightScanClub/MSCcodebase/ is in your path (especially the Utilities/ folder)

1) Regress local subcortical signal using subcort_regression.m
2) Run subnetwork detection at multiple thresholds using batch_MSC_Subnetworks.m
3) Test within-subnetwork homogeneity at each threshold vs random rotations using Subnetworks_pcahomogeneity_vsrandom.m (first run make_rotations.m to generate the random rotations)
4) Test variance explained in task activation maps vs random rotations using Subnetworks_activation_Rsq_vsrandom.m
5) Manually match subnetworks at an optimal threshold across individuals, as in the manuscript
6) Calculate matched subnetwork FC to large-scale networks using MSC_subnetwork_FC_tonetworks.m
7) Calculate matched subnetwork graph metrics using MSC_matched_subnetwork_PCandDeg.m
8) Calculate matched subnetwork time delays using MSC_subnetworks_tdmx.m (must have code from https://github.com/RaichleLab/lag-code/ in your path)

This code is provided "as is", without any warranty. We cannot provide technical support for this code. For general questions, please contact Evan Gordon (egordon@wustl.edu).