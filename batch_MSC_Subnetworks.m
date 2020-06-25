MSCnums = [1:10];
xdist = 30;
thresholds = [.0001 .0002 .0005 .001 .002 .005 .01 .02 .05];


%% Analyses to run

badverts = ft_read_cifti_mod('/data/nil-bluearc/GMT/Evan/Scripts/scripts_for_subnetworksppr/MSCmean_badverts_cort.dtseries.nii');
badverts = find(badverts.data);
    
for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    
    infomap_outfolder = ['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' MSCname '_infomap_subnetworks_subcortreg_ignoreverts/'];
    
    
    
    
    mkdir(infomap_outfolder);
    cd(infomap_outfolder);
    
    
    
    corrmatname = ['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegressedDconns/' MSCname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dconn.nii'];
    
    dmatname = ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/' MSCname '/fsaverage_LR32k/' MSCname '_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];

    Run_Infomap_nodespecificthresh(corrmatname, dmatname, xdist, thresholds, 0, infomap_outfolder, badverts,100,8,false);%, structure_indices);
    
    communities = modify_clrfile('simplify','rawassn.txt',10);
    regularized = regularize(communities);
    data = ft_read_cifti_mod(corrmatname);
    data.data = [];
    data.dimord = 'pos_time';
    data.data = regularized;
    ft_write_cifti_mod([MSCname '_rawassn_minsize10_regularized'],data);
    
    
end