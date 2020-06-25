MSCnums = [1:10];
matches = {'DMN'};
col = 4;
xdist = 30;
density = .05;
cortsizethresh = 20;
subcortsizethresh = 10;


othernetworks_toinclude = [3 7]; %FP and Language

DMN_subnetwork_IDnums = [1 7.5 16 11.5 10.8 1.5 3.75 15.5 10.5];



for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    disp(MSCname)
    
    DMN_subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' MSCname '_infomap_subcortreg_ignoreverts/' MSCname '_rawassn_minsize10_regularized_DMNmatch_v2_recolor.dtseries.nii']);
    thissub_networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V2/infomap/' MSCname '/' MSCname '_rawassn_minsize400_regularized_recolored_cleaned.dscalar.nii']);
    subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' MSCname '_infomap_subcortreg_ignoreverts/' MSCname '_rawassn_minsize10_regularized.dtseries.nii']);
    subnetworks.data = subnetworks.data(:,col);
    subnetwork_IDs = unique(subnetworks.data); subnetwork_IDs(subnetwork_IDs<1) = [];
    
    allclusters = [];
    prevstring = [];
    outputclusters = cell(length(subnetwork_IDs),1);
    counts = zeros(length(subnetwork_IDs),1);
    parfor IDnum_parfor = 1:length(subnetwork_IDs)
        string = ['clustering network ' num2str(IDnum_parfor) ' of ' num2str(length(subnetwork_IDs))];
        disp(string);
        ID = subnetwork_IDs(IDnum_parfor);
        
        outputclusters{IDnum_parfor} = cifti_cluster_surfacevolume(subnetworks,ID-.001,ID+.001,cortsizethresh,subcortsizethresh);
        counts(IDnum_parfor) = size(outputclusters{IDnum_parfor},2);
    end
    cum_num_total = [0 ; cumsum(counts)];
    allclusters = zeros(size(outputclusters{end},1),cum_num_total(end));
    for IDnum = 1:length(subnetwork_IDs)
        if cum_num_total(IDnum) < cum_num_total(IDnum+1)
            allclusters(:,(cum_num_total(IDnum)+1):cum_num_total(IDnum+1)) = outputclusters{IDnum};
            
        end
    end
    clear outputclusters

    data = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegTimeSeries/' MSCname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dtseries.nii']);
    
    distances = smartload(['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/' MSCname '/fsaverage_LR32k/' MSCname '_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat']);
    
    cluster_tcs = zeros(size(allclusters,2),size(data.data,2));
    IDs_forPC{MSCnum} = zeros(size(allclusters,2),1);
    subnetworkornot_IDs{MSCnum} = zeros(size(allclusters,2),1);
    centroids = zeros(size(allclusters,2),1);
    use = false(size(allclusters,2),1);
    cort = false(size(allclusters,2),1);
    for IDnum = 1:size(allclusters,2)
        clusterinds = find(allclusters(:,IDnum));
        cluster_tcs(IDnum,:) = mean(data.data(clusterinds,:),1);
        [~,centroidind] = min(sum(distances(clusterinds,clusterinds),2));
        centroids(IDnum) = clusterinds(centroidind);
        cort(IDnum) = all(clusterinds<=59412);

        if mode(DMN_subnetworks.data(clusterinds)) ~= 0
            IDs_forPC{MSCnum}(IDnum) = 1;
            subnetworkornot_IDs{MSCnum}(IDnum) = mode(DMN_subnetworks.data(clusterinds));
            use(IDnum) = true;
        elseif cort(IDnum) && any(mode(thissub_networks.data(clusterinds))==othernetworks_toinclude)
            IDs_forPC{MSCnum}(IDnum) = mode(thissub_networks.data(clusterinds));
            subnetworkornot_IDs{MSCnum}(IDnum) = 0;
            use(IDnum) = true;
        end
        
    end
    cluster_tcs = cluster_tcs(use,:);
    centroids = centroids(use);
    IDs_forPC{MSCnum} = IDs_forPC{MSCnum}(use);
    subnetworkornot_IDs{MSCnum} = subnetworkornot_IDs{MSCnum}(use);
    cort = cort(use);
    
    FCmat{MSCnum} = paircorr_mod(cluster_tcs');
    
    dmat{MSCnum} = distances(centroids,centroids);
    
    FCmat{MSCnum}(dmat{MSCnum}<xdist) = 0;
    FCmat{MSCnum}(~cort,~cort) = 0;
    
end
    

%%
subnetwork_PCs_bynet = ones(length(MSCnums),length(DMN_subnetwork_IDnums)) .* NaN;
subnetwork_WMdegZ_bynet = ones(length(MSCnums),length(DMN_subnetwork_IDnums)) .* NaN;
for MSCnum = MSCnums
    sorted = sort(FCmat{MSCnum}(:),'descend');
    thresh = sorted(ceil(length(sorted) .* density));
    connectionmat = FCmat{MSCnum} .* single(FCmat{MSCnum} > thresh);
    

    mean_PCs = PC_calc_simple(connectionmat,IDs_forPC{MSCnum});
    within_mod_degZ=module_degree_zscore(connectionmat,IDs_forPC{MSCnum});

    for n = 1:length(DMN_subnetwork_IDnums)
        if any(abs(subnetworkornot_IDs{MSCnum}-DMN_subnetwork_IDnums(n))<.01)
            subnetwork_PCs_bynet(MSCnum,n) = nanmean(mean_PCs(abs(subnetworkornot_IDs{MSCnum}-DMN_subnetwork_IDnums(n))<.01));
            subnetwork_WMdegZ_bynet(MSCnum,n) = nanmean(within_mod_degZ(abs(subnetworkornot_IDs{MSCnum}-DMN_subnetwork_IDnums(n))<.01));
        end
    end
    
    
end
order = [1 8 2 3 9 4 6 5 7];

make_network_bargraph(DMN_subnetwork_IDnums(order),nanmean(subnetwork_PCs_bynet(:,order),1),nanstd(subnetwork_PCs_bynet(:,order),1)./ sqrt(10),false)

make_network_bargraph(DMN_subnetwork_IDnums(order),nanmean(subnetwork_WMdegZ_bynet(:,order),1),nanstd(subnetwork_WMdegZ_bynet(:,order),1)./ sqrt(10),false)

%%
ps = zeros(4,2);
ts = zeros(4,2);
for i = 1:4
    [H,P,CI,STATS] = ttest(subnetwork_PCs_bynet(:,i+3),subnetwork_PCs_bynet(:,1));
    ps(i,1) = P;
    ts(i,1) = STATS.tstat;
    [H,P,CI,STATS] = ttest(subnetwork_WMdegZ_bynet(:,1),subnetwork_WMdegZ_bynet(:,i+3));
    ps(i,2) = P;
    ts(i,2) = STATS.tstat;
end
p_fdr = FDR(ps(:),.05);
disp(p_fdr)
disp(ps)
disp(ts)