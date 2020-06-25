%% Main lag script used to compute time delay matrix and lag projection
%
%
%NOTE: must have all scripts from https://github.com/RaichleLab/lag-code/ in your path



lag_lim = 4;    % lag limit (in seconds)
lags = -3:3;    % range of TR shifts; max(lags) = round(lag_lim/tr + 1)
tr = 2.2;

MSCnums = 1:10;

colors = [1 7.5 16 11.5 10.8 1.5 3.75 15.5 10.5];%[1 7.5 16 10.5 15.5 1.5 3.75 11.5 2.5];
order = [1 8 2 3 9 4 6 5 7];
colors = colors(order);

all_lagmat = nan(length(colors),length(colors),max(MSCnums));

for MSCnum = MSCnums
    
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    data = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegTimeSeries/' MSCname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dtseries.nii']);
    thissub_subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' MSCname '_infomap_subcortreg_ignoreverts/' MSCname '_rawassn_minsize10_regularized_DMNmatch_v2_recolor.dtseries.nii']);
    thissub_subnetworks = thissub_subnetworks.data;
    tmasklist = ['/data/nil-bluearc/GMT/Evan/MSC/subjects/' MSCname '_TMASKLIST_V2.txt'];
    [subjects,tmasks] = textread(tmasklist,'%s%s');
    
    subnetwork_IDs = unique(thissub_subnetworks); subnetwork_IDs(subnetwork_IDs<1) = [];
    
    
    num_nodes = length(colors);
    
    sessions = [];
    tmask_concat = [];
    inds_withintmaskedconcat = [];
    prevind = 0;
    for s = 1:numel(subjects)
        tmask = logical(load(tmasks{s}));
        sessions = [sessions ; repmat(s,size(tmask))];
        tmask_concat = [tmask_concat ; tmask];
        inds_withintmasked = zeros(length(tmask),1);
        inds_withintmasked(tmask) = [(prevind+1) : (prevind + nnz(tmask))];
        inds_withintmaskedconcat = [inds_withintmaskedconcat ; inds_withintmasked];
        prevind = prevind + nnz(tmask);
    end
    
    
    
    %% Setup
    % Set parameters
    
    %mkdir(outdir)
    
    
    % Specify data parameters
    
    
    motion_thresh = .2;    % important: must match motion criteria used during preproc
    
    min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)
    
    %% Loop over subjects
    % initialize group matrices for running sum
    grp_lags = single(nan(num_nodes));    % peak lags
    grp_ZL = grp_lags;      % zero-lag correlation
    grp_peak = grp_lags;    % peak correlation
    
    grp_lags_nans = single(zeros(num_nodes));
    grp_ZL_nans = grp_lags_nans;
    grp_peak_nans = grp_lags_nans;
    
    for s = 1:numel(subjects)
        tic
        subj = subjects{s};
        disp(['Processing ' subj]);
        
        % initialize subject matrices
        subj_lags = single(nan(num_nodes)); % peak lags
        subj_ZL = subj_lags;   % zero-lag correlation
        subj_peak = subj_lags; % peak correlation (correlation at optimal subnetwork_IDslag)
        
        format = logical(load(tmasks{s}));
        
        BOLD = zeros(length(format),length(colors));
        inds_withintmasked = inds_withintmaskedconcat(sessions==s);
        for IDnum = 1:length(colors)
            this_subnetwork_ID = colors(IDnum);
            if this_subnetwork_ID==1.5 && MSCnum==10
                this_subnetwork_ID=10.5;
            elseif this_subnetwork_ID==11.5 && (MSCnum==6 || MSCnum==7)
                this_subnetwork_ID=7.5;
            end
            
            BOLD(format,IDnum) = mean(data.data(thissub_subnetworks==this_subnetwork_ID,inds_withintmasked(inds_withintmasked>0)),1)';
        end
        
        good = true(length(colors),1); % read in spatial mask if desired
        
        % ignore pre-steady-state frames
        format(1:2) = false; % ignore first X frames
        
        FORMAT = create_blocks(format,min_block_durn,tr);
        
        %% Do the lagged correlation/covariance computation of TD matrices
        Cov = zeros([sum(good) sum(good) numel(lags)]);
        nblocks = numel(FORMAT);
        nframes = 0;
        
        % De-mean time series
        run_mean = nanmean(BOLD(format,:),1);
        BOLD = bsxfun(@minus,BOLD,run_mean);
        
        % Loop over blocks of contiguous frames
        for j = 1:numel(FORMAT)
            nframes = nframes + numel(FORMAT{j});
            FHCR = false(1,numel(format));
            FHCR(FORMAT{j}) = true;
            Cov = Cov + lagged_cov(BOLD(FHCR,good),BOLD(FHCR,good),max(lags));
        end
        
        % Normalize pairwise cross-covariance functions based on entire run
        for k = 1:numel(lags)
            Cov(:,:,k) = Cov(:,:,k)/(nframes - abs(lags(k))*nblocks);
        end
        
        % Parabolic interpolation to get peak lag/correlation
        [pl,pc] = parabolic_interp(Cov,tr);
        pl(abs(pl) > lag_lim) = nan; % Exclude long lags (generally occur when CCF is flat)
        
        % Get zero-lag correlation
        temp = Cov(:,:,lags==0);  % zero-lag correlation
        d = zeros(size(temp));
        d(logical(eye(length(temp)))) = sqrt(diag(temp));
        temp = d^(-1)*temp/d;
        temp = atanh(temp); % Fisher z transform
        temp(isnan(pl)) = nan;
        
        % Add to group running sum
        subj_lags(good,good) = pl;
        subj_ZL(good,good) = temp;
        subj_peak(good,good) = pc;
        
        grp_lags = cat(3,grp_lags,subj_lags);
        grp_lags = nansum(grp_lags,3);
        grp_ZL = cat(3,grp_ZL,subj_ZL);
        grp_ZL = nansum(grp_ZL,3);
        grp_peak = cat(3,grp_peak,subj_peak);
        grp_peak = nansum(grp_peak,3);
        
        % running sum of nans
        grp_lags_nans = grp_lags_nans + isnan(subj_lags);
        grp_ZL_nans = grp_ZL_nans + isnan(subj_ZL);
        grp_peak_nans = grp_peak_nans + isnan(subj_peak);
        
        toc
        
    end
    
    % Compute group averages
    grp_lags_mean = grp_lags ./ (numel(subjects) - grp_lags_nans);
    grp_peak_mean = grp_peak ./ (numel(subjects) - grp_peak_nans);
    grp_ZL_mean = grp_ZL ./ (numel(subjects) - grp_ZL_nans);
    grp_ZL_mean = tanh(grp_ZL_mean); % un-fisher z transform
    
    %% Sort group matrices
    
    assns = [1:num_nodes]; % import ROI network assignments for sorting TD matrix
    
    % Sort by matrices by lag
    [M,sorted_inds1] = sort(nanmean(grp_lags_mean));
    assns_sort = assns(sorted_inds1);
    
    grp_lags_temp = grp_lags_mean(sorted_inds1,sorted_inds1);
    grp_peak_temp = grp_peak_mean(sorted_inds1,sorted_inds1);
    grp_ZL_temp = grp_ZL_mean(sorted_inds1,sorted_inds1);
    
    % Sort by network
    [N,sorted_inds2] = sort(assns_sort);
    sorted_inds2 = sorted_inds2(find(N,1):end);
    
    grp_lags_mat = grp_lags_temp(sorted_inds2,sorted_inds2);
    grp_peak_mat = grp_peak_temp(sorted_inds2,sorted_inds2);
    grp_ZL_mat = grp_ZL_temp(sorted_inds2,sorted_inds2);
    
    %figure;imagesc(grp_lags_mat,[-.75,.75]);colorbar
    %figure;imagesc(grp_peak_mat,[-.7,.7]);colorbar
    %figure;imagesc(grp_ZL_mat,[-.7,.7]);colorbar
    
    %% Make group-level lag projection maps
    % Unweighted lag projection
    grp_lags_proj_unweighted = nanmean(grp_lags_mean);
    
    % Weighted lag projection (inversely weight lags by correlation magnitude
    % to reduce sampling error)
    lag_weights = tan((pi/2)*(1-abs(grp_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
    lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
    lag_weights(isnan(grp_lags_mean)) = nan;
    grp_lags_mean_wghtd = grp_lags_mean.*lag_weights;
    grp_lags_proj = nansum(grp_lags_mean_wghtd)./nansum(lag_weights);
    
    
    %[~,ordering,~] = intersect(colors,subnetwork_IDs);
    ordering = 1:length(colors);
    all_lagmat(ordering,ordering,MSCnum) = grp_lags_mat;
end
save('all_lagmat_structtogether.mat','all_lagmat','colors')

%% mean lag fig
load('all_lagmat_structtogether.mat')

meanorder = nanmean(nanmean(all_lagmat,1),3);
stdorder = nanstd(nanmean(all_lagmat,1),[],3);

neworder = [1 6 7 8 9];%[1 4 6 5 7];%[1 8 2 3 9 4 6 5 7];%

make_network_bargraph(colors(neworder),meanorder(neworder),stdorder(neworder)./sqrt(10),false)

%% mean lag fig core peripheral only
load('all_lagmat_structtogether.mat')

subnetworks_touse = [1 2 3 4 5 6 7 8 9];%[1 6 7 8 9];%

meanorder = nanmean(nanmean(all_lagmat(subnetworks_touse,subnetworks_touse,:),1),3);
stdorder = nanstd(nanmean(all_lagmat(subnetworks_touse,subnetworks_touse,:),1),[],3);



make_network_bargraph(colors(subnetworks_touse),meanorder,stdorder./sqrt(10),false)

%sigmat = zeros(length(subnetworks_touse));
%sigmat(1,2:4) = 1; sigmat(2:4,1) = 1;

avg_TDmat = nanmean(all_lagmat(subnetworks_touse,subnetworks_touse,:),3);
showmat = tril(ones(size(avg_TDmat)),-1);
parcel_correlmat_figmaker_v2_nolines_unsorted(avg_TDmat .* showmat,colors(subnetworks_touse),[-.4 .4],'test')

pmat = zeros(length(subnetworks_touse));
for i = 1:size(avg_TDmat,1)
    for j = (i+1):size(avg_TDmat,1)
        [H,P,CI,STATS] = ttest(all_lagmat(subnetworks_touse(i),subnetworks_touse(j),:));
        pmat(j,i) = P;
    end
end
sigmat_unc = single((pmat < .05) & (pmat>0));
parcel_correlmat_figmaker_v2_nolines_unsorted(avg_TDmat .* sigmat_unc,colors(subnetworks_touse),[-.4 .4],'test')
pval_corr = FDR(pmat(pmat>0),.05);
sigmat_corr = single((pmat < pval_corr) & (pmat>0));
parcel_correlmat_figmaker_v2_nolines_unsorted(avg_TDmat .* sigmat_corr,colors(subnetworks_touse),[-.4 .4],'test')

%% core v peripheral lag fig

load('all_lagmat_structtogether.mat')

corenum = 1;
peripheralnums = [6 7 8 9];

core_to_peripheral_cortex = squeeze(all_lagmat(corenum,peripheralnums,:));

make_network_bargraph(colors(peripheralnums),[nanmean(core_to_peripheral_cortex,2)'],[nanstd(core_to_peripheral_cortex,[],2)']./sqrt(10),false)