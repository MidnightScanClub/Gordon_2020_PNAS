MSCnums = 1:10;

iterations = 100;

columns = 9;

Homogeneity_sub = ones(max(MSCnums),columns) .* NaN;
Homogeneity_rot_sub = ones(max(MSCnums),iterations,columns) .* NaN;
Homogeneity_other_sub = ones(max(MSCnums),max(MSCnums),columns) .* NaN;


baddata = ft_read_cifti_mod('MSCmean_badverts_cort.dtseries.nii');
baddatainds = find(baddata.data); baddatainds(baddatainds>59412) = [];

rotations = ft_read_cifti_mod('Rotated_inds.dtseries.nii');


for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    disp(MSCname)
    
    corr = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegressedDconns/' MSCname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dconn.nii']);%['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN_MSM_V2/cifti_correlation_normalwall_native_MSM/' MSCname '.allses_LR_timeseries_concat_corr_s2.55.dconn.nii']);
    
    corr.data = corr.data(:,1:59412);
    
    covcorr = single(cov(corr.data));
    clear corr
    
    
    
    
    %% Subnetworks
    
    
    
    
    for col = 1:9;%columns
        fprintf('%s',['Subject ' MSCname ': Testing subnetwork Homogeneity for column ' num2str(col)])
        
        subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' MSCname '_infomap_subcortreg_ignoreverts/' MSCname '_rawassn_minsize10_regularized.dtseries.nii']);
        subnetworks.data(59413:end,:) = 0;
        subnetworks.data = subnetworks.data(:,col);
        
        
        
        subnetwork_IDs = unique(subnetworks.data); subnetwork_IDs(subnetwork_IDs<1) = [];
        
        subnetwork_IDs_temp = subnetwork_IDs;
        for IDnum = 1:length(subnetwork_IDs)
            if nnz(subnetworks.data==subnetwork_IDs(IDnum))<5
                subnetwork_IDs_temp(subnetwork_IDs_temp==subnetwork_IDs(IDnum)) = [];
            end
        end
        subnetwork_IDs = subnetwork_IDs_temp;
        
        all_real_homogeneity = zeros(length(subnetwork_IDs),1);
        pct_verts = zeros(length(subnetwork_IDs),1);
        for IDnum = 1:length(subnetwork_IDs)
            pct_verts(IDnum) = nnz(subnetworks.data==subnetwork_IDs(IDnum)) ./ nnz(subnetworks.data);
            all_real_homogeneity(IDnum) = PCA_reduction_cov_onecomp(covcorr(subnetworks.data==subnetwork_IDs(IDnum),subnetworks.data==subnetwork_IDs(IDnum)));
            
        end
        Homogeneity_sub(MSCnum,col) = nansum(all_real_homogeneity .* pct_verts);
        
        
        subnetworks_rot = rotations.data;
        subnetworks_rot(:) = 0;
        
        subnetworks_rot(rotations.data>0) = subnetworks.data(rotations.data(rotations.data>0));
        subnetworks_rot(baddatainds,:) = 0;
        
        
        
        disp(' ')
        
        
        fprintf('%s',['Subject ' MSCname ': Testing rotated subnetwork Homogeneity for column ' num2str(col)])
        
        
        all_rot_homogeneity = ones(length(subnetwork_IDs),iterations) .* NaN;
        pct_verts = zeros(length(subnetwork_IDs),iterations);
        prevstring = [];
        for iter = 1:iterations
            string = [': iteration ' num2str(iter) ' of ' num2str(iterations)];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            
            these_subnetworks = subnetworks_rot(:,iter);
            
            for IDnum = 1:length(subnetwork_IDs)
                pct_verts(IDnum,iter) = nnz(these_subnetworks==subnetwork_IDs(IDnum)) ./ nnz(these_subnetworks);
                all_rot_homogeneity(IDnum,iter) = PCA_reduction_cov_onecomp(covcorr(these_subnetworks==subnetwork_IDs(IDnum),these_subnetworks==subnetwork_IDs(IDnum)));
                
            end
            
        end
        disp(' ');
        
        sizescaled_rot_homogeneity = all_rot_homogeneity .* pct_verts;
        Homogeneity_rot_sub(MSCnum,:,col) = nansum(sizescaled_rot_homogeneity,1);
        
        
        disp(' ');
        
          fprintf('%s',['Subject ' MSCname ': Testing other subject subnetwork and rotated Homogeneity for column ' num2str(col)])
        
        prevstring = [];
        for MSCnum2 = 1:10
            if MSCnum~=MSCnum2
                
                
                subnetworks_other = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/reliability_correction/MSC' sprintf('%02i', MSCnum2) '_infomap_wacky2_subcortreg_ignoreverts/MSC' sprintf('%02i', MSCnum2) '_rawassn_minsize10_regularized.dtseries.nii']);
                subnetworks_other.data(59413:end,:) = 0;
                subnetworks_other.data = subnetworks_other.data(:,col);
                
                subnetwork_IDs = unique(subnetworks_other.data); subnetwork_IDs(subnetwork_IDs<1) = [];
                
                subnetwork_IDs_temp = subnetwork_IDs;
                for IDnum = 1:length(subnetwork_IDs)
                    if nnz(subnetworks_other.data==subnetwork_IDs(IDnum))<5
                        subnetwork_IDs_temp(subnetwork_IDs_temp==subnetwork_IDs(IDnum)) = [];
                    end
                end
                subnetwork_IDs = subnetwork_IDs_temp;
                
                all_othersub_homogeneity = ones(length(subnetwork_IDs),1) .* NaN;
                pct_verts = zeros(length(subnetwork_IDs),1);
                
                for IDnum = 1:length(subnetwork_IDs)
                    pct_verts(IDnum) = nnz(subnetworks_other.data==subnetwork_IDs(IDnum)) ./ nnz(subnetworks_other.data);
                    all_othersub_homogeneity(IDnum) = PCA_reduction_cov_onecomp(covcorr(subnetworks_other.data==subnetwork_IDs(IDnum),subnetworks_other.data==subnetwork_IDs(IDnum)));
                    
                end
                Homogeneity_other_sub(MSCnum,MSCnum2,col) = nansum(all_othersub_homogeneity .* pct_verts);
               
            end
            
        end
        disp(' ')
        
        save('Homogeneity_in_crosscol_subnetworks_wacky_subcortreg_ignoreverts.mat','Homogeneity_sub','Homogeneity_rot_sub','Homogeneity_other_sub');
    end
    
    
    

end

%% Individuals

 reverseorder = [9:-1:1];
 thresholdstrings = {'5.0%','4.0%','3.0%','2.0%','1.0%','.5%','.2%','.1%','.05%', '.02%','.01%'};
for s = 1:10;

figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])
set(gca,'FontSize',20)
hold on

plotSpread(squeeze(Homogeneity_rot_sub(s,:,reverseorder)./100),'distributionColors',[0 0 0],'MarkerSize',5)
plotSpread(squeeze(Homogeneity_other_sub(s,:,reverseorder)./100),'distributionColors',[0 1 0],'MarkerSize',20)
plotSpread(Homogeneity_sub(s,reverseorder)./100,'distributionColors',[1 0 0],'MarkerSize',40)


set(gca,'XTickLabel',thresholdstrings)
ylim([.25 .8])
set(gca,'YTickLabel',[.3 : .1 : .8])




end


%% Avg
thresholdstrings = {'5.0%','2.0%','1.0%','.5%','.2%','.1%','.05%', '.02%','.01%'};
 reverseorder = [9:-1:1];
 subinds = [1:10];
 delta_v_random = (Homogeneity_sub(subinds,reverseorder) - squeeze(nanmedian(Homogeneity_rot_sub(subinds,:,reverseorder),2))) ./ 100;
 
 figure
 set(gcf,'Position',[813 30 1102 805])
 set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
 set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])
set(gca,'FontSize',20)
hold on
h = errorbar(1:9,nanmean(delta_v_random,1),nanstd(delta_v_random,[],1)./sqrt(10),'r.');
set(h,'MarkerSize',40)
set(gca,'XTick',[1:9])
set(gca,'XTickLabel',thresholdstrings)

