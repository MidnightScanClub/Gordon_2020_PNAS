MSCnums = 1:10;

taskcols = [7 5 6];

iterations = 1000;

columns = 9;

Rsq_sub = ones(10,columns,length(taskcols)) .* NaN;
Rsq_other_sub = ones(10,10,columns,length(taskcols)) .* NaN;
Rsq_rot_sub = ones(10,columns,iterations,length(taskcols)) .* NaN;


sphere = gifti(['/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.coord.gii']); %any 32k_fsLR left hem sphere
sphereLcoords = sphere.vertices;
sphere = gifti(['/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.coord.gii']); %any 32k_fsLR right hem sphere
sphereRcoords = sphere.vertices;
allspherecoords = [sphereLcoords ; sphereRcoords];

rotations = ft_read_cifti_mod('/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/Rotated_inds.dtseries.nii');

for taskcolnum = 1:length(taskcols)
    taskcol = taskcols(taskcolnum);
    
    for MSCnum = MSCnums
        MSCname = ['MSC' sprintf('%02i',MSCnum)];
        disp(MSCname)
        
        taskdata = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC//Analysis_V2/Tasks/' MSCname '_allcontrasts_fs_LR_32k_LR_smooth2.55.dscalar.nii']);
        
        for col = 1:columns
            disp(['Subject ' MSCname ': Testing subnetwork and rotated subnetwork task Rsq for column ' num2str(col) ' in task number ' num2str(taskcolnum)])
            
            subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' MSCname '_infomap_subcortreg_ignoreverts/' MSCname '_rawassn_minsize10_regularized.dtseries.nii']);
            subnetworks.data(59413:end,:) = 0;
            subnetworks.data(subnetworks.data<0) = 0;
            subnetworks.data = subnetworks.data(:,col);
            
            
            subnetworks_rot = rotations.data;
            subnetworks_rot(:) = 0;
            
            subnetworks_rot(rotations.data>0) = subnetworks.data(rotations.data(rotations.data>0));
            
            
            subnetwork_inds_withvals = find(subnetworks.data>0);
            
            this_Rsq = 0;
            inds_withvals = setdiff(subnetwork_inds_withvals,find(isinf(taskdata.data(:,taskcol))));
            try [~,T,STATS,~]=anovan(taskdata.data(inds_withvals,taskcol),{subnetworks.data(inds_withvals)},'display','off');
                baseRsq = T{2,2} ./ T{end,2};
                adjustment_factor = (nnz(inds_withvals)-1) ./ (nnz(inds_withvals) - length(unique(subnetworks.data(inds_withvals))) + 1);
                this_Rsq = 1 - adjustment_factor + (adjustment_factor .* baseRsq);
            catch; end
            if this_Rsq>0
                Rsq_sub(MSCnum,col,taskcolnum) = this_Rsq;
            end
            
            
            
            
            parfor iter = 1:iterations
                
                these_subnetworks = subnetworks_rot(:,iter);
                
                inds_withvals_random = intersect(subnetwork_inds_withvals,setdiff(find(these_subnetworks),find(isinf(taskdata.data(:,taskcol)))));
                vals = these_subnetworks(inds_withvals_random);
                
                try [~,T,STATS,~]=anovan(taskdata.data(inds_withvals_random,taskcol),{vals},'display','off');
                    baseRsq = T{2,2} ./ T{end,2};
                    adjustment_factor = (nnz(inds_withvals_random)-1) ./ (nnz(inds_withvals_random) - length(unique(vals)) + 1);
                    this_Rsq = 1 - adjustment_factor + (adjustment_factor .* baseRsq);
                    if this_Rsq>0
                        Rsq_rot_sub(MSCnum,col,iter,taskcolnum) = this_Rsq;
                    end
                catch
                end
            end
            
            
            
            
            disp(['Subject ' MSCname ': Testing other subnetwork and rotated other subnetwork task Rsq for column ' num2str(col) ' in task number ' num2str(taskcolnum)])
            
            
            for MSCnum2 = 1:10
                if MSCnum~=MSCnum2
                    
                    subnetworks_other = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/MSC' sprintf('%02i', MSCnum2) '_infomap_subcortreg_ignoreverts/MSC' sprintf('%02i', MSCnum2) '_rawassn_minsize10_regularized.dtseries.nii']);
                    subnetworks_other.data(59413:end,:) = 0;
                    subnetworks_other.data = subnetworks_other.data(:,col);
                    
                    
                    inds_withvals = intersect(subnetwork_inds_withvals,setdiff(find(subnetworks_other.data > 0),find(isinf(taskdata.data(:,taskcol)))));
                    
                    try [~,T,STATS,~]=anovan(taskdata.data(inds_withvals,taskcol),{subnetworks_other.data(inds_withvals)},'display','off');
                        baseRsq = T{2,2} ./ T{end,2};
                        adjustment_factor = (nnz(inds_withvals)-1) ./ (nnz(inds_withvals) - length(unique(subnetworks_other.data(inds_withvals))) + 1);
                        this_Rsq = 1 - adjustment_factor + (adjustment_factor .* baseRsq);
                        Rsq_other_sub(MSCnum,MSCnum2,col,taskcolnum) = this_Rsq;
                    catch; end
                    
                    
                    subnetworks_rot = rotations.data;
                    subnetworks_rot(:) = 0;
                    
                    subnetworks_rot(rotations.data>0) = subnetworks_other.data(rotations.data(rotations.data>0));

                end
                
            end
            
            save('Task_activation_in_crosscol_subnetworks_subcortreg_ignoreverts.mat','Rsq_sub','Rsq_rot_sub','Rsq_other_sub')
        end
        
        
        
    end
end



%% Individuals

 taskcols = [1 2 3 ];
 %reverseorder = [6:-1:1];
 reverseorder = [9:-1:1];
 thresholdstrings = {'5.0%','2.0%','1.0%','.5%','.2%','.1%','.05%', '.02%','.01%'};
 for taskcol = taskcols;
for s = 1:10;

figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])
set(gca,'FontSize',20)
hold on

plotSpread(squeeze(Rsq_rot_sub(s,reverseorder,:,taskcol))','distributionColors',[0 0 0],'MarkerSize',5)
plotSpread(squeeze(Rsq_other_sub(s,:,reverseorder,taskcol)),'distributionColors',[0 1 0],'MarkerSize',20)
plotSpread(Rsq_sub(s,reverseorder,taskcol),'distributionColors',[0 0 1],'MarkerSize',40)
set(gca,'XTickLabel',thresholdstrings)
ylim([0 .8])
set(gca,'YTickLabel',[0 :.1 :.8])
export_fig(['MSC' sprintf('%02i',s) ' taskR2 task ' num2str(taskcol) '.png'],gca)




end
 end


 

%% Cross-task Avg
 thresholdstrings = {'5.0%','2.0%','1.0%','.5%','.2%','.1%','.05%', '.02%','.01%'};
 taskcols = [1 2 3];
 reverseorder = [9:-1:1];
 subinds = [1:10];

versionstring = version;
versionnum = str2num(versionstring(1));

disp('Rsq vs random: difference vs prev column')
delta_v_random = ((Rsq_sub(subinds,reverseorder,taskcols) -  squeeze(nanmean(Rsq_rot_sub(subinds,reverseorder,:,taskcols),3))));
for i = 1:8
    thisthresh = delta_v_random(:,i,:);
    nextthresh = delta_v_random(:,i+1,:);
    [H,P,CI,STATS] = ttest(nextthresh(:),thisthresh(:));
    disp([thresholdstrings{i+1} ' vs ' thresholdstrings{i} ': P=' num2str(P) '; T=' num2str(STATS.tstat) '; diff=' num2str(mean(nextthresh(:)-thisthresh(:)))])
end

figure
 set(gcf,'Position',[813 30 1102 805])
 set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
 set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])
set(gca,'FontSize',20)
hold on
if versionnum > 8
shadedErrorBar(1:9,nanmean(nanmean(delta_v_random,3),1),nanstd(nanmean(delta_v_random,3),[],1)./sqrt(10),'lineprops',{'r-','LineWidth',2})
else
h = errorbar(1:9,nanmean(nanmean(delta_v_random,3),1),nanstd(nanmean(delta_v_random,3),[],1)./sqrt(10),'r.');
set(h,'MarkerSize',40)
end
set(gca,'XTick',[1:9])
set(gca,'XTickLabel',thresholdstrings)


disp('Rsq vs other: difference vs prev column')
delta_v_other = ((Rsq_sub(subinds,reverseorder,taskcols) -  squeeze(nanmean(Rsq_other_sub(subinds,:,reverseorder,taskcols),2))));
for i = 1:8
    thisthresh = delta_v_other(:,i,:);
    nextthresh = delta_v_other(:,i+1,:);
    [H,P,CI,STATS] = ttest(nextthresh(:),thisthresh(:));
    disp([thresholdstrings{i+1} ' vs ' thresholdstrings{i} ': P=' num2str(P) '; T=' num2str(STATS.tstat) '; diff=' num2str(mean(nextthresh(:)-thisthresh(:)))])
end

figure
 set(gcf,'Position',[813 30 1102 805])
 set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
 set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])
set(gca,'FontSize',20)
hold on
if versionnum > 8
    shadedErrorBar(1:9,nanmean(nanmean(delta_v_other,3),1),nanstd(nanmean(delta_v_other,3),[],1)./sqrt(10),'lineprops',{'b-','LineWidth',2})
else
    h = errorbar(1:9,nanmean(nanmean(delta_v_other,3),1),nanstd(nanmean(delta_v_other,3),[],1)./sqrt(10),'b.');
    set(h,'MarkerSize',40)
end
set(gca,'XTick',[1:9])
set(gca,'XTickLabel',thresholdstrings)

