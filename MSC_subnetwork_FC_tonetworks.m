MSCnums = 1:10;

networks_totest = [1 3 7 9 10 2]; %DMN FP Language CON SM Vis


subnetwork_IDs = [1 7.5 16 11.5 10.8 1.5 3.75 15.5 10.5];



FC_strengths = zeros(max(MSCnums),length(subnetwork_IDs),length(networks_totest));

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
        
    data = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Scott/MSC_Subcortical/CorticalRegTimeSeries/' MSCname '_LR_surf_subcort_222_32k_fsLR_smooth2.55_subcortreg_20mm_regression.dtseries.nii']);
    
    thissub_subnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Subnetworks/' MSCname '_infomap_subcortreg_ignoreverts/' MSCname '_rawassn_minsize10_regularized_DMNmatch_v2_recolor.dtseries.nii']);
    thissub_networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V2/infomap/' MSCname '/' MSCname '_rawassn_minsize400_regularized_recolored_cleaned.dscalar.nii']);
    
    for n = 1:length(networks_totest)
        
        
        network_inds = thissub_networks.data==networks_totest(n);
        
        
        for s = 1:length(subnetwork_IDs)
            
            disp([MSCname ': network ' num2str(networks_totest(n)) '; subnetwork ' num2str(subnetwork_IDs(s))])
            
            this_subnetwork_ID = subnetwork_IDs(s);
            if this_subnetwork_ID==1.5 && MSCnum==10
                this_subnetwork_ID=10.5;
            elseif this_subnetwork_ID==11.5 && (MSCnum==6 || MSCnum==7)
                this_subnetwork_ID=7.5;
            end
                
            
            subnetwork_inds = abs(thissub_subnetworks.data-this_subnetwork_ID)<.01;%(thissub_subnetworks.data>=subnetwork_IDs(s)) & (thissub_subnetworks.data<(subnetwork_IDs(s)+1));
            
            these_networkinds = network_inds(1:59412) & ~(subnetwork_inds(1:59412));
            
            if any(subnetwork_inds) && any(these_networkinds)
            
            subnetwork_tc = mean(data.data(subnetwork_inds,:),1);
            
            
            
            network_tc = mean(data.data(these_networkinds,:),1);
            
            
            FC_strengths(MSCnum,s,n) = FisherTransform(paircorr_mod(subnetwork_tc',network_tc'));
            
            end
            
        end
    end
end
FC_strengths(FC_strengths==0) = NaN;
save('FC_tonetworks.mat','FC_strengths')
%%





            
            
            
            
            

    
    
