function residualized_timeseries = subcort_regression(cifti,dmat,distreg)

% Output = residual_timeseries: subcortical time series for each voxel
%          after regressing out cortical time series within "distreg" mm

% Input = cifti: path to cifti time series file 
%       = dmat: path to distance matrix 
%       = distreg: regress out cortical signal from this many mm

ciftistruct = ft_read_cifti_mod(cifti);
dmat = load(dmat);

% Regress cortical timecourses within Xmm of subcort structures
residualized_timeseries = ciftistruct.data;
brainstructind = ciftistruct.brainstructure(ciftistruct.brainstructure>0);
structs = [3:18];
subcortvox = find(ismember(brainstructind,structs));
neighbors=zeros(size(brainstructind));
for i=1:length(subcortvox)    
    neighbor_cort=find(dmat(subcortvox(i),:) <= distreg & brainstructind'<=2);
    if length(neighbor_cort)
        reg = mean(ciftistruct.data(neighbor_cort,:),1);
        [~,~,r]=regress(ciftistruct.data(subcortvox(i),:)',[reg' ones(size(reg))']);
        residualized_timeseries(subcortvox(i),:)=r;  
    end
end

end

