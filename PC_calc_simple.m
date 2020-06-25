function PCs = PC_calc_simple(thresh_corrmat,communities)

parcel_degree = sum(thresh_corrmat,2);
moduleconnectionssum = zeros(size(thresh_corrmat,1),1);
community_IDs = unique(communities);
community_IDs(community_IDs<1) = [];
for communitynum = 1:length(community_IDs)
    communityID = community_IDs(communitynum);
    communityindices = communities==communityID;
    parcel_community_degree = sum(thresh_corrmat(:,communityindices),2);
    parcel_community_ratio = (parcel_community_degree ./ parcel_degree);
    moduleconnectionssum = moduleconnectionssum + (parcel_community_ratio.^2);
end

PCs = 1-moduleconnectionssum;

