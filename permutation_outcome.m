function permutation_outcome(map1,nrot)
% permutation_outcome('mean_tmap_distribution_HGG+tlrc_BNA246_sqrt.txt',1000)
% permutation_outcome('mean_tmap_distribution_LGG+tlrc_BNA246_sqrt.txt',1000)

% ROI based spin test
    centroids_l = readtable('lh_BN_Atlas_210_centroids.txt','Delimiter',' ');
    centroids_r = readtable('rh_BN_Atlas_210_centroids.txt','Delimiter',' ');
    centroids_l = table2array(centroids_l);
    centroids_l = centroids_l(1:2:210,:);
    centroids_r = table2array(centroids_r);
    centroids_r = centroids_r(2:2:210,:);
    % shuffle cortical regions while keeping bilateral symmetry and spatial continuity
    cortex_perm_id = rotate_parcellation(centroids_l,centroids_r,nrot);

    % shuffle subcortical regions while keeping bilateral symmetry
    subcortical_perm_id = ones(36,nrot);
    for i=1:nrot
        subcortical_permID = randperm(36/2);
        lh_subcortical = [211:2:246];
        rh_subcortical = [212:2:246];
        lh_sub_permed = lh_subcortical(subcortical_permID);
        rh_sub_permed = rh_subcortical(subcortical_permID);
        subcortical_perm_id([1:2:36],i) = lh_sub_permed;
        subcortical_perm_id([2:2:36],i) = rh_sub_permed;
    end
    y=load(map1);
    if(size(y,1)~=246)
        y=y';
    end
    perm_id = [cortex_perm_id;subcortical_perm_id];
    y_perm = zeros(246,nrot);
    for r = 1:nrot
        for i = 1:246
            y_perm(i,r) = y(perm_id(i,r)); % permute y
        end
    end
    writematrix(y_perm,[map1 '_frequeny_perm.txt'],'Delimiter',',')
end