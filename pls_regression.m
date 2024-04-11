% the source of code:
% https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md
% and Mapping gene transcription and neurocognition  across human neocortex

% This script runs PLS on gene expression and neurosynth probability maps.
% Significance of latent variables is assessed using a
% spatial autocorrelation-preserving permutation test. Correlation between
% gene and term scores are cross-validated using a distance-based set
% assignment (see fcn_crossval_pls_brain_obvs.m). Most contributing terms
% are retained and scores are distributed among structural and functional
% networks.

% PLS code (pls_analysis.m) can be downloaded at
% http://pls.rotman-baycrest.on.ca/source/ ("Latest PLS Applications")

function pls_regression(group)
%% PLS analysis

% set up PLS analysis
X = csvread('/prot/bhb_Tumor/transcription_Image_Coupling/database/allAreasExpression_sorted.csv',1,1);
X = X(1:2:end,:)
Y_perm = readmatrix(['mean_tmap_' group '+tlrc_BNA246_sqrt.txt_frequeny_perm.txt']);
Y_perm = Y_perm(1:2:end,:);
Y = readmatrix(['mean_tmap_' group '+tlrc_BNA246_sqrt.txt']);
Y = Y';
Y = Y(1:2:end,:);

X = zscore(X);
Y = zscore(sqrt(Y));
Y_perm = zscore(sqrt(Y_perm));

nnodes = size(X,1);
nterms = size(Y,2);
ngenes = size(X,2);

% behav pls
option.method = 3;
option.num_boot = 1000;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
option.stacked_behavdata = Y;

exp{1} = X;

result = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
% save('result.mat','result')

%% spin test
% this code comes from pls_analysis.m and is modified to account for a
% spatial autocorrelation-preserving permutation test

nspins = 1000;                 % number of permutations ("spins")
s_spins = zeros(nterms,nspins); % singular values
option.method = 3;              % set up PLS
option.num_boot = 0;
option.num_perm = 0;
exp{1} = X;
for k = 1:nspins    
    option.stacked_behavdata = Y_perm(:,k);  % permute neurosynth matrix
    
    datamatsvd=rri_xcor(option.stacked_behavdata,exp{1},0); % refer to pls_analysis.m
    [r,c] = size(datamatsvd);
    if r <= c
        [pu, sperm, pv] = svd(datamatsvd',0);
    else
        [pv, sperm, pu] = svd(datamatsvd,0);
    end
    
    %  rotate pv to align with the original v
    rotatemat = rri_bootprocrust(result.v,pv);
 
    %  rescale the vectors
    pv = pv * sperm * rotatemat;

    sperm = sqrt(sum(pv.^2));
    
    s_spins(:,k) = sperm;
end

sprob = zeros(nterms,1); % p-value for each latent variable

for k = 1:nterms % get permuted (via spin test) p-values
    sprob(k) = (1+(nnz(find(s_spins(k,:)>=result.s(k)))))/(1+nspins);
end  
disp(sprob)

rvalue=result.u;
fid = fopen('/prot/bhb_Tumor/transcription_Image_Coupling/database/allAreasExpression_sorted.csv');
genename = strsplit(fgetl(fid), ',');
genename = genename(2:end);
genename = genename';
fclose(fid);
result = table(genename,rvalue);
writetable(sortrows(result,2,'descend'),[group '_gene_loadings_sqrt.txt']);
end