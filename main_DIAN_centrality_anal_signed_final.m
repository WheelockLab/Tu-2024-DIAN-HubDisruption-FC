%% main_DIAN_centrality_anal_signed
% reference: Rubinov & Sporns, 2011 NeuroImage
clear;close all;
%% Load data and correct for data amount
allROI = readtable('DIAN_Seitzman_246.xlsx');
exclusion_id = [72,145];% the two subjects that are the only ones for their sites

subjectdata = readtable('MRI_subjectlist_N207_221213.xlsx');% subject info
subjectdata(exclusion_id,:)=[];

savedir = './postcovbat_individual_signed_mst_0.05_Z_';
Centrality = load(fullfile(savedir,'Centrality.mat'));
Centrality = Util.excludesubjects(Centrality,exclusion_id); 
fs = fields(Centrality);
for i = 1:length(fs)
    eval([fs{i},'=','Centrality.',fs{i},';'])
end

bindata = load('Mutation_CDR_bins_NCmatched.mat');
bindata = Util.excludesubjects(bindata,exclusion_id);
fs = fields(bindata);
for i = 1:length(fs)
    eval([fs{i},'=','bindata.',fs{i},';'])
end

load(fullfile(savedir,'RealNets.mat'));
rmatGroupMean(:,:,exclusion_id)=[];

Nroi = size(S,2);
NPE=Nroi*(Nroi-1)/2;                    % number of possible edges
UDidx=find(triu(ones(Nroi),1)==1);      % indices of unique conns
Nsubj = size(S,1);

try
    normcoef = (Sneg./(Spos+Sneg));
    normcoef(isnan(normcoef)) = 0;
catch
    normcoef = 0;
end


% Correcting for num of minutes - N.B.the zscores of S and Scorr are the same
X = repmat(subjectdata.totalMinutes,1,Nroi);
model = fitlm(X(:),S(:));
Scorr = model.Residuals.raw+model.Coefficients.Estimate(1);
Scorr = reshape(Scorr,size(S));

model = fitlm(X(:),Pc(:));
Pccorr = model.Residuals.raw+model.Coefficients.Estimate(1);
Pccorr = reshape(Pccorr,size(Pc));
Pc = Pccorr;

FC = reshape(rmatGroupMean,[],Nsubj)';
X = repmat(subjectdata.totalMinutes,1,length(UDidx));
model = fitlm(X(:),Util.unroll(FC(:,UDidx)));
FCUD_corr = model.Residuals.raw+model.Coefficients.Estimate(1);
FCUD_corr = reshape(FCUD_corr,size(X));

mFC = mean(FC(:,UDidx));
DMNidx = find((IM.key(:,2)==5)*(IM.key(:,2)==5)');
DMNidx = intersect(UDidx,DMNidx);
mFC_DMN = mean(FC(:,DMNidx));

%% Intermodule and intramodule strength

M = double(IM.key(:,2)==IM.key(:,2)');
[withinS,betweenS,withinS_raw,betweenS_raw,Z] = deal(NaN(size(S)));

rmat0 = rmatGroupMean;
for isubj = 1:size(rmat0,3)
    tmp = rmat0(:,:,isubj);
    tmppos = tmp.*(M==1&tmp>0);
    tmpneg = -tmp.*(M==1&tmp<0);
    correctnanpos = sum(tmppos);correctnanpos(isnan(correctnanpos)) = 0;
    correctnanneg = (sum(tmpneg));correctnanneg(isnan(correctnanneg)) = 0;
    withinS_raw(isubj,:) = correctnanpos-normcoef(isubj,:).*correctnanneg;
    withinS(isubj,:) = withinS_raw(isubj,:)./sum(M==1&tmp>0);

    correctnanpos = BCT.module_degree_zscore(tmp.*(tmp>0),IM.key(:,2))';correctnanpos(isnan(correctnanpos)) = 0; % I think technically they got the same result because it is using within K only
    correctnanneg =BCT.module_degree_zscore(-tmp.*(tmp<0),IM.key(:,2))';correctnanneg(isnan(correctnanpos)) = 0;
    Z(isubj,:) =  correctnanpos-normcoef(isubj,:).*correctnanneg;
    
    tmppos = tmp.*(M==0&tmp>0);
    tmpneg = -tmp.*(M==0&tmp<0);
    correctnanpos = sum(tmppos);correctnanpos(isnan(correctnanpos)) = 0;
    correctnanneg = (sum(tmpneg));correctnanneg(isnan(correctnanneg)) = 0;
    betweenS_raw(isubj,:) = correctnanpos-normcoef(isubj,:).*correctnanneg;
    betweenS(isubj,:) = betweenS_raw(isubj,:)./sum(M==0&tmp>0);
end
withinS(isnan(withinS))=0;
betweenS(isnan(betweenS))=0;

%% plot figures
return
