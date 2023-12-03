function [zmatGroupMean,zmatGroupStd,grouplabel,bins] = group_fc_data(corrmat,subjectdata,groupby,nshuffle,nboot,nsubsample,indices,EYOedges)
% groupby options {'EYO','MCxSex','EYOxSex'}
% adapted from util_pipeline_DIAN where you can see the group counts and
% average age etc.
% updated 220302: the groups if seprated by EYO will be
% [-Inf,-15,-8,0,Inf], and if EYOxSex will be [-Inf,-15,0,Inf]. 
% the random seed is set to default for reproducibility for resampling methods
% addpath('/data/wheelock/data1/NLA/NLApublic/NLA_toolbox_070319/NLA_toolbox_v1.0/visualization');
%%
if ~exist('nshuffle','var')||isempty('nshuffle')
    nshuffle = 0;
end
if ~exist('nboot','var')||isempty('nboot')
    nboot = 0;
end
if ~exist('nsubsample','var')||isempty('nsubsample')
    nsubsample = 0;
end
if ~exist('indices','var')||isempty('indices')||~iscell(indices)
    indices = {};
end
%% Select groups

if ~exist('EYOedges','var')
    EYOedges = [-Inf,-15,-8,0,Inf];
end

switch groupby
    case 'EYO'
        % if separate by EYO only
        [~,edges,bins] = histcounts(subjectdata.MUT_PAR_EYO,EYOedges);
        bins(subjectdata.mutation==0 & bins>0)=bins(subjectdata.mutation==0 & bins>0)+max(bins);
        grouplabel = [sprintf('MC(EYO<%i)',EYOedges(2)),arrayfun(@(kk)sprintf('MC(%i<=EYO<%i)',EYOedges(kk),EYOedges(kk+1)),2:length(EYOedges)-2,'UniformOutput',false),sprintf('MC(EYO>=%i)',EYOedges(end-1))];
        grouplabel = [grouplabel,strrep(grouplabel,'MC','NC')];
    case 'Sex'
        % if separate by MCxSex
        [~,XEDGES,YEDGES,BINX,BINY] = histcounts2(subjectdata.mutation,subjectdata.SEX,[0,0.5,1],[1,1.5,2]);
        bins = zeros(size(BINX));
        bins(BINX~=0) = sub2ind([2,2],BINX(BINX~=0),BINY(BINX~=0));
        grouplabel = {'NC(Male)','NC(Female)','MC(Male)','MC(Female)'};
    case 'EYOxSex'
        % if separate by sex and EYO
        [~,XEDGES,YEDGES,BINX,BINY] = histcounts2(subjectdata.MUT_PAR_EYO,subjectdata.SEX,EYOedges,[1,1.5,2]);
        BINX(subjectdata.mutation==0 & BINX>0)= max(BINX)+1; % create bin4 manually which is NC and don't separate them by EYO
        bins = zeros(size(BINX));
        bins(BINX~=0) = sub2ind([max(BINX),2],BINX(BINX~=0),BINY(BINX~=0));
        grouplabel1 = [sprintf('MC(EYO<%i)_M',EYOedges(2)),arrayfun(@(kk)sprintf('MC(%i<=EYO<%i)_M',EYOedges(kk),EYOedges(kk+1)),2:length(EYOedges)-2,'UniformOutput',false),sprintf('MC(EYO>=%i)_M',EYOedges(end-1))];
        grouplabel2 = [sprintf('MC(EYO<%i)_M',EYOedges(2)),arrayfun(@(kk)sprintf('MC(%i<=EYO<%i)_M',EYOedges(kk),EYOedges(kk+1)),2:length(EYOedges)-2,'UniformOutput',false),sprintf('MC(EYO>=%i)_M',EYOedges(end-1))];
        grouplabel = [grouplabel1,'NC_M',grouplabel2,'NC_F'];
%         grouplabel = {'MC(EYO<-15)_M','MC(-15<=EYO<0)_M','MC(EYO>=0)_M','NC_M','MC(EYO<-15)_F','MC(-15<=EYO<0)_F','MC(EYO>=0)_F',...
%             'NC_F'};
    case 'CDR'
        [~,edges,bins] = histcounts(subjectdata.cdrglob,[0,0.4,0.9,Inf]);
        bins(subjectdata.mutation==0) = 4;
        grouplabel = {'MC(CDR=0)','MC(CDR=0.5)','MC(CDR>=1)','NC'};  
    case 'CDRxSex'
        [~,XEDGES,YEDGES,BINX,BINY] = histcounts2(subjectdata.cdrglob,subjectdata.SEX,[0,0.4,0.9,Inf],[1,1.5,2]);
        BINX(subjectdata.mutation==0 & BINX>0)= max(BINX)+1; % create bin4 manually which is NC and don't separate them by EYO
        bins  = zeros(size(BINX));
        bins(BINX~=0) = sub2ind([max(BINX),2],BINX(BINX~=0),BINY(BINX~=0));
        grouplabel = {'MC(CDR=0)_M','MC(CDR=0.5)_M','MC(CDR>=1)_M','NC_M'};  
        grouplabel = [grouplabel,strrep(grouplabel,'_M','_F')];
    otherwise
        error('Check your group name!');
end

if isempty(corrmat)
    zmatGroupMean = {};
    zmatGroupStd = {};
    return
end

%% Create group-average matrix
if (nshuffle>0) + (nboot>0) + (nsubsample>0) > 1
    error('only shuffle or bootstrap or , not both')
end

if nshuffle
    zmatGroupMean = {};
    zmatGroupStd = {};
    rng('default')
    for i = 1:nshuffle
        currbins = bins(randperm(length(bins)));
        zmatGroupMean = [zmatGroupMean,arrayfun(@(i)mean(corrmat(:,:,currbins==i),3),1:max(currbins),'UniformOutput',false)];
        zmatGroupStd = [zmatGroupStd,arrayfun(@(i)std(corrmat(:,:,currbins==i),[],3),1:max(currbins),'UniformOutput',false)];    
    end
elseif nboot
    zmatGroupMean = cell(1,max(bins)*nboot);
    zmatGroupStd = cell(1,max(bins)*nboot);
    rng('default')
    for i = 1:nboot
        zmatCurrMean = cell(1,max(bins));
        zmatCurrStd = cell(1,max(bins));
        for j = 1:max(bins)
            idx = randsample(find(bins==j),sum(bins==j),true);
            zmatCurrMean{j} = mean(corrmat(:,:,idx),3);
            zmatCurrStd{j} = std(corrmat(:,:,idx),[],3);
        end
        zmatGroupMean((max(bins)*(i-1)+1):max(bins)*i) = zmatCurrMean;%zmatGroupMean = [zmatGroupMean,zmatCurrMean];
        zmatGroupStd((max(bins)*(i-1)+1):max(bins)*i) = zmatCurrStd;%zmatGroupStd = [zmatGroupStd,zmatCurrStd];
    end
elseif nsubsample
    zmatGroupMean = {};
    zmatGroupStd = {};
    rng('default')
    mincounts = min(arrayfun(@(i)sum(bins==i),1:max(bins)));
    for i = 1:nsubsample
        for j = 1:max(bins)
            idx = randsample(find(bins==j),mincounts,false); % mincounts was 16
            zmatCurrMean{j} = mean(corrmat(:,:,idx),3);
            zmatCurrStd{j} = std(corrmat(:,:,idx),[],3);
        end
        zmatGroupMean = [zmatGroupMean,zmatCurrMean];
        zmatGroupStd = [zmatGroupStd,zmatCurrStd];
    end

elseif ~isempty(indices)
    zmatGroupMean = cellfun(@(C)mean(corrmat(:,:,C),3),indices,'UniformOutput',false);
    zmatGroupStd = cellfun(@(C)std(corrmat(:,:,C),[],3),indices,'UniformOutput',false);
else
    zmatGroupMean = arrayfun(@(i)nanmean(corrmat(:,:,bins==i),3),1:max(bins),'UniformOutput',false);
    zmatGroupStd = arrayfun(@(i)nanstd(corrmat(:,:,bins==i),[],3),1:max(bins),'UniformOutput',false);
end

end