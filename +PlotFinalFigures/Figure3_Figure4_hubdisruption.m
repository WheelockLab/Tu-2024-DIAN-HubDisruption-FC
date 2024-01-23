%% Hub disruption

cmap = BrBG_cmap(8);cmap = cmap([6:8,3:-1:1],:);

[p_MC,HDI.slope,HDI.intercept,HDI.p,HDI.R2,HDI.slopeupper,HDI.slopelower,HDI.interceptupper,HDI.interceptlower,HDI.r] = deal(NaN(4,1));
[ypred,yci] = deal(cell(1,8));

refgroup = bins==4;
refgroupstr = 'YoungNC';

S = Scorr;
Pc = Pccorr;

S_ref = mean(S(refgroup,:));
Pc_ref = mean(Pc(refgroup,:));
Z_ref = mean(Z(refgroup,:));
%% Plot 
clear ax
figure('position',[100 100 650 400]);
hubness = 'S';
count = 0;
for igroup = [1:6]
    count = count+1;
    switch hubness
        case 'S'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = zscore(S_ref);
            mstr2 = 'S(z-score)';
        case 'Pc'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = zscore(Pc_ref);
            mstr2 = 'Pc(z-score)';
        case 'Z'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = Z_ref;
            mstr2 = 'Z';
    end
    
    testgroup = bins==igroup;
    differencetest = (mean(m(testgroup,:))./mref-1)*100;
    
    outliers = false(size(differencetest));
    
    mdl = fitlm(m2(~outliers),differencetest(~outliers));
    [ypred{igroup},yci{igroup}] = predict(mdl,sort(m2'),'Simultaneous',true);

    [b_test,b_testInt,~,~,stats_MC] = regress(differencetest(~outliers)',[ones(size(mref,2)-sum(outliers),1),m2(~outliers)']);
    p_MC(igroup) = stats_MC(3);
    r = corr(m2',differencetest');
      
    HDI.slope(igroup) = b_test(2);HDI.intercept(igroup) = b_test(1);
    HDI.p(igroup) = stats_MC(3);HDI.R2(igroup) = stats_MC(1);
    HDI.slopeupper(igroup) = b_testInt(2,1);HDI.interceptupper(igroup) = b_testInt(1,1);
    HDI.slopelower(igroup) = b_testInt(2,2);HDI.interceptlower(igroup) = b_testInt(1,2);
    HDI.r(igroup) = r;
    HDI.F(igroup) = stats_MC(2);
    
    ax(igroup) = subplot(2,3,count);hold on;
    if igroup==4
        axis off
        continue
    end
    scatter(m2(~outliers),differencetest(~outliers),20,cmap(igroup,:),'d','filled');alpha(0.5);
    plot(sort(m2),ypred{igroup},'-','Color',cmap(igroup,:))
    p1 = patch([sort(m2),sort(m2,'descend')],[yci{igroup}(:,1);flipud(yci{igroup}(:,2))],cmap(igroup,:),'EdgeColor','none');
    p1.FaceVertexAlphaData = 0.2;p1.FaceAlpha = 'flat' ;
    set(gca,'FontSize',10,'FontWeight','Bold');
    title(grouplabel{igroup});
    text(0.5,0.8,{sprintf('\\kappa_{%s} = %1.1f',mstr,HDI.slope(igroup));sprintf('R^2 = %1.2f',HDI.R2(igroup))},'Units','Normalized','FontSize',12)
    xlabel(sprintf('<%s>,%s',grouplabel{4},mstr2));
end
linkaxes(ax,'xy');
xl = xlim; yl = ylim;
subplot(2,3,1);
ylabel(sprintf('%% %s difference',mstr));
% mkdir(['./HDI_results/',mstr,'_',mstr2])
% save(['./HDI_results/',mstr,'_',mstr2,'/HDI_CDRgroup_threshold',sprintf('%1.2f',thresholds),'.mat'],'HDI','thresholds','savedir')

%% test for interaction
differencetest = arrayfun(@(i)(mean(m(bins==i,:))./mref-1)*100,1:6,'UniformOutput',false);

tbl = table();
tbl.agebin = repelem([1:6],Nroi)';
tbl.difference = cell2mat(differencetest)';
tbl.ref = repmat(m2,1,6)';

idx = (tbl.agebin==2|tbl.agebin==3);
mdl = fitlm([tbl.ref(idx),tbl.agebin(idx)==3],tbl.difference(idx),'interactions')
anova_stats = anova(mdl)
partialetasq = anova_stats.SumSq(end-1)/(anova_stats.SumSq(end-1)+anova_stats.SumSq(end))
%% abeta + versus abeta-
abetaposCDR0=subjectdata.mutation==1 & subjectdata.cdrglob==0 & subjectdata.PIB_fSUVR_rsf_TOT_CORTMEAN>1.42;
abetanegCDR0=subjectdata.mutation==1 & subjectdata.cdrglob==0 & subjectdata.PIB_fSUVR_rsf_TOT_CORTMEAN<=1.42;
grouplabel{7} = 'A\beta-';grouplabel{8}='A\beta+';
clear ax
figure('position',[100 100 650 400]);
hubness = 'S';
count = 0;
for igroup = 7:8 % abeta+ or abeta-
    count = count+1;
    switch hubness
        case 'S'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = zscore(S_ref);
            mstr2 = 'S(z-score)';
        case 'Pc'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = zscore(Pc_ref);
            mstr2 = 'Pc(z-score)';
        case 'Z'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = Z_ref;
            mstr2 = 'Z';
    end
    if igroup==7
        testgroup = abetanegCDR0;
    else
        testgroup = abetaposCDR0;
    end
    differencetest = (mean(m(testgroup,:))./mref-1)*100;
    
    outliers = false(size(differencetest));
    
    mdl = fitlm(m2(~outliers),differencetest(~outliers));
    [ypred{igroup},yci{igroup}] = predict(mdl,sort(m2'),'Simultaneous',true);

    [b_test,b_testInt,~,~,stats_MC] = regress(differencetest(~outliers)',[ones(size(mref,2)-sum(outliers),1),m2(~outliers)']);
    p_MC(igroup) = stats_MC(3);
    r = corr(m2',differencetest');
      
    HDI.slope(igroup) = b_test(2);HDI.intercept(igroup) = b_test(1);
    HDI.p(igroup) = stats_MC(3);HDI.R2(igroup) = stats_MC(1);
    HDI.slopeupper(igroup) = b_testInt(2,1);HDI.interceptupper(igroup) = b_testInt(1,1);
    HDI.slopelower(igroup) = b_testInt(2,2);HDI.interceptlower(igroup) = b_testInt(1,2);
    HDI.r(igroup) = r;
    HDI.F(igroup) = stats_MC(2);
    
    ax(igroup) = subplot(2,3,count);hold on;
    scatter(m2(~outliers),differencetest(~outliers),20,cmap(1,:),'d','filled');alpha(0.5);
    plot(sort(m2),ypred{igroup},'-','Color',cmap(1,:))
    p1 = patch([sort(m2),sort(m2,'descend')],[yci{igroup}(:,1);flipud(yci{igroup}(:,2))],cmap(1,:),'EdgeColor','none');
    p1.FaceVertexAlphaData = 0.2;p1.FaceAlpha = 'flat' ;
    set(gca,'FontSize',10,'FontWeight','Bold');
    title(grouplabel{igroup});
    text(0.5,0.8,{sprintf('\\kappa_{%s} = %1.1f',mstr,HDI.slope(igroup));sprintf('R^2 = %1.2f',HDI.R2(igroup))},'Units','Normalized','FontSize',12)
    xlabel(sprintf('<%s>,%s',grouplabel{4},mstr2));
    ylim(yl);
    xlim(xl)
    yticks([-100,0,100]);
end
linkaxes(ax,'xy');
subplot(2,3,1);
ylabel(sprintf('%% %s difference',mstr));
% print(gcf,fullfile(savedir,['hub_disruption_index(ratio)_',mstr,'_',mstr2,'_AbetaCDR0groups']),'-dpdf');

differencetest = [mean(m(abetanegCDR0,:)./mref-1)';mean(m(abetaposCDR0,:)./mref-1)'];
mdl = fitlm([repmat(m2',2,1),[zeros(Nroi,1);ones(Nroi,1)]],differencetest,'interactions');
anova_stats = anova(mdl)


%% calculate individual HDI with the reference group defined    
[p_MC,HDI.slope,HDI.intercept,HDI.p,HDI.R2,HDI.slopeupper,HDI.slopelower,HDI.interceptupper,HDI.interceptlower,HDI.r] = deal(NaN(Nsubj,1));

for isubj = [1:size(S,1)] % change this
    switch hubness
        case 'S'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = zscore(S_ref);
            mstr2 = 'S(z-score)';
        case 'Pc'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = zscore(Pc_ref);
            mstr2 = 'Pc(z-score)';
        case 'Z'
            m = S;
            mref = S_ref;
            mstr = 'S';
            m2 = Z_ref;
            mstr2 = 'Z';
    end
    
    differencetest = [m(isubj,:)./mref-1]*100;
    [b_test,b_testInt,~,~,stats_MC] = regress(differencetest',[ones(size(mref,2),1),m2']);
    p_MC(isubj) = stats_MC(3);
    
    r = corr(m2',differencetest');
    
    HDI.slope(isubj) = b_test(2);HDI.intercept(isubj) = b_test(1);
    HDI.p(isubj) = stats_MC(3);HDI.R2(isubj) = stats_MC(1);
    HDI.slopeupper(isubj) = b_testInt(2,1);HDI.interceptupper(isubj) = b_testInt(1,1);
    HDI.slopelower(isubj) = b_testInt(2,2);HDI.interceptlower(isubj) = b_testInt(1,2);
    HDI.r(isubj) = r;
end
% mkdir(['./HDI_results/individuals/',mstr,'_',mstr2])
% save(['./HDI_results/individuals/',mstr,'_',mstr2,'/HDI_CDRgroup_threshold',sprintf('%1.2f',thresholds),'.mat'],'HDI','thresholds','savedir')
return
%%  Plot errorbar + scatter (jitterplot)
v  = arrayfun(@(ii)HDI.slope(bins==ii),1:max(bins),'UniformOutput',false);
y = cellfun(@mean,v);
e = cellfun(@std,v);


figure('units','inches','position',[10 10 5 3]);hold on;
clear h
for j = 1:max(bins)
    jj = rem(j,6);
    if jj==0
        jj = 6;
    end
    h(j) = errorbar(jj,y(j),e(j),'LineWidth',3,'CapSize',10,'LineStyle','none','Marker','d','MarkerFaceColor',cmap(jj,:));
    h(j).Color = cmap(jj,:);
    x = 0.05.*randn(sum(bins==j),1)+jj;
    scatter(x,v{j},10,'filled',...
        'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4);
end
xticks(1:max(bins));
xlim([0,max(bins)+1]);
xticklabels(grouplabel);
xtickangle(45)
ylabel('\kappa_S');
set(gca,'FontSize',15,'FontWeight','Bold')

comparewhat = 1;
switch comparewhat
    case 1
        % ANOVA + star across groups
        [P,ANOVATAB,STATS] = anova1(HDI.slope,bins,'off');
        [comparison,means] = multcompare(STATS,'display','off','ctype','bonferroni'); % default is 'tukey-kramer',best for one-way ANOVA
        p = comparison(:,6);
        p = mafdr(p/length(p),'BHFDR',true)
        Util.sigstar(mat2cell(comparison(:,1:2),ones(1,size(comparison,1)),2),p)
        legend(h,grouplabel,'location','northeastoutside');
%                 print(gcf,fullfile(savedir,['Kappa',mstr2,'_jitterplot(CDRyoungref)']),'-dtiff'); % comparing pairwise
    case 2
        % stars compare to zero
        [~,p] = cellfun(@ttest,v); % passed normality test so using ttest
        p = mafdr(p,'BHFDR',true)
        hl=refline(0,0);hl.LineStyle = '--';hl.LineWidth = 1;hl.Color = 'k';
        arrayfun(@(i)Plot.makeSignificanceBar([i,i],double(y(i)+e(i)+diff(ylim)*0.05),p(i)),1:max(bins),'UniformOutput',false);
        legend(h,grouplabel,'location','northeastoutside');
        title(mstr2)
%         print(gcf,fullfile(savedir,['Kappa',mstr2,'_jitterplot2(CDRyoungref)']),'-dtiff'); %compare to 0
end
%% Plot jitter (separate by mutation)
agebin = mod(bins,3);agebin(agebin==0)=3;
mval = HDI.slope;

v  = arrayfun(@(ii)mval(bins==ii,:),1:max(bins),'UniformOutput',false);
y = cellfun(@mean,v);
e = cellfun(@std,v);
% x-axis is age bins
xplot = [1:3];
[~,ANOVATAB,STATS] = anova1(mval(subjectdata.mutation==1),agebin(subjectdata.mutation==1),'off');
% post-hoc test for age bins for MC
[comparison,means] = multcompare(STATS,'display','off','ctype','lsd'); % default is 'tukey-kramer',best for one-way ANOVA
p = mafdr(comparison(:,6),'BHFDR',true); % my dirty way of using FDR instead of bonferroni
comparison = [comparison(:,1:2),p];
etasquared = ANOVATAB{2,2}/ANOVATAB{end,2}

clear h ax
figure('units','inches','position',[10 10 2.5 4]);
ax(1) = subplot(2,1,1);hold on;
for j = 1:3
    jj = rem(j,6);
    if jj==0
        jj = 6;
    end
    h(j) = errorbar(xplot(jj),y(j),e(j),'LineWidth',3,'CapSize',10,'LineStyle','none','Marker','d','MarkerFaceColor',cmap(jj,:));
    h(j).Color = cmap(jj,:);
    x = 0.05.*randn(sum(bins==j),1)+xplot(jj);
    scatter(x,v{j},10,'filled',...
        'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4);
end

Plot.sigstar(mat2cell(comparison(:,1:2),ones(1,size(comparison,1)),2),comparison(:,3))
xticks(1:3);
xlim([0,4]);

xticklabels(cellfun(@(C)C(4:end-1),grouplabel(1:3),'UniformOutput',false));
xtickangle(45)
ylabel(['\kappa_{',mstr,'}'])
set(gca,'FontSize',12,'FontWeight','Bold')

% legend(h,grouplabel(1:3),'location','northeastoutside');
% print(gcf,fullfile(savedir,['Kappa',mstr2,'_jitterplot3(CDRyoungref)_MConly']),'-dpdf'); % comparing pairwise

xplot = [1:3];
[~,ANOVATAB,STATS] = anova1(mval(subjectdata.mutation==0),agebin(subjectdata.mutation==0),'off');
% post-hoc test for age bins for MC
[comparison,means] = multcompare(STATS,'display','off','ctype','lsd'); % default is 'tukey-kramer',best for one-way ANOVA
p = mafdr(comparison(:,6),'BHFDR',true);
comparison = [comparison(:,1:2),p];

% report mean and sd for each group
for i = 0:1
    for j = 1:3
        tmp = mval(subjectdata.mutation==i & agebin==j);
        m = mean(tmp);
        s  = std(tmp);
        fprintf('%1.3f,%1.3f\n',m,s)
    end
end
d = Util.computeCohen_d(mval(subjectdata.mutation==1 & agebin==1),mval(subjectdata.mutation==1 & agebin==3))

clear h
ax(2) = subplot(2,1,2);hold on;
% figure('units','inches','position',[10 10 2.5 2]);hold on;
for j = 4:6
    jj = rem(j,6);
    if jj==0
        jj = 6;
    end
    h(j) = errorbar(xplot(jj-3),y(j),e(j),'LineWidth',3,'CapSize',10,'LineStyle','none','Marker','d','MarkerFaceColor',cmap(jj,:));
    h(j).Color = cmap(jj,:);
    x = 0.05.*randn(sum(bins==j),1)+xplot(jj-3);
    scatter(x,v{j},10,'filled',...
        'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4);
end
Plot.sigstar(mat2cell(comparison(:,1:2),ones(1,size(comparison,1)),2),comparison(:,3))
xticks(1:3);
xlim([0,4]);
xticklabels(cellfun(@(C)C(4:end),grouplabel(4:6),'UniformOutput',false));
xtickangle(30)
ylabel(['\kappa_{',mstr,'}'])
set(gca,'FontSize',12,'FontWeight','Bold')

linkaxes(ax);
% legend(h(4:6),grouplabel(4:6),'location','northeastoutside');
% print(gcf,fullfile(savedir,['Kappa',mstr2,'_jitterplot3(CDRyoungref)_NConly']),'-dpdf'); % comparing pairwise
% % print(gcf,fullfile(savedir,['Kappa',mstr,'_',mstr2,'_jitterplot3(CDRyoungref)_separatebymutation']),'-dpdf'); % comparing pairwise

%% Plot amyloid + and amyloid - MC CDR = 0 (checked with Util.normality test that the distribution is still normal)
abetaposCDR0=subjectdata.mutation==1 & subjectdata.cdrglob==0 & subjectdata.PIB_fSUVR_rsf_TOT_CORTMEAN>1.42;
abetanegCDR0=subjectdata.mutation==1 & subjectdata.cdrglob==0 & subjectdata.PIB_fSUVR_rsf_TOT_CORTMEAN<=1.42;
mval = HDI.slope;

[~,p] = cellfun(@ttest,{mval(abetaposCDR0),mval(abetanegCDR0)}); % passed normality test so using ttest
p = mafdr(p,'BHFDR',true)%p*length(p); % Bonferroni correction

figure('units','inches','position',[10 10 2.5 4]);
subplot(2,1,1);hold on;
h = errorbar(2,mean(mval(abetaposCDR0)),std(mval(abetaposCDR0)),'LineWidth',3,'CapSize',10,'LineStyle','none','Marker','d','MarkerFaceColor',cmap(1,:));
h.Color = cmap(1,:);
x = 0.05.*randn(sum(abetaposCDR0),1)+2;
scatter(x,mval(abetaposCDR0),10,'filled',...
    'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4);
if p(1)<0.001
    nstars = '***'
elseif p(2)<0.01
    nstars = '**'
elseif p(2)<0.05
    nstars = '*'
else
    nstars = ''
end 
% plot(1,mean(mval(abetaposCDR0))+1.2*std(mval(abetaposCDR0)),nstars,'Color','k');

h = errorbar(1,mean(mval(abetanegCDR0)),std(mval(abetanegCDR0)),'LineWidth',3,'CapSize',10,'LineStyle','none','Marker','d','MarkerFaceColor',cmap(1,:));
h.Color = cmap(1,:);
if p(2)<0.001
    nstars = '***'
elseif p(2)<0.01
    nstars = '**'
elseif p(2)<0.05
    nstars = '*'
else
    nstars = ''
end   
% plot(2,mean(mval(abetanegCDR0))+1.2*std(mval(abetanegCDR0)),nstars,'Color','k');
x = 0.05.*randn(sum(abetanegCDR0),1)+1;
scatter(x,mval(abetanegCDR0),10,'filled',...
    'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4);
% Plot.hline(0);
xlim([0,3]);
xticks([1,2]);
xticklabels({'A\beta-','A\beta+'});
ylabel(['\kappa_{',mstr,'}'])
% legend('MC(CDR=0)','location','northwestoutside');
set(gca,'FontSize',12,'FontWeight','Bold');
ylim([-40,40])

[h,p,ci,stats] = ttest2(mval(abetaposCDR0),mval(abetanegCDR0)) % this test for whether the Abeta+ and - groups are the same and they are not significantly different
d = Util.computeCohen_d(mval(abetaposCDR0),mval(abetanegCDR0),'independent')
% print(gcf,fullfile(savedir,['Kappa',mstr,'_',mstr2,'_jitterplot(CDRyoungref)_MC_CDR0_abeta']),'-dpdf'); % comparing pairwise

%% Plot e4 + and e4 - MC CDR = 0 (checked with Util.normality test that the distribution is still normal)
e4status = any(subjectdata.apoe==[24,34,44],2);
apoe4posCDR0=subjectdata.mutation==1 & subjectdata.cdrglob==0 & e4status;
apoe4negCDR0=subjectdata.mutation==1 & subjectdata.cdrglob==0 & (~e4status);

mval = HDI.slope;

[~,p] = cellfun(@ttest,{mval(apoe4posCDR0),mval(apoe4negCDR0)}); % passed normality test so using ttest
p = p*length(p); % Bonferroni correction

figure('units','inches','position',[10 10 5 3]);hold on;
h = errorbar(1,mean(mval(apoe4posCDR0)),std(mval(apoe4posCDR0)),'LineWidth',3,'CapSize',10,'LineStyle','none','Marker','d','MarkerFaceColor',cmap(1,:));
h.Color = cmap(1,:);
x = 0.05.*randn(sum(apoe4posCDR0),1)+1;
scatter(x,mval(apoe4posCDR0),10,'filled',...
    'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4);
if p(1)<0.001
    nstars = '***'
elseif p(2)<0.01
    nstars = '**'
elseif p(2)<0.05
    nstars = '*'
else
    nstars = ''
end 
text(1,mean(mval(apoe4posCDR0))+1.2*std(mval(apoe4posCDR0)),nstars,'Color','k');

h = errorbar(2,mean(mval(apoe4negCDR0)),std(mval(apoe4negCDR0)),'LineWidth',3,'CapSize',10,'LineStyle','none','Marker','d','MarkerFaceColor',cmap(1,:));
h.Color = cmap(1,:);
if p(2)<0.001
    nstars = '***'
elseif p(2)<0.01
    nstars = '**'
elseif p(2)<0.05
    nstars = '*'
else
    nstars = ''
end   
text(2,mean(mval(apoe4negCDR0))+1.2*std(mval(apoe4negCDR0)),nstars,'Color','k');
x = 0.05.*randn(sum(apoe4negCDR0),1)+2;
scatter(x,mval(apoe4negCDR0),10,'filled',...
    'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4);
Plot.hline(0);
xlim([0,3]);
xticks([1,2]);
xticklabels({'\epsilon4+','\epsilon4-'});
ylabel(['\kappa_{',mstr,'}'])
legend('MC(CDR=0)','location','northwestoutside');
set(gca,'FontSize',12,'FontWeight','Bold');

[h,p] = ttest2(mval(apoe4posCDR0),mval(apoe4negCDR0)) % this test for whether the Abeta+ and - groups are the same and they are not significantly different

print(gcf,fullfile(savedir,['Kappa',mstr,'_',mstr2,'_jitterplot(CDRyoungref)_MC_CDR0_apoe4']),'-dpdf'); % comparing pairwise
