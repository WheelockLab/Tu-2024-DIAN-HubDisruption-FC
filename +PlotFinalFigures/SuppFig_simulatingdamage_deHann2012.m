%% Simulation
% Figure 1B/1C

cd('/data/wheelock/data1/people/Cindy/DIAN'); % where the centrality results were saved
savedir = './postcovbat_individual_signed_complete_Z_/';
load(fullfile(savedir,'Centrality.mat'));
load('Mutation_CDR_bins_NCmatched.mat')

refgroup = bins==4; % NC match 1

meanS = sort(mean(S(refgroup,:)));
zmeanS = zscore(meanS);

% three possible scenarios:
rate_random = 0.5; % uniform rate
rate_targetedhubs = 0.5*(meanS-min(meanS))./(max(meanS)-min(meanS)); % larger rate for larger S
rate_targetednonhubs = 1-0.5*(meanS-min(meanS))./(max(meanS)-min(meanS)); % smaller rate for largerS

randomdegS = meanS.*(1-rate_random); % final strength with a uniform rate
actdegS = meanS.*(1-rate_targetedhubs); % final strength with a larger rate for larger S
actinvdegS = meanS.*(1-rate_targetednonhubs); % final strength with a smaller rate for larger S


%% Original HDI (Figure 1B)
brandom = regress([randomdegS-meanS]',[ones(size(meanS))',meanS']);
btargetedhubs = regress([actdegS-meanS]',[ones(size(meanS))',meanS']);
btargetednonhubs = regress([actinvdegS-meanS]',[ones(size(meanS))',meanS']);

figure('Units','inches','position',[10 10 3 2.3]);hold on;
plot(meanS,actdegS-meanS,'LineWidth',3,'color','r','LineStyle' ,'--');
plot(meanS,randomdegS-meanS,'LineWidth',3,'color','c','LineStyle' ,'--');
plot(meanS,actinvdegS-meanS,'LineWidth',3,'color','g','LineStyle' ,'--');

hline = refline(brandom(2),brandom(1));hline.Color = 'c';hline.LineStyle = '-';
hline = refline(btargetedhubs(2),btargetedhubs(1));hline.Color = 'r';hline.LineStyle = '-';
hline = refline(btargetednonhubs(2),btargetednonhubs(1));hline.Color = 'g';hline.LineStyle = '-';
text(0.7,0.7,sprintf('\\kappa_{s}=%1.3f',brandom(2)),'units','normalized','color','c','FontWeight','Bold');
text(0.7,0.8,sprintf('\\kappa_{s}=%1.3f',btargetedhubs(2)),'units','normalized','color','r','FontWeight','Bold');
text(0.7,0.9,sprintf('\\kappa_{s}=%1.3f',btargetednonhubs(2)),'units','normalized','color','g','FontWeight','Bold');

legend('Targeted attack at hubs','Random attack','Targeted attack at non-hubs','location','SE');legend('boxoff');
ylabel({'Difference in strength'; 'from reference'});xlabel('Reference strength');
set(gca,'FontWeight','Bold','FontSize',10);
title('Original');
 print('./Figures/hubdisruptiondemo_Achard','-dpdf');
% print('./Figures/hubdisruptiondemo_Achard','-dtiff','-r300');
%% Ratio based HDI (Figure 1C)
brandom = regress([(randomdegS./meanS-1)*100]',[ones(size(zmeanS))',zmeanS']);
btargetedhubs = regress([(actdegS./meanS-1)*100]',[ones(size(zmeanS))',zmeanS']);
btargetednonhubs = regress([(actinvdegS./meanS-1)*100]',[ones(size(zmeanS))',zmeanS']);

figure('Units','inches','position',[10 10 3 2.3]);hold on;
plot(zmeanS,(actdegS./meanS-1)*100,'LineWidth',2,'color','r','LineStyle' ,'--');
plot(zmeanS,(randomdegS./meanS-1)*100,'LineWidth',2,'color','c','LineStyle' ,'--');
plot(zmeanS,(actinvdegS./meanS-1)*100,'LineWidth',2,'color','g','LineStyle' ,'--');

hline = refline(brandom(2),brandom(1));hline.Color = 'c';hline.LineStyle = '-';
hline = refline(btargetedhubs(2),btargetedhubs(1));hline.Color = 'r';hline.LineStyle = '-';
hline = refline(btargetednonhubs(2),btargetednonhubs(1));hline.Color = 'g';hline.LineStyle = '-';
text(0.7,0.7,sprintf('\\kappa_{s}=%1.1f',brandom(2)),'units','normalized','color','c','FontWeight','Bold');
text(0.7,0.8,sprintf('\\kappa_{s}=%1.1f',btargetedhubs(2)),'units','normalized','color','r','FontWeight','Bold');
text(0.7,0.9,sprintf('\\kappa_{s}=%1.1f',btargetednonhubs(2)),'units','normalized','color','g','FontWeight','Bold');

legend('Targeted attack at hubs','Random attack','Targeted attack at non-hubs','location','SE');legend('boxoff');
ylabel({'% Difference in strength'; 'from reference'});
xlabel('Reference strength (zscore)');
set(gca,'FontWeight','Bold','FontSize',10);
title('New');
% print('./Figures/hubdisruptiondemo_Tu','-dpdf');
% print('./Figures/hubdisruptiondemo_Tu','-dtiff','-r300');
