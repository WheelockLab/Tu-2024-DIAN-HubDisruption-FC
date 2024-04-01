% Use corrected S
S = Scorr;

% run the first block of main_DIAN_centrality_anal_signed first
BrainViewCmap = jet(101);%jet(100);% hot(100);
% BrainViewCmap = interp1(linspace(0,100,10),Plot.redbluecmap(10),linspace(0,100,100));

T = readtable('DIAN_Seitzman_246.xlsx');

load('IM_13nets_246_newcolor_MNI.mat'); % load parcellation based on Seitzman 2020 300 ROI

load('MNI_coord_meshes_32k.mat');
Anat.CtxL=MNIl;Anat.CtxR=MNIr;
clear MNIl MNIr

disp(savedir)

%%
for gp = [4]%[1:4,6:8]   
    %% if looking at raw mean S
    S_means = NaN(max(bins),Nroi);
    for ii = 1:max(bins)
        S_means(ii,:) = mean(S(bins==ii,:));
    end
    S_baseline = S_means(gp,:);
    maptounity = @(v)(v-min(S_means(:)))/range(S_means(:));
    S_baseline = floor(maptounity(S_baseline)*100)+1;
    ROI.coord = IM.ROIxyz(IM.key(:,2)>0,:);
    ROI.color = BrainViewCmap(S_baseline(IM.key(:,2)>0),:);
    ROI.radius = repmat(4,sum(IM.key(:,2)>0),1);
    % Plot 
    Anat.ctx = 'std'
    Anat.alpha = 0.8; % was 1
    figure('Units','Inches','position',[1,1,4,3]);
    subplot(5,1,1:2);
    Anat.view = 'lat';
    Plot.Draw_ROIs_on_Cortex(Anat,ROI)
    subplot(5,1,3:4);
    Anat.view = 'med';
    Plot.Draw_ROIs_on_Cortex(Anat,ROI)
    subplot(5,1,5);
    colormap(BrainViewCmap)

    hbar = colorbar('location','south');caxis([min(S_means(:)),max(S_means(:))]);hbar.Label.String = 'Strength';hbar.FontWeight = 'bold';hbar.FontSize = 12;
    hbar.Label.Position(2) = -hbar.Label.Position(2);
    axis('off');
    print(gcf,fullfile(savedir,['Sraw_gp',num2str(gp)]),'-dtiff','-r300');
    close all;
%%
for     cutoff = 85:10:95 %75:10:95 %2023.09.21 found a bug in my code that this cutoff should not use the normalized data because I got rid of the intercept when scaling to 0-100
    %% Plot with threshold
    S_means = NaN(max(bins),Nroi);
    for ii = 1:max(bins)
        S_means(ii,:) = mean(S(bins==ii,:));
    end
    S_baseline = S_means(gp,:);
    idx = S_baseline>prctile(S_baseline,cutoff);
    T.gyrus(idx)
    ROI.coord = IM.ROIxyz(idx,:);
    ROI.color = repmat([0 0 0],sum(idx),1); %BrainViewCmap(S_baseline(idx),:);
    ROI.radius = repmat(4,sum(idx),1);
    Anat.ctx = 'std'
    Anat.alpha = 0.8;% was 1
    figure('Units','Inches','position',[1,1,4,3]);
    subplot(5,1,1:2);
    Anat.view = 'lat';
    Plot.Draw_ROIs_on_Cortex(Anat,ROI)
    subplot(5,1,3:4);
    Anat.view = 'med';
    Plot.Draw_ROIs_on_Cortex(Anat,ROI)
    subplot(5,1,5);
    text(0.3,0.5,[num2str(cutoff),' th percentile'],'Units','normalized','FontSize',12,'FontWeight','Bold');
    axis off
    print(gcf,fullfile(savedir,['Sraw_gp',num2str(gp),'_hubs_cutoff_',num2str(cutoff)]),'-dtiff','-r300');

end
end
close all


return
%% Plot colorbar separately
figure('position',[100 100 400 400]);
colormap(BrainViewCmap);axis off;
hbar = colorbar('south');
caxis([min(S_means(:)),max(S_means(:))]);hbar.Label.String = 'Strength';hbar.FontWeight = 'bold';hbar.FontSize = 12;
set(gca,'FontSize',15,'FontWeight','Bold');
filestr = ['./Figures/StrengthColorbar']
print(gcf,filestr,'-dpdf');
