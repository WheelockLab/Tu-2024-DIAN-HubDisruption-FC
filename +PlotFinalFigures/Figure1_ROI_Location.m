load('IM_13nets_246_newcolor_MNI.mat'); % load parcellation based on Seitzman 2020 300 ROI

load('MNI_coord_meshes_32k.mat');
Anat.CtxL=MNIl;Anat.CtxR=MNIr;
clear MNIl MNIr
Anat.alpha = 0.8;
    
IM.name = '';

Params.radius = 4;
figure('Position',[100,100,500,300]);
View_ROI_Modules(IM,Anat,IM.ROIxyz,Params);

% set(gcf, 'renderer', 'painters')
% print(gcf,'./Figures/Seitzman246ROIlocation','-dtiff','-r300');
%% Separately plot ROI legend
N = length(IM.Nets);
% figure('Units','inches','position',[10 10 2,3]);%[10 10 5 2]
figure('Units','inches','position',[10 10 2,2]);%[10 10 5 2]
h = gscatter(ones(1,N),ones(1,N),IM.Nets,IM.cMap,'o',50);
for i = 1:N
    set(h(i),'Color','None','MarkerFaceColor',IM.cMap(i,:));
end
legend(IM.Nets,'interpreter','none','FontSize',10,'FontWeight','Bold','location','best','Orientation','horizontal','NumColumns',2);
xlim([10,11]);
axis('off')
% print(['./Figures/Seitzman_networks_Legend'],'-dpdf');