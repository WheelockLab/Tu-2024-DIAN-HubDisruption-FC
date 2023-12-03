%% Plot mean Fc
figure('units','inches','position',[10 10 6 4]);
count = 0;
for i = 1:length(unique(bins))
    count = count+1;
    subplot(2,3,count);
    rmat = mean(rmatGroupMean(IM.order,IM.order,bins==i),3);
    Matrix_Org3(rmat,...
    IM.key,10,[-0.3,0.3],IM.cMap,0);
    title(grouplabel{i});
    set(gca,'FontSize',11);
end
%print(gcf,[savedir,'/Fc_CDRgroups_mean'],'-dpdf');
%% Plot mean and SD
figure('units','inches','position',[10 10 6 4]);
count = 0;
for i = 1:length(unique(bins))
    count = count+1;
    subplot(2,3,count);
    rmat_mean = mean(rmatGroupMean(IM.order,IM.order,bins==i),3);
    rmat_sd = std(rmatGroupMean(IM.order,IM.order,bins==i),[],3);
    Matrix_Org3(tril(rmat_mean)+triu(rmat_sd),...
    IM.key,10,[-0.3,0.3],IM.cMap,0);
    title(grouplabel{i});
    set(gca,'FontSize',11);
end
% print(gcf,[savedir,'/Fc_CDRgroups_mean_SD'],'-dpdf');

%% Example Fc from one subject
figure;
rmat = rmatGroupMean(IM.order,IM.order,1);
Matrix_Org3(rmat,...
    IM.key,10,[-0.3,0.3],IM.cMap,0);
% print(gcf,[savedir,'/Fc_example'],'-dpdf');

%% Plot Spring-Embedded Graph
G = graph(rmat);
Lwidths = 1*G.Edges.Weight/max(G.Edges.Weight);
figure;
h = plot(G,'ko-','layout','force','NodeCData',IM.key(:,2),'NodeColor','flat','MarkerSize',5,'LineWidth',Lwidths);
colormap(IM.cMap)
% print(gcf,[savedir,'/Graph_example'],'-dtiff','-r300')