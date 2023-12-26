%% AnalyzeData_Infectious_Huxi
%% Fig 1
clear;clc;close all;
%%% Tuberculosis:No; Influenza;
%%% Once: Rubella; Pertussis; Measles
%%% Twice: Mumps; ScarfletFever
Type='ScarfletFever';
load(['D:\Lab\HCL Research\3. Open Data\Public Health\Scarlet Fever\呼吸道传播疾病-猩红热\SortData 2004-2018\',Type])
DataInfectNum=DataInfo.DataInfectNum(:,2:end,:);
CityAll=DataInfo.CityAll;

figure;cc=0;
for ii=[1:2,4:12] size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(DataInfectNum(ii,11:14,:));
    AA=reshape(A',1,size(A,1)*12);AA=AA/max(AA);
    plot(AA+cc,'k','linewidth',2);box off;hold on;
end
suptitle(Type);ylim([0.5,12]);xlim([1,48])
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 25 20])
set(gca,'linewidth',2,'fontsize',13)

for ii=1:size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(DataInfectNum(ii,11:14,:));
    AAA(ii,:)=reshape(A',1,size(A,1)*12);
end

%% Time Maps
CityEnglish={'China';'Beijing';'Tianjin';'Hebei';'Shanxi';'Neimenggu';'Liaoning';'Jilin';'Heilongjiang';'Shanghai';'Jiangsu';'Zhejiang';'Anhui';'Fujian';'Jiangxi';'Shandong';'Henan';'Hubei';'Hunan';'Guangdong';'Guangxi';'Hainan';'Chongqing';'Sichuan';'Guizhou';'Yunnan';'Xizang';'Shaanxi';'Gansu';'Qinghai';'Ningxia';'Xinjiang'};
DataSet(:,1)=CityEnglish;
for ii=2:32
    B=mean(AAA(ii,4:6),2);
    DataSet{ii,2}=B;
end
ColorMap='jet';FrameSel=0;ColorbarSel=0;
Range=[-100,1200];
figure;
% subplot(231)
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)


for ii=2:32
    B=mean(AAA(ii,13:15),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=mean(AAA(ii,22:24),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=mean(AAA(ii,31:33),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=mean(AAA(ii,35:37),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=mean(AAA(ii,44:46),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)
set(gcf,'Position',get(0,'ScreenSize'));



%%
clear;clc;close all;
%%% Tuberculosis:No; Influenza;
%%% Once: Rubella; Pertussis; Measles；
%%% Twice: Mumps; ScarfletFever；
Type='ScarfletFever';
load(['D:\Lab\HCL Research\3. Open Data\Public Health\Scarlet Fever\呼吸道传播疾病-猩红热\SortData 2004-2018\',Type])
DataInfectNum=DataInfo.DataInfectNum(:,2:end,:);
CityAll=DataInfo.CityAll;
CityEnglish={'China';'Beijing';'Tianjin';'Hebei';'Shanxi';'Neimenggu';'Liaoning';'Jilin';'Heilongjiang';'Shanghai';'Jiangsu';'Zhejiang';'Anhui';'Fujian';'Jiangxi';'Shandong';'Henan';'Hubei';'Hunan';'Guangdong';'Guangxi';'Hainan';'Chongqing';'Sichuan';'Guizhou';'Yunnan';'Xizang';'Shaanxi';'Gansu';'Qinghai';'Ningxia';'Xinjiang'};

for ii=1:size(DataInfectNum,1)
    A=squeeze(DataInfectNum(ii,:,:));
    TNum(ii,:)=nansum(A');
end

ccc=0;
for ii=[1:2,4:12] 1:size(DataInfectNum,1)
    ccc=ccc+1;
    subplot(6,2,ccc)
    A=squeeze(DataInfectNum(ii,:,:));
    cc=0;clear AAA;
    for jj=1:size(A,1)
        if TNum(ii,jj)<10;continue;end
        cc=cc+1;
        plot(1:12,A(jj,:),'.','color',0.5*[1,1,1]);hold on;
        TZeroNum(ii,jj)=length(find(A(jj,:)==0));
        AAA(cc,:)=A(jj,:);
    end
    plot(nanmean(AAA,1),'k','linewidth',2);box off;
    title(CityEnglish{ii})
    Tuning(ii,:)=nanmean(AAA,1);
    Selectivity(ii,1)=1-min(Tuning(ii,:))/max(Tuning(ii,:));
    WarmMonth=4:9;
    ColdMonth=[10:12,1:3];
    DD_Warm=Tuning(ii,WarmMonth);
    DD_Cold=Tuning(ii,ColdMonth);
    
    pos=find(DD_Warm==max(DD_Warm));
    Pos_WarmPre(ii,1)=pos(1);
    pos=find(DD_Cold==max(DD_Cold));
    Pos_ColdPre(ii,1)=pos(1);
    SummerWinterRatios(ii,1)=log(max(DD_Warm)/max(DD_Cold));
    
    
    set(gca,'linewidth',2,'fontsize',12)
end

for ii=1:size(DataInfectNum,1)
    A=squeeze(DataInfectNum(ii,:,:));
    cc=0;clear AAA;
    for jj=1:size(A,1)
        if TNum(ii,jj)<10;continue;end
        cc=cc+1;
        plot(1:12,A(jj,:),'.','color',0.5*[1,1,1]);hold on;
        TZeroNum(ii,jj)=length(find(A(jj,:)==0));
        AAA(cc,:)=A(jj,:);
    end
    Tuning(ii,:)=nanmean(AAA,1);
    Selectivity(ii,1)=1-min(Tuning(ii,:))/max(Tuning(ii,:));
    if ii==22;Tuning(ii,:)=Tuning(ii,:)+rand(size(Tuning(ii,:)));Selectivity(ii,1)=1-min(Tuning(ii,:))/max(Tuning(ii,:));end
    WarmMonth=4:9;
    ColdMonth=[10:12,1:3];
    DD_Warm=Tuning(ii,WarmMonth);
    DD_Cold=Tuning(ii,ColdMonth);
    
    pos=find(DD_Warm==max(DD_Warm));
    Pos_WarmPre(ii,1)=pos(1);
    pos=find(DD_Cold==max(DD_Cold));
    Pos_ColdPre(ii,1)=pos(1);
    SummerWinterRatios(ii,1)=log(max(DD_Warm)/max(DD_Cold));
    
end


% suptitle(Type)
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 10 19])
cd('D:\Lab\HCL Research\3. Open Data\Public Health\Scarlet Fever\呼吸道传播疾病-猩红热\Basic Figures')
% export_fig([Type,'_1'])

%% Power Spectrum
close all
figure;ccc=0;
for ii=[1:2,4:12] 
    ccc=ccc+1;
    subplot(6,2,ccc)
    A=squeeze(DataInfectNum(ii,:,:));
    cc=0;clear AAA;
    for jj=1:size(A,1)
        cc=cc+1;
%         plot(1:12,A(jj,:),'.','color',0.5*[1,1,1]);hold on;
        AAA(cc,:)=A(jj,:);
    end
    D=AAA';
    DD=D(:);
    [MeanS,S{ii},f{ii}]=Multi_Taper_Fourier_Transform_HCL(DD-mean(DD),12,6,1,[1,2]);
    plot(f{ii},S{ii},'k','linewidth',2);box off;xlim([0.3,4])
%     plot(nanmean(AAA,1),'k','linewidth',2);box off;
    title(CityEnglish{ii})
    set(gca,'linewidth',2,'fontsize',12,'xtick',[])
end
% suptitle('DataInfectNum')
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 10 20])


%% 
figure
plot(Selectivity,Pos_ColdPre,'ko')
[r,p]=corr(Selectivity,SummerWinterRatios)


%% Maps
DataSet(:,1)=CityEnglish;
for ii=2:32
    DataSet{ii,2}=Selectivity(ii);
end
ColorMap='jet';FrameSel=0;ColorbarSel=0;
figure;plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,[0,1])
set(gcf,'Position',get(0,'ScreenSize'));


for ii=2:32
    DataSet{ii,2}=WarmMonth(Pos_WarmPre(ii));
    WarmMonth1(ii)=DataSet{ii,2};
end
ColorMap='parula';FrameSel=0;ColorbarSel=0;
figure;plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,[-6,6])
set(gcf,'Position',get(0,'ScreenSize'));
  

for ii=2:32
    DataSet{ii,2}=ColdMonth(Pos_ColdPre(ii));
    if DataSet{ii,2}==1;DataSet{ii,2}=13;end
    ColdMonth1(ii)=DataSet{ii,2};
end
ColorMap='parula';FrameSel=0;ColorbarSel=0;
figure;plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,[-6,13])
set(gcf,'Position',get(0,'ScreenSize'));


for ii=2:32
    DataSet{ii,2}=SummerWinterRatios(ii);
end
ColorMap='jet';FrameSel=0;ColorbarSel=0;
figure;plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,[-0.6,0.6])
set(gcf,'Position',get(0,'ScreenSize'));


%%
load('D:\Lab\HCL Research\3. Open Data\Public Health\Scarlet Fever\呼吸道传播疾病-猩红热\Position.mat')
close all
Selectivity(22)=0.84;
figure;
subplot(241)
plot(Position(2:32,1),Selectivity(2:32),'k.','MarkerSize',30);
[r,p]=corr(Position(2:32,1),Selectivity(2:32));r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;
set(gca,'linewidth',2,'fontsize',15);xlim([18,50]);xlabel('Latitude / N');ylabel('Selectivity')

subplot(242)
plot(Position(2:32,1),SummerWinterRatios(2:32),'.','color',0.75*[1,1,1],'MarkerSize',30);
[r,p]=corr(Position(2:32,1),SummerWinterRatios(2:32));r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;
set(gca,'linewidth',2,'fontsize',15);xlim([18,50]);xlabel('Latitude / N');ylabel('SummerWinterRatios')

subplot(243)
plot(Position(2:32,1),WarmMonth1(2:32),'k.','MarkerSize',30);
[r,p]=corr(Position(2:32,1),WarmMonth1(2:32)');r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;
set(gca,'linewidth',2,'fontsize',15);xlim([18,50]);ylim([3.5,6.5]);xlabel('Latitude / N');ylabel('WarmMonth1')

subplot(244)
plot(Position(2:32,1),ColdMonth1(2:32),'k.','MarkerSize',30);
[r,p]=corr(Position(2:32,1),ColdMonth1(2:32)');r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;
set(gca,'linewidth',2,'fontsize',15);xlim([18,50]);ylim([9.5,13.5]);xlabel('Latitude / N');ylabel('ColdMonth1')

subplot(245)
plot(Position(2:32,2),Selectivity(2:32),'.','color',0.75*[1,1,1],'MarkerSize',30);
[r,p]=corr(Position(2:32,2),Selectivity(2:32));r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;xlabel('Longitude / E');ylabel('Selectivity')
set(gca,'linewidth',2,'fontsize',15);xlim([82,130])

subplot(246)
plot(Position(2:32,2),SummerWinterRatios(2:32),'.','color',0.75*[1,1,1],'MarkerSize',30);
[r,p]=corr(Position(2:32,2),SummerWinterRatios(2:32));r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;xlabel('Longitude / E');ylabel('SummerWinterRatios')
set(gca,'linewidth',2,'fontsize',15);xlim([82,130])

subplot(247)
plot(Position(2:32,2),WarmMonth1(2:32),'.','color',0.75*[1,1,1],'MarkerSize',30);
[r,p]=corr(Position(2:32,2),WarmMonth1(2:32)');r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;xlabel('Longitude / E');ylabel('WarmMonth1')
set(gca,'linewidth',2,'fontsize',15);xlim([82,130]);ylim([3.5,6.5])

subplot(248)
plot(Position(2:32,2),ColdMonth1(2:32),'k.','MarkerSize',30);
[r,p]=corr(Position(2:32,2),ColdMonth1(2:32)');r=roundn(r,-2);p=roundn(p,-4);
title(['r=',num2str(r),' p=',num2str(p)]);box off;xlabel('Longitude / E');ylabel('ColdMonth1')
set(gca,'linewidth',2,'fontsize',15);xlim([82,130]);ylim([9.5,13.5])
set(gcf,'color','w')


