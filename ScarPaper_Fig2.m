%% AnalyzeData_Infectious_Huxi
%% Fig 2
clear;clc;close all;
%%% Tuberculosis:No; Influenza;
%%% Once: Rubella; Pertussis; Measles
%%% Twice: Mumps; ScarfletFever
Type='ScarfletFever';
load(['D:\Lab\HCL Research\Infectious Disease Open Data-Province level\呼吸道传播疾病\SortData 2004-2018\',Type])
DataInfectNum=DataInfo.DataInfectNum(:,2:end,:);
CityAll=DataInfo.CityAll;

load('D:\Lab\HCL Research\Infectious Disease Open Data-Province level\呼吸道传播疾病\SortData 2004-2018\NaturalFactor.mat')
Temp=Temp(:,2:end,:);
RainFall=RainFall(:,2:end,:);
Humid=Humid(:,2:end,:);
SunLight=SunLight(:,2:end,:);

DataInfectNum=Temp;

figure;cc=0;
for ii=[1:2,4:12] size(DataInfectNum,1)
    cc=cc+1;
%     subplot(6,2,cc)
    A=squeeze(DataInfectNum(ii,11:14,:));
    AAA=reshape(A',1,size(A,1)*12);
    AA=reshape(A',1,size(A,1)*12);AA=AA/max(AA);
    tt=1:48;
    x=tt;y=AA;pos=find(isnan(AA)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    AA_interp = interp1(x, y, tt);
    plot(AA_interp+cc,'k','linewidth',2);box off;hold on;
end
suptitle('Temp');ylim([0.5,12]);xlim([1,48])
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 15 20])
set(gca,'linewidth',2,'fontsize',13)

for ii=1:size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(Temp(ii,11:14,:));
    AAA(ii,:)=reshape(A',1,size(A,1)*12);
end

%%
CityEnglish={'China';'Beijing';'Tianjin';'Hebei';'Shanxi';'Neimenggu';'Liaoning';'Jilin';'Heilongjiang';'Shanghai';'Jiangsu';'Zhejiang';'Anhui';'Fujian';'Jiangxi';'Shandong';'Henan';'Hubei';'Hunan';'Guangdong';'Guangxi';'Hainan';'Chongqing';'Sichuan';'Guizhou';'Yunnan';'Xizang';'Shaanxi';'Gansu';'Qinghai';'Ningxia';'Xinjiang'};
DataSet(:,1)=CityEnglish;
for ii=2:32
    B=nanmean(AAA(ii-1,6:8),2);
    DataSet{ii,2}=B;
end
ColorMap='jet';FrameSel=0;ColorbarSel=0;
Range=[-20,30];
figure;
% subplot(231)
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)


for ii=2:32
    B=nanmean(AAA(ii-1,23:26),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=nanmean(AAA(ii-1,42:44),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)



%%

DataInfectNum=RainFall;

figure;cc=0;
for ii=[1:2,4:12] size(DataInfectNum,1)
    cc=cc+1;
%     subplot(6,2,cc)
    A=squeeze(DataInfectNum(ii,11:14,:));
    AA=reshape(A',1,size(A,1)*12);AA=AA/max(AA);
    tt=1:48;
    x=tt;y=AA;pos=find(isnan(AA)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    AA_interp = interp1(x, y, tt);
    plot(AA_interp+cc,'k','linewidth',2);box off;hold on;
end
suptitle('RainFall');ylim([0.5,12]);xlim([1,48])
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 15 20])
set(gca,'linewidth',2,'fontsize',13)

for ii=1:size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(DataInfectNum(ii,11:14,:));
    AAA(ii,:)=reshape(A',1,size(A,1)*12);
    
    tt=1:48;
    x=tt;y=AAA(ii,:);pos=find(isnan(AA)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    AA_interp(ii,:) = interp1(x, y, tt);
end


CityEnglish={'China';'Beijing';'Tianjin';'Hebei';'Shanxi';'Neimenggu';'Liaoning';'Jilin';'Heilongjiang';'Shanghai';'Jiangsu';'Zhejiang';'Anhui';'Fujian';'Jiangxi';'Shandong';'Henan';'Hubei';'Hunan';'Guangdong';'Guangxi';'Hainan';'Chongqing';'Sichuan';'Guizhou';'Yunnan';'Xizang';'Shaanxi';'Gansu';'Qinghai';'Ningxia';'Xinjiang'};
DataSet(:,1)=CityEnglish;
for ii=2:32
    B=nanmean(AA_interp(ii-1,6:8),2);
    DataSet{ii,2}=B;
end
ColorMap='jet';FrameSel=0;ColorbarSel=0;
Range=[0,300];
figure;
% subplot(231)
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)


for ii=2:32
    B=nanmean(AA_interp(ii-1,21:27),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=nanmean(AA_interp(ii-1,41:45),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)




%%
close all
DataInfectNum=Humid;

figure;cc=0;
for ii=[1:2,4:12] size(DataInfectNum,1)
    cc=cc+1;
%     subplot(6,2,cc)
    A=squeeze(DataInfectNum(ii,11:14,:));
    AA=reshape(A',1,size(A,1)*12);AA=AA/max(AA);
    tt=1:48;
    x=tt;y=AA;pos=find(isnan(AA)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    AA_interp = interp1(x, y, tt);
    plot(AA_interp+cc,'k','linewidth',2);box off;hold on;
end
suptitle('Humid');ylim([0.5,12]);xlim([1,48])
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 15 20])
set(gca,'linewidth',2,'fontsize',13)

for ii=1:size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(DataInfectNum(ii,11:14,:));
    AAA(ii,:)=reshape(A',1,size(A,1)*12);
end


for ii=1:size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(DataInfectNum(ii,11:14,:));
    AAA(ii,:)=reshape(A',1,size(A,1)*12);
    
    tt=1:48;
    x=tt;y=AAA(ii,:);pos=find(isnan(AA)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    AA_interp(ii,:) = interp1(x, y, tt);
end

CityEnglish={'China';'Beijing';'Tianjin';'Hebei';'Shanxi';'Neimenggu';'Liaoning';'Jilin';'Heilongjiang';'Shanghai';'Jiangsu';'Zhejiang';'Anhui';'Fujian';'Jiangxi';'Shandong';'Henan';'Hubei';'Hunan';'Guangdong';'Guangxi';'Hainan';'Chongqing';'Sichuan';'Guizhou';'Yunnan';'Xizang';'Shaanxi';'Gansu';'Qinghai';'Ningxia';'Xinjiang'};
DataSet(:,1)=CityEnglish;
for ii=2:32
    B=nanmean(AA_interp(ii-1,6:8),2);
    DataSet{ii,2}=B;
end
ColorMap='jet';FrameSel=0;ColorbarSel=0;
Range=[20,100];
figure;
% subplot(231)
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)


for ii=2:32
    B=nanmean(AA_interp(ii-1,26:28),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=nanmean(AA_interp(ii-1,44:47),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)




%%

DataInfectNum=SunLight;

figure;cc=0;
for ii=[1:2,4:12] size(DataInfectNum,1)
    cc=cc+1;
%     subplot(6,2,cc)
    A=squeeze(DataInfectNum(ii,11:14,:));
    AA=reshape(A',1,size(A,1)*12);AA=AA/max(AA);
    tt=1:48;
    x=tt;y=AA;pos=find(isnan(AA)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    AA_interp = interp1(x, y, tt);
    plot(AA_interp+cc,'k','linewidth',2);box off;hold on;
end
suptitle('SunLight');ylim([0.5,12]);xlim([1,48])
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 15 20])
set(gca,'linewidth',2,'fontsize',13)

for ii=1:size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(DataInfectNum(ii,11:14,:));
    AAA(ii,:)=reshape(A',1,size(A,1)*12);
end
for ii=1:size(DataInfectNum,1)
    cc=cc+1;
    A=squeeze(DataInfectNum(ii,11:14,:));
    AAA(ii,:)=reshape(A',1,size(A,1)*12);
    
    tt=1:48;
    x=tt;y=AAA(ii,:);pos=find(isnan(AA)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    AA_interp(ii,:) = interp1(x, y, tt);
end
CityEnglish={'China';'Beijing';'Tianjin';'Hebei';'Shanxi';'Neimenggu';'Liaoning';'Jilin';'Heilongjiang';'Shanghai';'Jiangsu';'Zhejiang';'Anhui';'Fujian';'Jiangxi';'Shandong';'Henan';'Hubei';'Hunan';'Guangdong';'Guangxi';'Hainan';'Chongqing';'Sichuan';'Guizhou';'Yunnan';'Xizang';'Shaanxi';'Gansu';'Qinghai';'Ningxia';'Xinjiang'};
DataSet(:,1)=CityEnglish;
for ii=2:32
    B=nanmean(AA_interp(ii-1,6:8),2);
    DataSet{ii,2}=B;
end
ColorMap='jet';FrameSel=0;ColorbarSel=0;
Range=[0,350];
figure;
% subplot(231)
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)


for ii=2:32
    B=nanmean(AA_interp(ii-1,23:26),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)

for ii=2:32
    B=nanmean(AA_interp(ii-1,44:47),2);
    DataSet{ii,2}=B;
end
figure;
plotChinaMap_fun(DataSet,ColorMap,FrameSel,ColorbarSel,Range)


%%
%%
clear;clc;close all;
%%% Tuberculosis:No; Influenza;
%%% Once: Rubella; Pertussis; Measles；
%%% Twice: Mumps; ScarfletFever；
Type='ScarfletFever';
load(['D:\Lab\HCL Research\Infectious Disease Open Data-Province level\呼吸道传播疾病\SortData 2004-2018\',Type])
DataInfectNum=DataInfo.DataInfectNum(:,2:end,:);
CityAll=DataInfo.CityAll;
CityEnglish={'China';'Beijing';'Tianjin';'Hebei';'Shanxi';'Neimenggu';'Liaoning';'Jilin';'Heilongjiang';'Shanghai';'Jiangsu';'Zhejiang';'Anhui';'Fujian';'Jiangxi';'Shandong';'Henan';'Hubei';'Hunan';'Guangdong';'Guangxi';'Hainan';'Chongqing';'Sichuan';'Guizhou';'Yunnan';'Xizang';'Shaanxi';'Gansu';'Qinghai';'Ningxia';'Xinjiang'};

for ii=1:size(DataInfectNum,1)
    A=squeeze(DataInfectNum(ii,:,:));
    TNum(ii,:)=nansum(A');
end
load('D:\Lab\HCL Research\Infectious Disease Open Data-Province level\呼吸道传播疾病\SortData 2004-2018\NaturalFactor.mat')
Temp=Temp(:,2:end,:);
RainFall=RainFall(:,2:end,:);
Humid=Humid(:,2:end,:);
SunLight=SunLight(:,2:end,:);

DataInfectNum=SunLight;
ccc=0;
for ii=[12:-1:7] 1:size(DataInfectNum,1)
    ccc=ccc+1;
    subplot(6,2,2*ccc-1)
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
    
    
    set(gca,'linewidth',2,'fontsize',12,'xtick',[])
end
% suptitle(Type)
set(gcf,'color','w')
% set(gcf,'Units','centimeters','Position',[0 1 10 19])
cd('D:\Lab\HCL Research\Infectious Disease Open Data-Province level\呼吸道传播疾病\Basic Figures')
% export_fig([Type,'_1'])

%%% Power Spectrum

ccc=0;
for ii=[12:-1:7] 
    ccc=ccc+1;
    subplot(6,2,2*ccc)
    A=squeeze(DataInfectNum(ii,:,:));
    cc=0;clear AAA;
    for jj=1:size(A,1)
        cc=cc+1;
%         plot(1:12,A(jj,:),'.','color',0.5*[1,1,1]);hold on;
        AAA(cc,:)=A(jj,:);
    end
    D=AAA';
    DD=D(:);
    
    tt=1:168;
    x=tt;y=DD;pos=find(isnan(DD)==1);
    if ~isempty(pos)
        x(pos)=[];y(pos)=[];
    end
    DD_interp = interp1(x, y, tt);
    
    
    [MeanS,S{ii},f{ii}]=Multi_Taper_Fourier_Transform_HCL(DD_interp-mean(DD_interp),12,6);
    plot(f{ii},S{ii},'k','linewidth',2);box off;xlim([0.3,4])
%     plot(nanmean(AAA,1),'k','linewidth',2);box off;
    title(CityEnglish{ii})
    set(gca,'linewidth',2,'fontsize',12,'xtick',[])
end
% suptitle('DataInfectNum')
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 1 10 20])


