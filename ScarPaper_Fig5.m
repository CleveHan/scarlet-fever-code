%% ScarPaper_Fig5
clear;clc;close all;

%%% Tuberculosis:No; Influenza;
%%% Once: Rubella; Pertussis; Measles
%%% Twice: Mumps; ScarfletFever
Type='ScarfletFever';
load(['D:\Lab\HCL Research\Infectious Disease Open Data-Province level\呼吸道传播疾病\SortData 2004-2018\',Type])
DataInfectNum=DataInfo.DataInfectNum(2:end,2:end,:);
CityAll=DataInfo.CityAll(2:32);

load('D:\Lab\HCL Research\Infectious Disease Open Data-Province level\呼吸道传播疾病\SortData 2004-2018\NaturalFactor.mat')
Temp=Temp(:,2:end,:);
RainFall=RainFall(:,2:end,:);
Humid=Humid(:,2:end,:);
SunLight=SunLight(:,2:end,:);

% DataInfectNum=Temp;
%% Time Sample
close all
for Province=28
    DD1=reshape(squeeze(DataInfectNum(Province,:,:))',1,14:12);
    TT1=reshape(squeeze(Temp(Province,:,:))',1,14:12);
    RR1=reshape(squeeze(RainFall(Province,:,:))',1,14:12);
    HH1=reshape(squeeze(Humid(Province,:,:))',1,14:12);
    SS1=reshape(squeeze(SunLight(Province,:,:))',1,14:12);
    
%     x=1:14*12;y=DD1;tt=1:14*12;pos=find(isnan(y)==1);
%     if ~isempty(pos);x(pos)=[];y(pos)=[];end
%     DD1_interp = interp1(x, y, tt);
    DD1_interp=reshape(DD1',12,14);
    
    x=1:14*12;y=TT1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    TT1_interp = interp1(x, y, tt);TT1_interp=reshape(TT1_interp',12,14);
    
    x=1:14*12;y=RR1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    RR1_interp = interp1(x, y, tt);RR1_interp=reshape(RR1_interp',12,14);
    
    x=1:14*12;y=HH1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    HH1_interp = interp1(x, y, tt);   HH1_interp=reshape(HH1_interp',12,14); 
    
    x=1:14*12;y=SS1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    SS1_interp = interp1(x, y, tt);    SS1_interp=reshape(SS1_interp',12,14);
    
    for yy=1:14
        dd=DD1_interp(:,yy);Num_DD(yy)=length(find(dd==0));dd1=[dd-mean(dd)]/std(dd);DD1_Norm(:,yy)=dd1;
        dd=TT1_interp(:,yy);Num_TT(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);TT1_Norm(:,yy)=dd1;
        dd=RR1_interp(:,yy);Num_RR(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);RR1_Norm(:,yy)=dd1;
        dd=HH1_interp(:,yy);Num_HH(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);HH1_Norm(:,yy)=dd1;
        dd=SS1_interp(:,yy);Num_SS(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);SS1_Norm(:,yy)=dd1;
    end
    Year=1:14;
    
    subplot(4,5,1:5)
    plot(DD1_interp(:),'color',0.25*[1,1,1],'linewidth',2);xlim([0,170]);box off;set(gca,'fontsize',12,'linewidth',2,'xtick',[]);
    ylim([0,350])
    hold on;
    for ii=1:13
    plot([12*ii,12*ii],[0,350],'--','color',0*[1,1,1],'linewidth',2);
    end
    
    subplot(4,5,6:10)
    plot(DD1_Norm(:),'color',0*[1,0.5,0],'linewidth',2);xlim([0,170]);box off;set(gca,'fontsize',12,'linewidth',2,'xtick',[])
    
end
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[0 0 20 21])



cc=0;cc1=0;
for Province=1:31
    if Province==14;continue;end
    if Province==21;continue;end
    if Province==30;continue;end
    cc=cc+1;
    DD1=reshape(squeeze(DataInfectNum(Province,:,:))',1,14:12);
    TT1=reshape(squeeze(Temp(Province,:,:))',1,14:12);
    RR1=reshape(squeeze(RainFall(Province,:,:))',1,14:12);
    HH1=reshape(squeeze(Humid(Province,:,:))',1,14:12);
    SS1=reshape(squeeze(SunLight(Province,:,:))',1,14:12);
    
%     x=1:14*12;y=DD1;tt=1:14*12;pos=find(isnan(y)==1);
%     if ~isempty(pos);x(pos)=[];y(pos)=[];end
%     DD1_interp = interp1(x, y, tt);
    DD1_interp=reshape(DD1',12,14);
    
    x=1:14*12;y=TT1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    TT1_interp = interp1(x, y, tt);TT1_interp=reshape(TT1_interp',12,14);
    
    x=1:14*12;y=RR1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    RR1_interp = interp1(x, y, tt);RR1_interp=reshape(RR1_interp',12,14);
    
    x=1:14*12;y=HH1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    HH1_interp = interp1(x, y, tt);   HH1_interp=reshape(HH1_interp',12,14); 
    
    x=1:14*12;y=SS1;tt=1:14*12;pos=find(isnan(y)==1);
    if ~isempty(pos);x(pos)=[];y(pos)=[];end
    SS1_interp = interp1(x, y, tt);    SS1_interp=reshape(SS1_interp',12,14);
    
    for yy=1:14
        dd=DD1_interp(:,yy);Num_DD(yy)=length(find(dd==0));dd1=[dd-mean(dd)]/std(dd);DD1_Norm(:,yy)=dd1;
        dd=TT1_interp(:,yy);Num_TT(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);TT1_Norm(:,yy)=dd1;
        dd=RR1_interp(:,yy);Num_RR(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);RR1_Norm(:,yy)=dd1;
        dd=HH1_interp(:,yy);Num_HH(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);HH1_Norm(:,yy)=dd1;
        dd=SS1_interp(:,yy);Num_SS(yy)=length(find(isnan(dd)==1));%dd1=[dd-mean(dd)]/std(dd);SS1_Norm(:,yy)=dd1;
    end
    Year=1:14;
    
    SelYear_DD=Year(find(Num_DD==0));
    SelYear_TT=Year(find(Num_TT==0));
    SelYear_RR=Year(find(Num_RR==0));
    SelYear_HH=Year(find(Num_HH==0));
    SelYear_SS=Year(find(Num_SS==0));
    
    disp([num2str(Province),' ',num2str(length(SelYear_DD))])
    if length(SelYear_DD)==14
        cc1=cc1+1;
        DD1_Norm_All(cc1,:)=DD1_Norm(:)';
    end
    
    if Province==21;SelYear_DD=11:14;end
    if Province==26;SelYear_DD=1:14;end
    
    SelYear=intersect(SelYear_DD,SelYear_TT);
    dd1=TT1_interp(:,SelYear);dd2=DD1_Norm(:,SelYear);
%     subplot(221);plot(dd1,dd2,'ko');ylim([-3,4]);
    [r,p]=corr(dd1(:),dd2(:));r=roundn(r,-2);p=roundn(p,-4);%title(['r=',num2str(r),' p=',num2str(p)]);
    r_TT(Province)=r;p_TT(Province)=p;
    
    DataLength{Province,1}=SelYear_DD;   %% 14, 21, 30 can not do this analysis
    
    pos1=find(SelYear_DD<=7);
    TuningSmall(cc,:)=mean(dd2(:,pos1),2);
    pos1=find(SelYear_DD>=8);
    TuningLarge(cc,:)=mean(dd2(:,pos1),2);
    
    
    WarmMonth=2:7;
    ColdMonth=[8:12,1];
    DD_Warm=TuningSmall(cc,WarmMonth);
    DD_Cold=TuningSmall(cc,ColdMonth);
    pos=find(DD_Warm==max(DD_Warm));
    Pos_WarmPre_Small(cc,1)=pos(1);
    pos=find(DD_Cold==max(DD_Cold));
    Pos_ColdPre_Small(cc,1)=pos(1);
    
    DD_Warm=TuningLarge(cc,WarmMonth);
    DD_Cold=TuningLarge(cc,ColdMonth);
    pos=find(DD_Warm==max(DD_Warm));
    Pos_WarmPre_Large(cc,1)=pos(1);
    pos=find(DD_Cold==max(DD_Cold));
    Pos_ColdPre_Large(cc,1)=pos(1);
    
end
subplot(4,5,11:15)
imagesc(DD1_Norm_All)
set(gca,'linewidth',2,'fontsize',12)

% plot(mean(TuningSmall),'o-','color',[0,0.45,0.75],'linewidth',2);hold on;
% plot(mean(TuningLarge),'o-','color',[0.3,0.75,0.9],'linewidth',2);
% box off;
subplot(4,5,16:18)
for mm=1:12
    [h,pp(mm)]=ttest(TuningSmall(:,mm),TuningLarge(:,mm));
    plot(mm*4-2+randn(1,28)*0.1,TuningSmall(:,mm),'o','color',[0.3,0.75,0.9],'MarkerSize',10);hold on;
    plot(mm*4-1+randn(1,28)*0.1,TuningLarge(:,mm),'o','color',[1,0.45,0.2],'MarkerSize',10);hold on;
end
set(gca,'fontsize',12,'linewidth',2,'xtick',[])
legend('2004-2011','2011-2018')
xlim([0,48]);box off;set(gca,'linewidth',2,'fontsize',12)
set(gcf,'color','w')

subplot(4,5,19)
plot(1+randn(1,28)*0.1,(Pos_WarmPre_Small),'o','color',[0.3,0.75,0.9],'MarkerSize',10);hold on;
plot(2+randn(1,28)*0.1,(Pos_WarmPre_Large),'o','color',[1,0.6,0.4],'MarkerSize',10);ylim([2,6]);box off;
set(gca,'linewidth',2,'fontsize',12);xlim([0.5,2.5])
ylabel('Preferred month / warm season')
subplot(4,5,20)
plot(1+randn(1,28)*0.1,(Pos_ColdPre_Small),'o','color',[0,0.45,0.75],'MarkerSize',10);hold on;
plot(2+randn(1,28)*0.1,(Pos_ColdPre_Large),'o','color',[1,0.3,0.2],'MarkerSize',10);ylim([2,6]);box off;
box off;set(gca,'linewidth',2,'fontsize',12);xlim([0.5,2.5])
ylabel('Preferred month / cold season')
set(gcf,'color','w')

% figure;
% subplot(121)
% set(gca,'color',[0.3,0.75,0.9]);alpha(0.5)
% subplot(122)
% set(gca,'color',[1,0.45,0.2]);




