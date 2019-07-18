%%%%%% Matlab script 
% /home/jr3830/proj/auditory/codes/dogear/bfdata/ProcessDogearNMAvg.m
%%%% Latex Doc
% /home/jr3830/proj/auditory/docs/2014-streaming-model/jr-streaming-model.tex
function ProcPFData(fname)
%% 
% Code to accompany the paper Byrne, Rinzel and Rankin (2019)
% Entrainment of stream segregation in a dynamic environment
% Contact: james.rankin@gmail.com
% If you use or adapt this code acknowledge us by citing our paper
%
% Euler-Muryama scheme to timestep model

make_colors

ifexport=0;

% set where to save data/plots based on which pc this runs on
[~,pcname]=system('hostname');
if strcmp(strtrim(pcname),'jrankin-dt')
datadir='/home/jr3830/Insync/proj/auditory/codes/auditory-streaming/pfdata/';
plotdir='/home/jr3830/Insync/RankinExpt/DFMod/HearRes/makefigs/';
else; svdir=pwd;plotdir=pwd;
end
cd(datadir)
% fnames={
% 'SDIB_pfp4x1x72_PFperiod0_20_null0_0_g_sig_15Jul2019_160654.mat';
% };
synd=0;
ifinterp=0;
ifsave=0;
ifnormalise=0;
par2int=1;
meandur=[5.3,8.05];
xlims={[0 16];[-3,3]};
if ifnormalise
ylims={[0,3.5];[3.5,11]};
else
ylims={[0,30];[3.5,11]};
end
l=1;
xlabels={'$\phi$ (st)';'``Attend Group"  ||``Attend Split"  ||'};
% if ~exist('SwitchData','var')
load(fname);
% end

MinWindow=0.05;
t0=0.5;
Tmin=.5;
Tmax=60;
ifplot=0;

NumRepeats=size(SwitchData,3);
NumConditions=size(SwitchData,1);
NumCases=size(SwitchData,2);
PropIntMtx=zeros(NumRepeats,NumConditions);
DurationsIntByCond={};
DurationsSegByCond={};
SwTimesByCond={};
parvals=par1vals;
DFIdx=1:length(par1vals);
propFirstInt=zeros(NumConditions,1);
meanFirstInt=zeros(NumConditions,1);
meanFirstSeg=zeros(NumConditions,1);
period=1;

for i=1:NumConditions
    for k=1:NumCases
    NumFstIntTmp=0;
    NumFstSegTmp=0;
    meanFirstIntTmp=[];
    meanFirstSegTmp=[];
    DurationsIntTmp=[];
    DurationsSegTmp=[];
    SwTimesByCondTmp=[];
    for j=1:NumRepeats
        if  isfield(SwitchData,'uFiltCoarse')
            midval=5;
            uFilt=SwitchData(i,k,j).uFiltCoarse;
            vFilt=SwitchData(i,k,j).vFiltCoarse;
            trsc=tintCoarse/rsctime;
        else
            midval=12;
            uFilt=SwitchData(i,k,j).uFilt;
            vFilt=SwitchData(i,k,j).vFilt;
            trsc=tint/rsctime;
        end
        [NumberOfSwitches,First,Durations,DurationsInt,DurationsSeg,PropIntegrated,...
            SwitchTimes,SwDirections]=...
            ProcSwTimesT0(trsc,t0,MinWindow,uFilt,vFilt,ifplot);
            PropIntMtx(j,i)=PropIntegrated;
            %         figure;plot(trsc,uFilt,trsc,vFilt)
            DurationsInt(DurationsInt<Tmin)=[];
            DurationsInt(DurationsInt>Tmax)=[];
            DurationsSeg(DurationsSeg<Tmin)=[];
            DurationsSeg(DurationsSeg>Tmax)=[];
            
            DurationsIntTmp=[DurationsIntTmp;DurationsInt];
            DurationsSegTmp=[DurationsSegTmp;DurationsSeg];
            SwTimesByCondTmp=[SwTimesByCondTmp;SwitchTimes',SwDirections'];
%         if i==size(DurationsIntByCond,2)
%             DurationsIntByCond{i}=[DurationsIntByCond{i};DurationsInt];
%         else
%             DurationsIntByCond{i}=DurationsInt;
%         end
%         if i==size(DurationsSegByCond,2)
%             DurationsSegByCond{i}=[DurationsSegByCond{i};DurationsSeg];
%         else
%             DurationsSegByCond{i}=DurationsSeg;
%         end
        if First(2)==1
            NumFstIntTmp=NumFstIntTmp+1;
            meanFirstIntTmp=[meanFirstIntTmp;First(1)];
        else
            NumFstSegTmp=NumFstIntTmp+1;
            meanFirstSegTmp=[meanFirstSegTmp;First(1)];
        end
            
    end
    if ~isempty(DurationsIntTmp)
    DurationsIntByCond{i,k}=DurationsIntTmp;
    else        
    DurationsIntByCond{i,k}=[];
    end
    if ~isempty(DurationsSegTmp)
    DurationsSegByCond{i,k}=DurationsSegTmp;
    else
    DurationsSegByCond{i,k}=[];
    end
    if ~isempty(SwTimesByCondTmp)
    SwTimesByCond{i,k}=SwTimesByCondTmp;
    else
    SwTimesByCond{i,k}=[];
    end
%     meanFirstIntTmp;
%     meanFirstSegTmp;
%     meanFirstInt(i)=mean(meanFirstIntTmp);
%     meanFirstSeg(i)=mean(meanFirstSegTmp);
%     propFirstInt(i)=NumFstIntTmp/NumRepeats;
    end
end


meanInt=zeros(size(DurationsIntByCond));
meanSeg=zeros(size(DurationsSegByCond));
stdInt=zeros(size(DurationsIntByCond));
stdSeg=zeros(size(DurationsSegByCond));
meanAll=zeros(size(DurationsIntByCond));
stdAll=zeros(size(DurationsIntByCond));

%%  
normIdx=[3,4];
clear allDurations
for i=1:size(DurationsIntByCond,1)
    for j=1:size(DurationsIntByCond,2)
allDurations{i,j}=[];
    end
end

for i=1:size(DurationsIntByCond,1)
    for j=1:size(DurationsIntByCond,2)
    meanInt(i,j)=mean(DurationsIntByCond{i,j});
    meanSeg(i,j)=mean(DurationsSegByCond{i,j});
    meanAll(i,j)=mean([DurationsIntByCond{i,j};DurationsSegByCond{i,j}]);
    stdInt(i,j)=std(DurationsIntByCond{i,j});
    stdSeg(i,j)=std(DurationsSegByCond{i,j});
    stdAll(i,j)=std([DurationsIntByCond{i,j};DurationsSegByCond{i,j}]);
    allDurations{i,j}=[allDurations{i,j};DurationsIntByCond{i,j};DurationsSegByCond{i,j}];
    end
end
% 
meanIntNorm=0;
meanSegNorm=0;

disp('Mean Int:')
disp(meanInt)
disp('Mean Seg:')
disp(meanSeg)

if ifnormalise
meanInt=meanInt/mean(allDurations{normIdx(1),normIdx(2)});
meanSeg=meanSeg/mean(allDurations{normIdx(1),normIdx(2)});
end
% return

%%  
figure(1);clf;hold on
set(gcf,'units','centimeters','position',[30,0,8,6])
set(gcf,'color','w');

figure(2);clf;hold on
set(gcf,'units','centimeters','position',[20,0,8,6])
set(gcf,'color','w')
NumIterations=size(SwitchData,3)

plotidx=[4];

phans1=[];
phans2=[];

colorsS2I={[0.1 0.1 0.5],[0.1 0.1 0.9],[0.2 0.5 0.8],[0.4 0.7 0.7]};
colorsI2S={[0.3 0 0.1],[1 0 0],[1 0.5 0],[0.7 0.1 0.3]};

binstr=[];nobins=20;
plotIdx=[1 3];legText={'$T_\textrm{mod}$ off','$T_\textrm{mod}=10\,s$'};
% plotIdx=[1:4];legText={'$T_\textrm{mod}$ off','$T_\textrm{mod}=5\,s$','$T_\textrm{mod}=10\,s$','$T_\textrm{mod}=20\,s$'};
% plotIdx=[2:3];legText={'$T_\textrm{mod}=5\,s$','$T_\textrm{mod}=10\,s$'};nobins=50;

for i=plotIdx%1:4%[1,3]
    t=tintCoarse;
    trsc=tintCoarse/rsctime;
    pfp=par1vals(i);
    if pfp==0
        pfp=10;
    end
    FF=1/pfp;% in s
    phiFcn=@(t)DFmin+0.5*(DFmax-DFmin)*(cos(2*pi*(t/100)*FF)+1);
    bins=DFmin:DFmax;
    SwtimesTmp=SwTimesByCond{i};
    SwtimesSegToInt=SwtimesTmp(SwtimesTmp(:,2)==1,1);
    SwtimesIntToSeg=SwtimesTmp(SwtimesTmp(:,2)==-1,1);
    
    h=2*pi/nobins;
    bins=h/2:h:2*pi-h/2;
    
    figure(1);
    if nobins>20;
        plot([pi/2,pi/2+1/10*2*pi],[0.49,0.49],'-','color',colorsS2I{2})
        plot([pi/2,pi/2+1/20*2*pi],[0.48,0.48],'-','color',colorsS2I{3})
    end
    PhaseStoI=mod(SwtimesSegToInt,pfp)/pfp*2*pi;
    [nStoI]=histnorm(PhaseStoI,bins);
    phans1(i)=plot(bins,nStoI,'color',colorsS2I{i});
    plot(linspace(0,2*pi),0.25+0.25*sin(linspace(0,2*pi)+pi/2),'color',grey)
    plot([pi pi],[0,0.5],'--','color',grey,'linewidth',1.5);
    set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
    if i==3
    set(gca,'XTickLabels',{'${\pi}/{2}$','$\pi$', '${3\pi}/{2}$', '0', '${\pi}/{2}$'})
    else
    end
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
    if i==plotIdx(end)
    lh=legend([phans1(plotIdx)],legText,'box','off','location','northeast');
    set(lh,'interpreter','latex','fontsize',10)
    end
    set(gca,'ylim',[0,0.5])    
     morestr='';
    if nobins>20
        morestr='-morebins';
    end
    svname=[plotdir,'Fig1CD-StoI',morestr,'-',num2str(length(plotIdx)),'-pfp-rsr2015-',num2str(period)];
    if ifexport && i==plotIdx(end)
        export_fig([svname,'.pdf'],'-transparent');
    end
    
    figure(2);
     if nobins>20;
        plot([pi/2,pi/2+1/10*2*pi],[0.49,0.49],'-','color',colorsI2S{2})
        plot([pi/2,pi/2+1/20*2*pi],[0.48,0.48],'-','color',colorsI2S{3})
    end
    
    PhaseItoS=mod(SwtimesIntToSeg-pfp/2,pfp)/pfp*2*pi;
    [nItoS]=histnorm(PhaseItoS,bins);
    phans2(i)=plot(bins,nItoS,'color',colorsI2S{i});
    plot(linspace(0,2*pi),0.25+0.25*cos(linspace(0,2*pi)+pi),'color',grey)
    plot([pi pi],[0,0.5],'-','color',grey,'linewidth',1.5);
    set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])    
    set(gca,'XTickLabels',{'${3\pi}/{2}$', '0', '${\pi}/{2}$','$\pi$','${3\pi}/{2}$'})
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
    if i==plotIdx(end)
    lh=legend([phans2(plotIdx)],legText,'box','off','location','northeast');
    set(lh,'interpreter','latex','fontsize',10)
    end
    set(gca,'ylim',[0,0.5])
    morestr='';
    if nobins>20
        morestr='-morebins';
    end
    svname=[plotdir,'Fig1CD-ItoS',morestr,'-',num2str(length(plotIdx)),'-pfp-rsr2015-',num2str(period)];
    if ifexport && i==plotIdx(end)
        export_fig([svname,'.pdf'],'-transparent');
    end
end
 
end