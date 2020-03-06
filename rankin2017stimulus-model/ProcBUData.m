function ProcBUData(fname)


if nargin==0;fname= 'buildup_data.mat';end

plotdir=pwd;
datadir=pwd;

load(fname);
make_colors

vfac=2;
ifexport=0;% requires export_fig; download from Matlab Central
par2int=1;
MinWindow=0.05;
T0=0.5;
Tmin=.5;
Tmax=60;
ifplot=0;

NumRepeats=size(SwitchData,3);
NumConditions=size(SwitchData,1);
NumCases=size(SwitchData,2);
PropIntMtx=zeros(NumRepeats,NumConditions);
DurationsIntByCond={};
DurationsSegByCond={};
parvals=par1vals;
DFIdx=1:length(par1vals);
propFirstInt=zeros(NumConditions,1);
meanFirstInt=zeros(NumConditions,1);
meanFirstSeg=zeros(NumConditions,1);

timeVals=[1.2,2.8,4];

spcount=0;
binwidth=0.1;
t0=0;
tf=7;
bincs=[t0:binwidth:tf-binwidth;t0+binwidth:binwidth:tf]';

sliceCell={};labstrsP1={};
buildupCell={};
make_colors
clrs={purple,green,grey};
MtxFstInt=zeros(length(par1vals),length(par2vals),NumIterations);
MtxScdSeg=zeros(length(par1vals),length(par2vals),NumIterations);
for p2idx=1:length(par2vals)
BuildUpMtx=zeros(length(par1vals),size(bincs,1));
for i=1:length(par1vals)
    t=tintCoarse;
    trsc=tintCoarse/rsctime;
    
    PropIntBinned=zeros(size(bincs,1),1);
    TimesIntToSeg=[];
    TimesSegToInt=[];
    spcount=spcount+1;
    
    for j=1:size(bincs,1)
        IdxBinned=find(trsc>=bincs(j,1)&trsc<bincs(j,2));
        PropIntCount=0;
        for k=1:NumIterations
            uFilt=SwitchData(i,(p2idx),k).uFiltCoarse;
            vFilt=SwitchData(i,(p2idx),k).vFiltCoarse;
            if j==1
            [NumberOfSwitches,First,Durations,DurationsInt,DurationsSeg,PropIntegrated,...
                SwitchTimes,SwDirections]=...
                ProcSwTimesT0(trsc,T0,MinWindow,uFilt,vFilt,ifplot);
            if numel(Durations>=2) && First(2)==1
                MtxFstInt(i,p2idx,k)=Durations(1);
                MtxScdSeg(i,p2idx,k)=Durations(2);
            end
            end
            PropIntCount=PropIntCount+numel(find(vfac*vFilt(IdxBinned)>uFilt(IdxBinned)));        
        end
        PropIntBinned(j)=PropIntCount/numel(IdxBinned)/NumIterations;
    end
    BuildUpMtx(i,:)=PropIntBinned;
end

bincenters=0.15+(bincs(:,2)+bincs(:,1))/2;
sliceCell{p2idx}=sliceMtx;
buildupCell{p2idx}=BuildUpMtx;
end
MtxFstInt(MtxFstInt==0)=NaN;
MtxScdSeg(MtxScdSeg==0)=NaN;


%%

make_colors
figure(1);clf;hold on
set(gca,'ylim',[0,1],'xlim',[0,tfin/rsctime])
df0Mtx=buildupCell{1};
% clrs={black,grey,light_grey};
clrs={green,red,blue,black}; 
clrs={light_grey,grey,black};
plot(bincenters,df0Mtx(1,:),'--','color',clrs{1});
plot(bincenters,df0Mtx(2,:),'--','color',clrs{2});
l1=plot(bincenters,df0Mtx(3,:),'--','color',clrs{3});
ylabel('Proportion Segregated')
xlabel('time (s)')
df0p8Mtx=buildupCell{2};
plot(bincenters,df0p8Mtx(1,:),'-','color',clrs{1},'linewidth',1.5);
plot(bincenters,df0p8Mtx(2,:),'-','color',clrs{2},'linewidth',1.5);
l2=plot(bincenters,df0p8Mtx(3,:),'-','color',clrs{3},'linewidth',1.5);

text(bincenters(end)-0.7,df0p8Mtx(1,end)-0.07,'DF=4');
text(bincenters(end)-0.7,df0p8Mtx(2,end)-0.06,'DF=7');
text(bincenters(end)-0.7,df0p8Mtx(3,end)+0.05,'DF=10');

legend([l1,l2],{'Input static','Input adapting'},'location','southeast','box','off')

tripStr={'3trips','7trips','10trips'};
for tvIdx=1:length(timeValsTmp);
    plot([timeValsTmp(tvIdx),timeValsTmp(tvIdx)],[0,1],'k--','linewidth',1);
    text(timeValsTmp(tvIdx)-0.3,1.02,tripStr(tvIdx))
end
% th=title('Model: build-up function','fontsize',16,'fontweight','normal');
% set(th,'position',[3.5,1.05,0]);


set(gcf,'units','centimeters','position',[0,0,10,4])
set(gca,'position',[0.1200    0.1550    0.83    0.79])

fig=gcf;
set(findall(fig,'-property','Fontname'),'Fontname','helvetica')
set(findall(fig,'-property','FontSize'),'FontSize',8)

if ifexport
    export_fig([plotdir,'model-build-up-cns.pdf']);
end

%%

clrs={light_grey,grey,black};
figure(2);clf;hold on
set(gca,'ylim',[0,1],'xlim',[par1vals(1)-0.5,par1vals(end)+0.5])

df0slMtx=sliceCell{2};
plot(par1vals,df0slMtx(:,1),'.-','color',clrs{1},'linewidth',1,'markersize',18);
plot(par1vals,df0slMtx(:,2),'.-','color',clrs{2},'linewidth',1,'markersize',18);
l2=plot(par1vals,df0slMtx(:,3),'.-','color',clrs{3},'linewidth',1,'markersize',18);


tripStr={'3trips (1.2s)','7trips (2.8s)','10trips (4s)'};
yshift=[-.35,-0.20,0.06];
for tvIdx=1:length(timeValsTmp);
    text(7.8,df0slMtx(end,tvIdx)+yshift(tvIdx),tripStr(tvIdx),'color',clrs{tvIdx})
end
set(gca,'xtick',par1vals)


rawdir='/home/jr3830/proj/auditory/codes/bdisdev/rawexpdata/';
rawname='gaps_pseg_exp1.mat';

load([rawdir,rawname]);
condIdx=1:3;
NumSubj=8;
NumDF=3;
offset=[0.15,0.15,0.15]
for i=1:length(condIdx)
   gapDataTmp=reshape(gaps_pseg_exp1(condIdx(i),:,:),NumDF,NumSubj)/100;
   errorbar(par1vals+offset(i),mean(gapDataTmp'),std(gapDataTmp')/sqrt(NumSubj),'--d'...
       ,'color',clrs{i},'markerfacecolor',clrs{i},'linewidth',1,'markersize',5)
end

set(gcf,'units','centimeters','position',[0,0,5.5,4])
set(gca,'position',[0.1400    0.1550    0.83    0.79])
fig=gcf;
set(findall(fig,'-property','Fontname'),'Fontname','helvetica')
set(findall(fig,'-property','FontSize'),'FontSize',8)

if ifexport
    export_fig([plotdir,'model-snapshot.pdf']);
end
