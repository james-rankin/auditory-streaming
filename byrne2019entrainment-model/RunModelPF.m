function RunModelPF(bfflag,numit_in)
%% 
% Code to accompany the paper Byrne, Rinzel and Rankin (2019)
% Entrainment of stream segregation in a dynamic environment
% Contact: james.rankin@gmail.com
% If you use or adapt this code acknowledge us by citing our paper
%
% To produce Fig 1B:
% RunModelPF
% or 
% RunModelPF('single')
%
% To produce quick  (less 'trials') version of Fig 1C--D:
% RunModelPF('pfp')
% or the versions exactly as in our paper
% RunModelPF('pfp',320)
%
% If you don't have a Parallel Computing Toolbox license replace 'parfor'
% with 'for'

ifexport=0; % don't save figure

% plotflag='short'; % plot first 20s
plotflag='long'; % plot full 4 minutes

if nargin==0;bfflag='single';end

%%  forcing parameters
DFmin=2;
DFmax=8;


%% set up brute force or single computation
switch bfflag
    case 'pfp' 
        NumIterations=36;
        par1name='PFperiod';
        par1vals=[0,5,10,20];
        N=length(par1vals);
        par2name='null';
        M=1;par2vals=0; 
    case 'single'
        NumIterations=1;
        par1name='PFperiod';
        par1vals=[10]; % set Tmod (0 for no forcing)
        N=length(par1vals);
        par2name='null';
        M=1;par2vals=0; 
        rng(1020) % noise instantiation used in paper
    otherwise
end


%% Intrinsic parameters
taua=140; % 1.4s;  adaptation timescale
taunm=7; % 70ms; recurrent excitation timescale
taud=300; % 3s; synaptic depression timescale; unused when synd=0;
betai=0.3; % inhibition strength
g=0.11; % adaptation strength
kf=12; tf=0.2; % for firing rate function
betae=0.85;% recurrent excitation strength
synd=0.25;% synaptic depression strength


%% input and noise parameters
PeriodFlag='ABA'; % AB or ABA or sync or A
PR=1000/125;
alpha1=1.5; % 15ms rise time 
alpha2=8.25; % 82.5ms rise time 
per=(100/PR)*2; % parameters for input 'impulse' function

taux=10; % noise timescale 100ms

gamma=0.075;% noise strength
tphi=4.25;% spatial decay constant for input over tonotopy 
Iamp=0.47;% input amplitude

%% Initial conditions
tstep=0.5;tfin=24000; % time discretisation
u0=zeros(12,1);u0(10:12)=1;
u0(1)=0.2+0.5*rand(1);
u0(7:9)=u0(1:3)+0.3;


%% Firing rate function, connectivity, impulse function

% Input dependence on DF
InputDF=@(phi) exp(-phi/tphi);
f=@(u)1/(1+exp((tf-u)*kf));
heaviside=@(x) (x==0).*0.5 + (x>0).*1;

%% Define Input
per_b=2*per;
switch PeriodFlag
    case 'AB'
        phA=@(t)t-floor(t/per)*per;% phase for A
        phB=@(t)t-floor((t-per/2)/per)*per-per/2;% phase for B same period
    case 'ABA'
        phA=@(t)t-floor(t/per)*per;
        phB=@(t)t-floor((t-per/2)/per_b)*per_b-per/2;% phase for B half period
end

DILong=@(t)...
            heaviside(t).*exp(2)./(alpha1).^2.*t.^2 .* exp(-2/(alpha1).*abs(t))+...
            1/6*heaviside(t).*exp(2)./(alpha2).^2.*t.^2 .* exp(-2/(alpha2).*abs(t));        

DI=DILong;
DIA=@(t)DI(phA(t));
DIB=@(t)DI(phB(t));
%% override parameter values 
if nargin>1; if ~isempty(numit_in);NumIterations=numit_in;end;end
        
%% Start simulation loop(s)
SwitchData=struct([]);

% set where to save data/plots based on which pc this runs on
[~,pcname]=system('hostname');
if strcmp(strtrim(pcname),'jrankin-dt')
svdir='/home/jr3830/Insync/proj/auditory/codes/auditory-streaming/pfdata/';
plotdir='/home/jr3830/Insync/RankinExpt/DFMod/HearRes/makefigs/';
else; svdir=pwd;plotdir=pwd;
end

svname=['SDIB_',bfflag,num2str(N),'x',num2str(M),'x',num2str(NumIterations),'_',...
    par1name,num2str(par1vals(1)),'_',num2str(par1vals(end)),'_',...
    par2name,num2str(par2vals(1)),'_',num2str(par2vals(end)),'_',...
    'dfpm',num2str((DFmax-DFmin)/2),'_',...
    'ts',num2str(tstep),'_',...
    'g',num2str(g),'_',...
    'gam',num2str(gamma)...
    ];
svname=strrep(svname,'.','p');
svname=strrep(svname,'-','m');
disp(svname)

for i=1:length(par1vals)
    eval([par1name,'=',num2str(par1vals(i)),';']);
    
    for j=1:length(par2vals)
        eval([par2name,'=',num2str(par2vals(j)),';']);
        disp(['par1 value ',num2str(i),' of ',...
            num2str(length(par1vals)),'  ',par1name,'=',num2str(par1vals(i))])
        disp(['par2 value ',num2str(j),' of ',...
            num2str(length(par2vals)),'  ',par2name,'=',num2str(par2vals(j))])
        
        opts=odeset('reltol',10e-5,'maxstep',2.5);
        
        yn=[];% noise will be generated in sde_dogear_loadinp... if yn is empty

        tint=0:tstep:tfin;
        t=tint;
        rsctime=100;trsc=t/rsctime;
        
        if PFperiod>0
            FF=1/PFperiod;
            phiFcn=@(t)DFmin+0.5*(DFmax-DFmin)*(cos(2*pi*(t/100)*FF)+1);
        else % no modulation input case
            phiFcn=@(t)(DFmin+DFmax)/2;
        end
        
        % main 3 equations
        SargInLoad=@(t,u)[betae*u(7)*u(10)-betai*u(1)-betai*(u(2)+u(3))-g*u(4);
            betae*u(8)*u(11)-betai*u(2)-2*betai*u(1)-betai*u(3)-g*u(5);
            betae*u(9)*u(12)-betai*u(3)-2*betai*u(1)-betai*u(2)-g*u(6)];
        % 9 timescale eqs defined in SDERun to allow input pre-allocation 
        
        % pre-allocating the inputs gives significant speed up in SDERun
        InVals=zeros(length(tint),3);
        InVals(:,1)=InputDF(phiFcn(t)/2).*Iamp.*(DIA(t)+DIB(t));
        InVals(:,2)=Iamp.*DIA(t)+InputDF(phiFcn(t)).*Iamp.*DIB(t);
        InVals(:,3)=Iamp.*DIB(t)+InputDF(phiFcn(t)).*Iamp.*DIA(t);
        
        ycell=cell(NumIterations,1);
        tic % stopwatch whilst running the model
        parfor k=1:NumIterations
            [ycell{k}]=SDERun(tint,f,SargInLoad,u0,taua,taunm,taux,taud,synd,abs(gamma),InVals,yn);
        end
        toc
        
        tic % stopwatch for post-processing output for saving
        tintCoarse=0:tstep*10:tfin;
        MinWindow=0.25;ifplot=0;
        SlidingWidthT=50;SlidingWidthIdx=SlidingWidthT/tstep;
        parfor k=1:NumIterations
            ycell{k};
            uFilt{k}=slidingavg(ycell{k}(:,1),SlidingWidthIdx);
            vFilt{k}=slidingavg((ycell{k}(:,2)+ycell{k}(:,3))/2,SlidingWidthIdx);
            [NumberOfSwitches{k},First{k},Durations{k},DurationsInt{k},DurationsSeg{k}...
                ,PropIntegrated{k},SwitchTimes{k},SwDirections{k}]=...
                ProcSwTimes(trsc,MinWindow,uFilt{k},vFilt{k},ifplot);
            
            if k==1
                disp(['Iteration ',num2str(k),' of ',...
                    num2str(NumIterations),'; Number of switches: ',num2str(NumberOfSwitches{k})])
            end
            uFiltCoarse{k}=interp1(tint,uFilt{k},tintCoarse);
            vFiltCoarse{k}=interp1(tint,vFilt{k},tintCoarse);
        end
        
        for k=1:NumIterations
            PropIntegratedMtx(i,j,k)=PropIntegrated{k};
            SepMtx(i,j,k)=mean(abs(uFilt{k}-vFilt{k}));
            SwitchData(i,j,k).uFiltCoarse=uFiltCoarse{k};
            SwitchData(i,j,k).vFiltCoarse=vFiltCoarse{k};
            SwitchData(i,j,k).First=First{k};
            SwitchData(i,j,k).Durations=Durations{k};
            SwitchData(i,j,k).DurationsInt=DurationsInt{k};
            SwitchData(i,j,k).DurationsSeg=DurationsSeg{k};
            SwitchData(i,j,k).SwitchTimes=SwitchTimes{k};
            SwitchData(i,j,k).SwDirections=SwDirections{k};
            SwitchData(i,j,k).Ipamp=Iamp;
            if NumberOfSwitches{k}==0
                NumberOfSwitchesMtx(i,j,k)=sign(uFilt{k}(end)-vFilt{k}(end));
            else
                NumberOfSwitchesMtx(i,j,k)=NumberOfSwitches{k};
            end
        end
        toc
    end
end


if strcmp(bfflag,'pfp')
    datestring=datestr(floor(now));
    datestring=strrep(datestring,'-','');
    timestring=datestr(rem(now,1),13);
    timestring=strrep(timestring,':','');
    clear ycell y ifplot ifexport
    
    savename=[svdir,svname,'_',datestring,'_',timestring,'.mat'];
    save(savename);
    disp(savename);
  
    ProcPFData(savename)% plot the saved data from many simulations
    return% don't run plotting commands below
end
%% Everything below here is plotting for bfflag='single' case
make_colors
y=ycell{k};

%% define colours and linestyles
make_colors;
docolour=1;
fontsize=12;
fsax=8;
fsleg=10;
linewidth=1.5;
panelW=9;
panelH=6;
switch docolour
    case 1
        clrs={blue,red,green};clrs_per={blue,red};
        linewi={linewidth,linewidth,linewidth};linest={'-','--','-'};
    case 0
        clrs={grey,dark_grey,black};clrs_per={light_grey,dark_grey};
        linewi={linewidth,linewidth,linewidth};linest={'-','--','-'};
end
 
%% Post-process
MinWindow=0.55;
t0=0.5;
Tmin=.5;
Tmax=240;
SlidingWidthT=50;SlidingWidthIdx=SlidingWidthT/tstep;
uFilt=slidingavg(y(:,1),SlidingWidthIdx);
vFilt=slidingavg((y(:,2)+y(:,3))/2,SlidingWidthIdx);
[NumberOfSwitches,First,Durations,DurationsInt,DurationsSeg,PropIntegrated,...
    SwitchTimes,SwDirections]=...
    ProcSwTimesT0(trsc,t0,MinWindow,uFilt,vFilt,0);
DurationsInt(DurationsInt<Tmin)=[];
DurationsInt(DurationsInt>Tmax)=[];
DurationsSeg(DurationsSeg<Tmin)=[];
DurationsSeg(DurationsSeg>Tmax)=[];
tl0=15;
tlocal=25;
tlocalidx=find(trsc>tlocal,1);
tlIdx=tl0:tlocalidx;
disp(['Number of switches:',num2str(NumberOfSwitches)])
disp(['Proportion Integrated:',num2str(PropIntegrated)]);
disp(['Mean Duration:',num2str(mean(Durations))]);
disp(['Mean Integrated:',num2str(mean(DurationsInt))]);
disp(['Mean Segregated:',num2str(mean(DurationsSeg))]);
First
switch plotflag
    case 'short'
        xend=20;
        timelab=16;
    case 'long'
        xend=240;
        timelab=165;
end

%%
figure(2);clf;hold on
 

h1=plot(trsc,y(:,1),linest{1},'color',clrs{1},'linewidth',linewi{1});
h2=plot(trsc,y(:,2),linest{2},'color',clrs{2},'linewidth',linewi{2});
h3=plot(trsc,y(:,3),linest{3},'color',clrs{3},'linewidth',linewi{3});
plot([trsc(1),trsc(end)],[tf,tf],'k--','linewidth',linewidth)
ylabhan=ylabel('Activity','Interpreter','latex','fontname','helvetica','fontsize',fontsize);
for i=1:length(SwitchTimes)
   plot([SwitchTimes(i),SwitchTimes(i)],[-.2,1.1],'k-','linewidth',linewidth)
end

plot([tl0],[-2 2],'k-')

set(gca,'xlim',[0 xend],'ylim',[-.05 1.05]);
set(gca,'fontname','helvetica','fontsize',fsax);

set(gca,'linewidth',linewidth);set(gcf,'color','w');box on

switch plotflag
    case 'short'
        text(-1.3,1.0,'\textbf{A}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.15,'time (s)');
        set(than,'Interpreter','latex','fontname','helvetica','fontsize',fontsize);
        set(gca,'xtick',[0:5:30]);
    case 'long'
        text(-16,1.,'\textbf{*}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.15,'time (s)');
        set(than,'Interpreter','latex','fontname','helvetica','fontsize',fontsize);
end
lh=legend([h1,h2,h3],'$r_{AB}$','$r_{A}$','$r_{B}$');
set(lh,'Interpreter','latex');

p1windows=[];p2windows=[];
SwitchTimes=[0;SwitchTimes'];
SwDirections=[-SwDirections(1);SwDirections];


for i=1:length(SwitchTimes)-1
   if SwDirections(i)==1
       p1windows=[p1windows;SwitchTimes(i),SwitchTimes(i+1)];
   else
       p2windows=[p2windows;SwitchTimes(i),SwitchTimes(i+1)];
   end
end

%%
figure(12);clf;hold on
box on
xmax=90;
PresentationLengthInSeconds=trsc(end);
set(gca,'xlim',[0,PresentationLengthInSeconds]);
set(gcf,'units','centimeters','position',[3,6,16,4]);
set(gca,'xlim',[0 xmax],'ylim',[-.05 1.75]);
set(gca,'activepositionproperty','position');
set(gca,'position',[0.09,0.23,0.90,0.7]);
set(gca,'xtick',0:20:80)

make_colors
clrs_per={clrs{1},clrs{2}};
for i=1:size(p1windows,1)
    XY=[p1windows(i,1),0.5;
        p1windows(i,1),1;
        p1windows(i,2),1;
        p1windows(i,2),0.5;];
   patch(XY(:,1),XY(:,2),clrs_per{1});
end

for i=1:size(p2windows,1)
    XY=[p2windows(i,1),0;
       p2windows(i,1),.5;
       p2windows(i,2),.5;
       p2windows(i,2),0;];
   patch(XY(:,1),XY(:,2),clrs_per{2});
end
set(gca,'TickLabelInterpreter','latex')
if PFperiod==10
set(gca,'ytick',[0.25,0.75,1.1,1.6],'yticklabel',{'Seg','Int','DF=2','DF=8'}); 
else
set(gca,'ytick',[0.25,0.75,1.35],'yticklabel',{'Seg','Int','DF=5'}); 
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)

if PFperiod>0
FFUpIdx=find(diff([0,phiFcn(t)])>=0);
FFDnIdx=find(diff([0,phiFcn(t)])<0);

FFUCoraseIdx=1:10:length(FFUpIdx);
FFDCoraseIdx=1:10:length(FFDnIdx);
plot(trsc,1.1+.5*(phiFcn(t)-DFmin)/(DFmax-DFmin),'color',grey,'linewidth',1.5)

TMod=PFperiod;tend=xmax;
fMaxTimes=(0:TMod:tend);
fMinTimes=(TMod/2:TMod:tend);

for i=1:length(fMaxTimes)
    plot([fMaxTimes(i),fMaxTimes(i)],[1.1,1.6],'-','color',grey,'linewidth',1);
end

for i=1:length(fMinTimes)
    plot([fMinTimes(i),fMinTimes(i)],[1.1,1.6],'--','color',grey,'linewidth',1);
end
end
than=text(65,-0.33,'Time (s)');

set(than,'Interpreter','latex','fontname','helvetica','fontsize',12);

svname=[plotdir,'Fig1B-new-percept',...
    '-pfp',num2str(PFperiod)];
svname=strrep(svname,'.','p');
set(gcf,'renderer','painters')
if ifexport
% export_fig([svname,'.pdf'],'-transparent');
saveas(gcf,[svname,'.svg'])
end

if strcmp(plotflag,'long')
    return
end

%%
figure(4);clf;hold on
set(gca,'xlim',[0 xend],'ylim',[-.05 1.05]);
set(gca,'fontname','helvetica','fontsize',fsax);

set(gca,'linewidth',linewidth);set(gcf,'color','w');box on
h4=plot(trsc,y(:,4),linest{1},'color',clrs{1},'linewidth',linewi{1});
h5=plot(trsc,y(:,5),linest{2},'color',clrs{2},'linewidth',linewi{2});
h6=plot(trsc,y(:,6),linest{3},'color',clrs{3},'linewidth',linewi{3});

ylabhan=ylabel('Activity','Interpreter','latex','fontname','helvetica','fontsize',fontsize);

for i=1:length(SwitchTimes)
   plot([SwitchTimes(i),SwitchTimes(i)],[-2 2],'k-','linewidth',linewidth)
end
lh2=legend([h4,h5,h6],{'$a_{AB}$','$a_{A}$','$a_{B}$'});
set(lh2,'Interpreter','latex');

switch plotflag
    case 'short'
        text(-1.3,1.0,'\textbf{C}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.15,'time (s)');
        set(gca,'xtick',[0:5:30]);
    case 'long'
        text(-16,1.0,'\textbf{F}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.15,'time (s)');
end
set(than,'Interpreter','latex','fontname','helvetica','fontsize',fontsize);

%%
figure(13);clf;hold on
set(gca,'xlim',[0 xend],'ylim',[-.05 1.05]);
set(gca,'fontname','helvetica','fontsize',fsax);

set(gca,'linewidth',linewidth);set(gcf,'color','w');box on
h7=plot(trsc,y(:,7),linest{1},'color',clrs{1},'linewidth',linewi{1});
h8=plot(trsc,y(:,8),linest{2},'color',clrs{2},'linewidth',linewi{2});
h9=plot(trsc,y(:,9),linest{3},'color',clrs{3},'linewidth',linewi{3});

ylabhan=ylabel('Activity','Interpreter','latex','fontname','helvetica','fontsize',fontsize);

for i=1:length(SwitchTimes)
   hst=plot([SwitchTimes(i),SwitchTimes(i)],[-2 2],'k-','linewidth',linewidth);
end
lh3=legend([h7,h8,h9],{'$e_{AB}$','$e_{A}$','$e_{B}$'});
set(lh3,'Interpreter','latex');


switch plotflag
    case 'short'
        text(-1.3,1.0,'\textbf{B}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.15,'time (s)');
        set(gca,'xtick',[0:5:30]);
    case 'long'
        text(-16,1.0,'\textbf{E}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.15,'time (s)');
end
set(than,'Interpreter','latex','fontname','helvetica','fontsize',fontsize);


return



