function RunModelNM
%% 
% Code to accompany the paper Rankin, Sussman and Rinzel (2015)
% Neuromechanistic Model of Auditory Bistability
% PLoS Computational Biology DOI:10.1371/journal.pcbi.1004555
% Contact: james.rankin@gmail.com
% 
% If you use or adapt this code acknowledge us by citing our paper
%
% To produce Fig 4.
% RunModelNM

clear,close all
make_colors
%% Change this string to reproduce fig 4 or a simulation for model version 
% - with fixed nmda excitation and local inhibition (efix)
% - with dynamic nmda exication (with synaptic depression) and global inhibition (edyn)
caseflag='fig4';DF=5;
% caseflag='efix';phi=5; % DF value
% caseflag='edyn';phi=5; % DF value

plotflag='short'; % plot first 20s
% plotflag='long'; % plot full 4 minutes
%% Intrinsic parameters
taua=140; % 1.4s;  adaptation timescale
taunm=7; % 70ms; recurrent excitation timescale
taud=300; % 3s; synaptic depression timescale; unused when synd=0;
betai=0.3; % inhibition strength
g=0.065; % adaptation strength
kf=12; tf=0.2; % for firing rate function

% Parameters changed in paper
switch caseflag
    case {'fig4','efix'}
        betanm=0.4;
        synd=0.0;
    case 'edyn'
        betanm=0.55;
        synd=0.25;
    otherwise
        error('specify your own parameter values in otherwise statements')
end

%% input and noise parameters
PR=8;
alpha1=1.5; % 15ms rise time 
alpha2=8.25; % 82.5ms rise time 
per=(100/PR)*2; % some parameters for input 'impulse' function
Ibase=0; %  not used
taux=10; % noise timescale 100ms


% Parameters changed in paper
switch caseflag
    case 'fig4'
        tphi=5;% reported as 10 in paper
        A=0.05/2; % stricly speaking this was done with A=0.5 but it makes no difference at DF=5;
        B=5; % reported as 10 in the paper 
        Ipamp=0.525;
    case 'efix'
        tphi=4;% reported as 8 in paper
        A=0.05; % excitation is local only  
        Ipamp=0.525;  
        B=6; % reported as 10 in the paper 
    case 'edyn'
        tphi=4.25;% reported as 8.5 in paper
        A=0.05; % excitation is local only  
        Ipamp=0.47;  
        B=100; % inhibition is global
    otherwise
        
        
end



%% noise pars
xidx=[1,2,3]; % add noise to 1st 3 eqns
loadnoise=1;
switch caseflag
    case 'fig4'
        sigma=0.08;%%%%%%%%%%%%%%%%%
        loadnoise=1;
    case 'efix'
        sigma=0.08; % reported as 0.075 in paper
        loadnoise=0;
    case 'edyn'
        sigma=0.075;
    otherwise
        
end

%% Initial conditions
tstep=0.5;tfin=24000; % time discretisation
u0=zeros(12,1);u0(10:12)=1;
switch caseflag
    case 'fig4'
        u0(1)=0.630418733103864;
        u0(7)=u0(1);
    case {'efix','edyn'}
        u0(1)=0.2+0.5*rand(1);
        u0(7:9)=u0(1:3)+0.3;
    otherwise
        
end


%% Firing rate function, connectivity, impulse function
% Functions controlling DF-dependent Ex/Inh
Ce=betai+betanm;Ci=betai;
CE=@(DF) Ce*exp(-(DF.^2)/(2*A^2));
CI=@(DF) Ci*exp(-(DF.^2)/(2*B^2));
conn=@(phi) Ce*exp(-(phi.^2)/(2*A^2))-Ci*exp(-(phi.^2)/(2*B^2));

% Input dependence on DF
InputDF=@(phi) exp(-phi/tphi);
f=@(u)1/(1+exp((tf-u)*kf));
heaviside=@(x) (x==0).*0.5 + (x>0).*1;

%% Define Input
per_b=2*per;

%     case 'ABA'
phA=@(t)t-floor(t/per)*per;
phB=@(t)t-floor((t-per/2)/per_b)*per_b-per/2;% phase for B half period


DILong=@(t)...
            heaviside(t).*exp(2)./(alpha1).^2.*t.^2 .* exp(-2/(alpha1).*abs(t))+...
            1/6*heaviside(t).*exp(2)./(alpha2).^2.*t.^2 .* exp(-2/(alpha2).*abs(t)); 
DI=DILong;
DIA=@(t)DI(phA(t));
DIB=@(t)DI(phB(t));


%% Start simulation
SargInLoad=@(t,u)[(CE(0)*u(7)*u(10)+CE(DF/2)*(u(8)+u(9))...Excitation NMDA rAB
    -CI(0)*u(1)-(1)*CI(DF/2)*(u(2)+u(3))-g*u(4));
    (CE(0)*(u(8)*u(11))+CE(DF/2)*u(7)+CE(DF)*u(9)...Excitation NMDA rA
    -CI(0)*u(2)-(1)*2*CI(DF/2)*u(1)-CI(DF)*u(3)-g*u(5));
    (CE(0)*(u(9)*u(12))+CE(DF/2)*u(7)+CE(DF)*u(8)...Excitation NMDA rB
    -CI(0)*u(3)-(1)*2*CI(DF/2)*u(1)-CI(DF)*u(2)-g*u(6))];

opts=odeset('reltol',10e-5,'maxstep',2.5);

yn=[];
if loadnoise
    load('fig4noise.mat','ynloaded');
    yn=ynloaded;
end
tint=0:tstep:tfin;
t=tint;

InVals=zeros(length(tint),3);
InVals(:,1)=InputDF(DF/2)*(2*Ibase+Ipamp*(DIA(t)+DIB(t)));
InVals(:,2)=(Ibase+Ipamp*DIA(t))+InputDF(DF)*(Ibase+Ipamp*DIB(t));
InVals(:,3)=(Ibase+Ipamp*DIB(t))+InputDF(DF)*(Ibase+Ipamp*DIA(t));

tic
[y,yn]=SDERun(tint,f,SargInLoad,u0,taua,taunm,taux,taud,synd,abs(sigma),InVals,yn);
toc

rsctime=100;trsc=t/rsctime;

%% Everything below here is plotting 
% define colours and linestyles
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

% Make 30s TH main variables


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
% plot([tlocal],[-2 2],'k-')
set(gca,'xlim',[0 xend],'ylim',[-.05 1.05]);
set(gca,'fontname','helvetica','fontsize',fsax);
% set(gca,'activepositionproperty','position');
% set(gca,'position',[0.07,0.14,0.91,0.84]);
set(gca,'linewidth',linewidth);set(gcf,'color','w');box on
% set(gca,'ylim',[30,51])
% set(gcf,'units','centimeters','position',[3,12,18,4]);
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


% ThrowIdx=find(Durations<MinWindow)
% SwitchTimes(ThrowIdx)=[];
% disp(['Throwing away ',num2str(numel(ThrowIdx)),' short durations']);
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
PresentationLengthInSeconds=trsc(end);
% set(gca,'xlim',[0,xend]);
% set(gcf,'units','centimeters','position',[3,6,18,3]);
box on
set(gca,'xlim',[0 xend],'ylim',[-.05 1.05]);
set(gca,'xtick',[0:40:240])
set(gca,'fontname','helvetica','fontsize',fsax);
% set(gca,'activepositionproperty','position');
% set(gca,'position',[0.07,0.23,0.91,0.7]);
set(gca,'linewidth',linewidth);set(gcf,'color','w');box on

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
set(gca,'ytick',[0.25,0.75],'yticklabel',{'Seg','Int'});

switch plotflag
    case 'short'
        text(-1.3,1.05,'\textbf{D}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.28,'time (s)');
        set(gca,'xtick',[0:5:30]);
    case 'long'
        text(-16,1.05,'\textbf{E}','interpreter','latex','fontsize',fontsize)
        than=text(timelab,-0.28,'time (s)');
end
set(than,'Interpreter','latex','fontname','helvetica','fontsize',fontsize);


%%
figure(4);clf;hold on
set(gca,'xlim',[0 xend],'ylim',[-.05 1.05]);
set(gca,'fontname','helvetica','fontsize',fsax);
% set(gca,'activepositionproperty','position');
% set(gca,'position',[0.07,0.14,0.91,0.84]);
set(gca,'linewidth',linewidth);set(gcf,'color','w');box on
h4=plot(trsc,y(:,4),linest{1},'color',clrs{1},'linewidth',linewi{1});
h5=plot(trsc,y(:,5),linest{2},'color',clrs{2},'linewidth',linewi{2});
h6=plot(trsc,y(:,6),linest{3},'color',clrs{3},'linewidth',linewi{3});

ylabhan=ylabel('Activity','Interpreter','latex','fontname','helvetica','fontsize',fontsize);
% hi=plot(trsc,y(:,7),'b--',trsc,y(:,8),'r--',trsc,y(:,9),'g--');
% set(hi,'linewidth',2)
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

% set(gcf,'units','centimeters','position',[3,4,18,4]);

%%
figure(13);clf;hold on
set(gca,'xlim',[0 xend],'ylim',[-.05 1.05]);
set(gca,'fontname','helvetica','fontsize',fsax);
% set(gca,'activepositionproperty','position');
% set(gca,'position',[0.07,0.14,0.91,0.84]);
set(gca,'linewidth',linewidth);set(gcf,'color','w');box on
h7=plot(trsc,y(:,7),linest{1},'color',clrs{1},'linewidth',linewi{1});
h8=plot(trsc,y(:,8),linest{2},'color',clrs{2},'linewidth',linewi{2});
h9=plot(trsc,y(:,9),linest{3},'color',clrs{3},'linewidth',linewi{3});

ylabhan=ylabel('Activity','Interpreter','latex','fontname','helvetica','fontsize',fontsize);
% hi=plot(trsc,y(:,7),'b--',trsc,y(:,8),'r--',trsc,y(:,9),'g--');
% set(hi,'linewidth',2)
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

% set(gcf,'units','centimeters','position',[3,4,18,4]);


return

