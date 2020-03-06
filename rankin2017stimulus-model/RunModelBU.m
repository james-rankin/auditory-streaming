function RunModelBU(bfflag)
%% 
% Code to accompany the paper Rankin, Osborn Popp, Rinzel 2017
% Stimulus Pauses and Perturbations Differentially Delay or Promote the
% Segregation of Auditory Objects: Psychoacoustics and Modeling
% Frontiers in Neuroscience https://doi.org/10.3389/fnins.2017.00198
% Contact: james.rankin@gmail.com
% If you use or adapt this code acknowledge us by citing our paper
%
% To produce Fig 1c:
% RunModelBU
% or 
% RunModelPF('fig1c')
%
% To produce quick  data for Fig 1C--D:
% RunModelPF('fig1de') to produce buildup_data.mat
% plot with >> ProcBUData('buildupdata.mat')

ifexport=0; % don't save figure; requires export_fig; download from Matlab Central

if nargin==0;bfflag='fig1c';end
par1name='DF';
N=3;par1vals=[4,7,10]; 
simpleFilenameFlag=0;
switch bfflag
    case 'fig1de'
        par2name='DFFac';
        M=2;par2vals=[0,0.9];
        NumIterations=750;
        rngStart=2046;
        simpleFilenameFlag=1;
    case 'fig1c' % Note this does not produce the _exact_ figure from the paper
        N=1;par1vals=[4];
        par2name='null';
        M=1;par2vals=0;
        NumIterations=1;
        rngStart=2046+11;
end

%% Intrinsic parametersDF=5;
DF=5; % overwritten if par1name is 'DF'
taua=140;
taunm=7;
taud=300;

g=0.045;
kf=12; tf=0.2; % for firing rate function

synd=0.0;

betai=0.3;
betae=betai+0.35;
 %% input and noise parameters
PeriodFlag='ABAfix'; % AB or ABA or sync or A
PR=10;
dr=1.5;

per=(100/PR)*2;

taux=10;

gamma=0.0875;
tpConstExp=14.7;
tpSlopeExp=9;

Ibase=0.03;
IpConst=0.6;
IpSlope=0;
PpFrac=.5; % because B is at half PR of A

%% A1 adaptation parameters
DFFac=0.9;
MichFac=2.5;
phidec=50;

%% Firing rate function
f=@(u)1/(1+exp((tf-u)*kf));
heaviside=@(x) (x==0).*0.5 + (x>0).*1; 

%% Initial conditions and set up noise
tstep=0.5;tfin=700;
ICFac=3;
rng(rngStart)
icRnd=rand(3,NumIterations);
rngSeeds=rngStart:rngStart+NumIterations;

%% Preallocate data structures
NumberOfSwitchesMtx=zeros(length(par1vals),length(par2vals),NumIterations);
PropIntegratedMtx=zeros(length(par1vals),length(par2vals),NumIterations);
TSegMtx=zeros(length(par1vals),length(par2vals),NumIterations);
TIntMtx=zeros(length(par1vals),length(par2vals),NumIterations);
SepMtx=zeros(length(par1vals),length(par2vals),NumIterations);

SwitchData=struct([]);

svdir=[pwd,'/'];

if simpleFilenameFlag
    svname='buildup_data.mat';
else
    svname=['buildup_',bfflag,num2str(N),'x',num2str(M),'x',num2str(NumIterations),'_',...
        par1name,num2str(par1vals(1)),'_',num2str(par1vals(end)),'_',...
        par2name,num2str(par2vals(1)),'_',num2str(par2vals(end)),'_',...
        'g',num2str(g),'_',...
        'betae',num2str(betae),'_',...
        'betai',num2str(betai),'_',...
        'icf',num2str(ICFac),'_',...
        'mf',num2str(MichFac),'_',...
        'df',num2str(DFFac)%,'_',...
        ];svname=strrep(svname,'.','p');
    svname=strrep(svname,'-','m');
end

for i=1:length(par1vals)
    eval([par1name,'=',num2str(par1vals(i)),';']);
    
    for j=1:length(par2vals)
        eval([par2name,'=',num2str(par2vals(j)),';']);
        disp(['par1 value ',num2str(i),' of ',...
        num2str(length(par1vals)),'  ',par1name,'=',num2str(par1vals(i))])
        disp(['par2 value ',num2str(j),' of ',...
            num2str(length(par2vals)),'  ',par2name,'=',num2str(par2vals(j))])
        IpFunA=@(PR)(IpSlope*PR+IpConst);
        IpFunB=@(PR)(PpFrac*IpSlope)*PR+IpConst;
        
        IpampA=IpFunA(PR);
        IpampB=IpFunB(PR);

        per_b=2*per;
        sqwave=@(t,per1)(heaviside(t)-heaviside(t-per1));
        sqwaveper=@(t,per1,per2)sqwave(t-floor(t/per2)*per2,per1);
        DILong=@(t)...
            heaviside(t).*exp(2)./(dr).^2.*t.^2 .* exp(-2/(dr).*abs(t))+...
            1/6*heaviside(t).*exp(2)./(dr*5.5).^2.*t.^2 .* exp(-2/(dr*5.5).*abs(t));
       
        phA1=@(t)(t-floor(t/per)*per).*sqwaveper(t,per/2,per*2);
        phA2=@(t)(t-floor(t/per)*per).*sqwaveper(t-per,per,per*2);
        phB=@(t)(t-floor((t-per/2)/per_b)*per_b-per/2).*sqwaveper(t-per/2,per/2,2*per);% phase for B half period
        DIA=@(t)DILong(phA1(t))+DILong(phA2(t));
        DIB=@(t)DILong(phB(t));
                
        % Input dependence on DF
        tpFun=@(PR)tpConstExp*exp(-PR/tpSlopeExp);
        tphi=tpFun(PR);
        InputDF=@(phi) exp(-abs(phi)/tphi);
        
        SargInLoad=@(t,u)[betae*u(7)*u(10)-betai*u(1)-betai*(u(2)+u(3))-g*u(4);%...Excitation NMDA rAB
            betae*u(8)*u(11)-betai*u(2)-2*betai*u(1)-betai*u(3)-g*u(5);%...Excitation NMDA rA
            betae*u(9)*u(12)-betai*u(3)-2*betai*u(1)-betai*u(2)-g*u(6)];%...Excitation NMDA rB
      
        % DF increases from DF*(1-DFFac) to DF on timescale phidec
        DFFcn=@(t)DF-DFFac*DF*exp(-(1/phidec)*t);
        % Input amplitude decays from Ipamp*(1+MichFac) to Ipamp on
        % timescale phidec
        % if DFFac is zero, also set MichFac to zero; i.e. no amplitude adaptation
        IpampinA=@(t)IpampA+(DFFac~=0)*MichFac*IpampA*exp(-(1/phidec)*t);
        IpampinB=@(t)IpampB+(DFFac~=0)*MichFac*IpampB*exp(-(1/phidec)*t);
        for k=1:NumIterations 
            u0=zeros(12,1);
            u0(1:3)=.1*icRnd(:,k);
            u0(1)=ICFac*u0(1);
            u0(7:9)=u0(1:3);
            u0(10:12)=1;
            tint=0:tstep:tfin;
            t=tint;
            rsctime=100;trsc=t/rsctime;
            InVals=zeros(length(tint),3);
            InVals(:,1)=Ibase+InputDF(DFFcn(t)).*(IpampinA(t).*DIA(t)+IpampinB(t).*DIB(t));
            InVals(:,2)=Ibase+(IpampinA(t).*DIA(t))+InputDF(2*DFFcn(t)).*(IpampinB(t).*DIB(t));
            InVals(:,3)=Ibase+(IpampinB(t).*DIB(t))+InputDF(2*DFFcn(t)).*(IpampinA(t).*DIA(t));
            rng(rngSeeds(k));
            tic
            
            [y,yn]=SDERun(tint,f,SargInLoad,u0,taua,taunm,taux,taud,synd,abs(gamma),InVals,[]);
          
            MinWindow=0.25;ifplot=0;        
            SlidingWidthT=50;SlidingWidthIdx=SlidingWidthT/tstep;
            uFilt=slidingavg(y(:,1),SlidingWidthIdx);
            vFilt=slidingavg((y(:,2)+y(:,3))/2,SlidingWidthIdx);
            
            [NumberOfSwitches,First,Durations,DurationsInt,DurationsSeg,PropIntegrated,...
                SwitchTimes,SwDirections]=...
                ProcSwTimes(trsc,MinWindow,uFilt,vFilt,ifplot);
            SepMtx(i,j,k)=mean(abs(uFilt-vFilt));
            if k==1
                disp(['Iteration ',num2str(k),' of ',...
                    num2str(NumIterations),'; Number of switches: ',num2str(NumberOfSwitches)])
            end
            PropIntegratedMtx(i,j,k)=PropIntegrated;
            tintCoarse=0:tstep*10:tfin;
            uFiltCoarse=interp1(tint,uFilt,tintCoarse);
            vFiltCoarse=interp1(tint,vFilt,tintCoarse);
            SwitchData(i,j,k).uFiltCoarse=uFiltCoarse;
            SwitchData(i,j,k).vFiltCoarse=vFiltCoarse;
            SwitchData(i,j,k).First=First;
            SwitchData(i,j,k).Durations=Durations;
            SwitchData(i,j,k).DurationsInt=DurationsInt;
            SwitchData(i,j,k).DurationsSeg=DurationsSeg;
            SwitchData(i,j,k).SwitchTimes=SwitchTimes;
            SwitchData(i,j,k).SwDirections=SwDirections;
            if NumberOfSwitches==0
                NumberOfSwitchesMtx(i,j,k)=sign(uFilt(end)-vFilt(end));
            else
                NumberOfSwitchesMtx(i,j,k)=NumberOfSwitches;
            end
        end
    end
end

if strcmp(bfflag,'single')
    figure(1);clf;hold on
    set(gcf,'units','centimeters','position',[0,0,10,4])
    set(gca,'position',[0.1200    0.1200    0.83    0.79])
    addidx=6;
    adapaddidx=3;
    make_colors;
    clrs={blue,red,green};
    plot(trsc,y(:,1+addidx),'color',clrs{1})
    plot(trsc,y(:,2+addidx),'color',clrs{2})
    plot(trsc,y(:,3+addidx),'color',clrs{3})
    plot(trsc,y(:,1+adapaddidx),'--','color',clrs{1})
    plot(trsc,y(:,2+adapaddidx),'--','color',clrs{2})
    plot(trsc,y(:,3+adapaddidx),'--','color',clrs{3})
    lh=legend('AB','A','B','location','west');
    set(lh,'box','off','fontname','helvetica','fontsize',8)
    tf=0.2;
    plot(trsc,tf*ones(size(trsc)),'k--');
    for i=1:length(SwitchTimes)
        plot([SwitchTimes(i),SwitchTimes(i)],[-1,2],'k-');
    end
    set(gca,'ylim',[0,1],'xlim',[0,tfin/rsctime])
    set(gca,'ytick',[0:0.5:1])
    set(gca,'fontname','helvetica','fontsize',8)
    drawnow
    plotdir=pwd;
    if ifexport
        export_fig([plotdir,'model-time-history-',num2str(addidx),'.pdf']);
    end
    
    return
end

if simpleFilenameFlag
    savename=[svdir,svname];
else
    datestring=datestr(floor(now));
    datestring=strrep(datestring,'-','');
    timestring=datestr(rem(now,1),13);
    timestring=strrep(timestring,':','');
    clear y
    savename=[svdir,svname,'_',datestring,'_',timestring,'.mat'];
end
save(savename);
disp(savename)

