clear all
close all

data = readtable('data-all-switches.csv');
% Loads table with columns: trialnum,subject,percept,TMod,dur,ontime,offtime,direction,phase,SwTGroup

TMods = [0,5,10,20];
NumSubj = 17;
make_colors

%% Single trial plots A-D
substr='s17';
trial=[1 7 15 28];
DF0=5;
AmpMod=1.5;
tend=90;

tripletlength = 0.5; 
tripTime=0:tripletlength:(tend+10)/tripletlength;

for i = 1:length(trial)
    trialind = ['t',num2str(trial(i))];
    
    indI = find(strcmp(trialind,data{:,1})&strcmp(substr,data{:,2})&strcmp('int',data{:,3}));
    indS = find(strcmp(trialind,data{:,1})&strcmp(substr,data{:,2})&strcmp('seg',data{:,3}));
    
    OnTimesI = data{indI,6};
    OffTimesI = data{indI,7};
    OnTimesS = data{indS,6};
    OffTimesS = data{indS,7};
    
    figure
    box on
    set(gca,'ytick',[0.25,0.75,1.1,1.6],'yticklabel',{'Seg','Int','DF=3.5','DF=6.5'});
    set(gca,'TickLabelInterpreter','latex')
    
    hold on
    for j=1:length(OnTimesI)
        XY=[OnTimesI(j),.5;
            OnTimesI(j),1;
            OffTimesI(j),1;
            OffTimesI(j),.5;];
        fill(XY(:,1),XY(:,2),blue,'linewidth',1);
    end
    for k=1:length(OnTimesS)
        XY=[OnTimesS(k),0;
            OnTimesS(k),.5;
            OffTimesS(k),.5;
            OffTimesS(k),0;];
        fill(XY(:,1),XY(:,2),red,'linewidth',1);
    end
    
    if TMods(i)~=0
        fMaxTimes=(TMods(i)/4:TMods(i):tend);
        fMinTimes=(3*TMods(i)/4:TMods(i):tend);
        
        for m=1:length(fMaxTimes)
            plot([fMaxTimes(m),fMaxTimes(m)],[0,1.6],'-','color',grey,'linewidth',1);
        end
        
        for m=1:length(fMinTimes)
            plot([fMinTimes(m),fMinTimes(m)],[0,1.6],'--','color',grey,'linewidth',1);
        end
        than=title(['$T_\textrm{mod} = ',num2str(TMods(i)),'\,s$'],'fontweight','normal');
        DFArray=DF0+AmpMod*sin(2*pi*tripTime./TMods(i));
    else
        DFArray=DF0*ones(size(tripTime));
        than=title(['Modulation off'],'fontweight','normal');
    end
    
    dfmin=DF0-AmpMod;
    dfmax=DF0+AmpMod;
    
    plot(tripTime,(DFArray-3.5)/6+1.1,'color',grey,'linewidth',1.5)
    
    set(than,'interpreter','latex')
    set(gcf,'units','centimeters','position',[3,6,20,8]);
    box on
    set(gca,'xlim',[0 tend],'ylim',[-0.05 1.75]);
    set(gca,'fontname','Times','fontsize',20);
    set(gca,'linewidth',1.5);set(gcf,'color','w');box on
    set(gca,'xtick',0:40:tend)
    
end




%% Phase histograms E & F
nobins=20;
h=2*pi/nobins;
bins=h/2:h:2*pi-h/2;
vonmises = @(p,x)exp(p(1)*cos(x-p(2)))/(2*pi*besseli(0,p(1)));
theta =linspace(0,2*pi,100);

for i = 1:length(TMods)
    
    indSI = find(strcmp('S2I',data{:,8})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
    indIS = find(strcmp('I2S',data{:,8})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
    
    % Recentre data around forcing max/min
    PhaseIS = mod(data{indIS,9}+pi/2,2*pi);
    PhaseSI = mod(data{indSI,9}-pi/2,2*pi);
    
    [N1,X1]=histnorm(PhaseIS,bins);
    [N2,X2]=histnorm(PhaseSI,bins);
    vonmises = @(p,x)exp(p(1)*cos(x-p(2)))/(2*pi*besseli(0,p(1)));
    p0IS = [std(PhaseIS),mean(PhaseIS)];
    p0SI = [std(PhaseSI),mean(PhaseSI)];
    coefIS = nlinfit(X1,N1,vonmises,p0SI);
    coefSI = nlinfit(X2,N2,vonmises,p0SI);

    colorsS2I={[0.1 0.1 0.9],[0.2 0.5 0.8],[0.4 0.7 0.7],[0.1 0.1 0.5]};
    figure(5)
    hold on
    if TMods(i) == 20
        plot(linspace(0,2*pi),0.25+0.25*sin(linspace(0,2*pi)+pi/2),'k')
        plot([pi pi],[0,1],'k-','linewidth',1.5);
        plot([0 0],[0,1],'k--','linewidth',1.5);
        plot(theta,vonmises(coefSI,theta),'--','Color',colorsS2I{1},'linewidth',1.5)
        g1=plot(X2,N2,'Color',colorsS2I{1},'linewidth',1.5);
        legend([g4,g3,g2,g1],{'TMod0','TMod5','TMod10','TMod20'})
        set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
        set(gca,'XTickLabels',{'$\frac{\pi}{2}$','$\pi$', '$\frac{3\pi}{2}$', '0', '$\frac{\pi}{2}$'})
        set(gca,'TickLabelInterpreter','latex')
        set(gca,'linewidth',1.5,'fontsize',16,'fontname','Helvetica')
    elseif TMods(i) == 10
        plot(theta,vonmises(coefSI,theta),'--','Color',colorsS2I{2},'linewidth',1.5)
        g2=plot(X2,N2,'Color',colorsS2I{2},'linewidth',1.5);
    elseif TMods(i) == 5
        plot(theta,vonmises(coefSI,theta),'--','Color',colorsS2I{3},'linewidth',1.5)
        g3=plot(X2,N2,'Color',colorsS2I{3},'linewidth',1.5);
    else
        plot(theta,vonmises(coefSI,theta),'--','Color',colorsS2I{4},'linewidth',1.5)
        g4=plot(X2,N2,'Color',colorsS2I{4},'linewidth',1.5);
    end
    axis tight
    set(gca,'ylim',[0,0.51])
    
    colorsI2S={[1 0 0],[1 0.5 0],[0.7 0.1 0.3],[0.3 0 0.1]};
    figure(6)
    hold on
    if TMods(i) == 20
        plot(linspace(0,2*pi),0.25+0.25*sin(linspace(0,2*pi)-pi/2),'k')
        plot([0 0],[0,1],'k-','linewidth',1.5);
        plot([pi pi],[0,1],'k--','linewidth',1.5);
        plot(theta,vonmises(coefIS,theta),'--','Color',colorsI2S{1},'linewidth',1.5)
        x1=plot(X1,N1,'Color',colorsI2S{1},'linewidth',1.5);
        legend([x4,x3,x2,x1],{'TMod0','TMod5','TMod10','TMod20'})
        set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
        set(gca,'XTickLabels',{'$\frac{3\pi}{2}$', '0', '$\frac{\pi}{2}$', '$\pi$','$\frac{3\pi}{2}$'})
        set(gca,'TickLabelInterpreter','latex')
        set(gca,'linewidth',1.5,'fontsize',16,'fontname','Helvetica')
    elseif TMods(i) == 10
        plot(theta,vonmises(coefIS,theta),'--','Color',colorsI2S{2},'linewidth',1.5)
        x2=plot(X1,N1,'Color',colorsI2S{2},'linewidth',1.5);
    elseif TMods(i) == 5
        plot(theta,vonmises(coefIS,theta),'--','Color',colorsI2S{3},'linewidth',1.5)
        x3=plot(X1,N1,'Color',colorsI2S{3},'linewidth',1.5);
    else
        plot(theta,vonmises(coefIS,theta),'--','Color',colorsI2S{4},'linewidth',1.5)
        x4=plot(X1,N1,'Color',colorsI2S{4},'linewidth',1.5);
    end
    axis tight
    set(gca,'ylim',[0,0.51])
end


%% Boxplots H

kappaIS = zeros(length(TMods),NumSubj);
thetaIS = zeros(length(TMods),NumSubj);
kappaSI = zeros(length(TMods),NumSubj);
thetaSI = zeros(length(TMods),NumSubj);

for j = 1:NumSubj
    for i = 1:length(TMods)
        SubjectString=['s',num2str(j)];
        
        ind = find(strcmp(SubjectString,data{:,2})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
        
        kappaIS(i,j) = data{ind(1),12};
        kappaSI(i,j) = data{ind(1),13};
        
        thetaIS(i,j) = data{ind(1),14};
        thetaSI(i,j) = data{ind(1),15};
        
    end
end

make_colors

figure(7)
boxplot(kappaIS',TMods,'orientation','horizontal','colors',red)
set(findall(gcf,'-property','Fontname'),'Fontname','Times')
xlabel('Shape $\kappa$ I to S','interpreter','latex')
figure(7)
set(gca,'yticklabel',{'$T_\textrm{m}$ off','$T_\textrm{m}\!\!=\!\!5$'...
    ,'$T_\textrm{m}\!\!=\!\!10$','$T_\textrm{m}\!\!=\!\!20$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'xlim',[0,8])
set(findall(gcf,'-property','FontSize'),'FontSize',10)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)

% S2I only
figure(8)
boxplot(kappaSI',TMods,'orientation','horizontal','colors',blue)
set(findall(gcf,'-property','Fontname'),'Fontname','Times')
xlabel('Shape $\kappa$ S to I','interpreter','latex')
figure(8)
set(gca,'yticklabel',{'$T_\textrm{m}$ off','$T_\textrm{m}\!\!=\!\!5$'...
    ,'$T_\textrm{m}\!\!=\!\!10$','$T_\textrm{m}\!\!=\!\!20$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'xlim',[0,8])
set(findall(gcf,'-property','FontSize'),'FontSize',8)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
