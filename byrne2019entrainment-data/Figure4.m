clear all
close all

data = readtable('data-all-switches.csv');
% Loads table with columns: trialnum,subject,percept,TMod,dur,ontime,offtime,direction,phase,SwTGroup

% Experimental parameters
NumSubj = 17;                                   % Known number of subjects
TrialLength = 180;                              % Known trial length
TMods = [0 5 10 20];                            % Modulation periods
trials = {1:6,7:12,13:21,22:33};                % Trial indices for each modulation period

% Temporal resolution (ms) to analyse switch times at
dt = 0.1;
t = 0:dt:TrialLength;

Tcycle = TMods/dt;
bins = linspace(0,2*pi,51);
PhaseDur = zeros(length(bins),100);
x = 1:length(t)-1;


DataForMap = zeros(length(bins),length(x));


%% Duration histograms A 

figure(1);clf;hold on

set(gcf,'units','centimeters','position',[6,12,15,14]);
legstrs={'Modulation off','TMod = 5s','TMod = 10s','TMod = 20s'};

colors = {[0.2 0.7 1],[1 0.1 0.15],[1 0.55 0.1],[0.4 0.85 0.2]};
for i = 1:length(TMods)
    
    ind = find(strcmp(['TMod',num2str(TMods(i))],data{:,4}));
    
    Dur = data{ind,5};
    
    % Duration histograms
    bins=0:0.5:180;
    [N,X]=histnorm(Dur,bins);
    plot(X,N,'Color',colors{i},'LineWidth',2);
    xlabel('Duration (s)')
    ylabel('Probability')
    axis square
    axis([0 15 0 0.4])

    

end
%%
set(findall(gcf,'-property','Fontname'),'Fontname','Times')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gca,'linewidth',1.5)
legend(legstrs,'box','off','FontSize',16);


%% Duration boxplot B
MeanDur = zeros(length(TMods),NumSubj);
MeanDurI = zeros(length(TMods),NumSubj);
MeanDurS = zeros(length(TMods),NumSubj);
SwType = {};

for i = 1:length(TMods)
    for j = 1:NumSubj
        substr=['s',num2str(j)];
        
        indI = find(strcmp(substr,data{:,2})&strcmp('int',data{:,3})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
        indS = find(strcmp(substr,data{:,2})&strcmp('seg',data{:,3})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
        
        MeanDur(i,j) = mean(data{[indI;indS],5});
        MeanDurI(i,j) = mean(data{indI,5});
        MeanDurS(i,j) = mean(data{indS,5});
    
        SwType{j} = data{indI(1),10}{1};
    end
end

figure(2)
boxplot(MeanDur',TMods,'orientation','horizontal')
figure(2)
axis square
axis([0 30 0.5 4.5])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
xlabel('Duration (s)','Interpreter', 'Latex','FontSize',24)
ylabel('$T_{\rm mod}$','Interpreter', 'Latex','FontSize',24)
set(gca,'fontsize',20,'fontname','Times')



%% Duration scatter C

figure(3)
clf
hold on
set(gcf,'units','centimeters','position',[3,6,14.5,14]);
make_colors
axis equal
axis square

FastInd = strcmp(SwType,'fast');
MediumInd = strcmp(SwType,'medium');
SlowInd = strcmp(SwType,'slow');

plot(MeanDurI(1,FastInd),MeanDurS(1,FastInd),'o','color',black,'markerfacecolor',red,'markersize',10); 
plot(MeanDurI(1,MediumInd),MeanDurS(1,MediumInd),'o','color',black,'markerfacecolor',blue,'markersize',10);
plot(MeanDurI(1,SlowInd),MeanDurS(1,SlowInd),'o','color',black,'markerfacecolor',green,'markersize',10);

xlabel('Mean duration I')
ylabel('Mean duration S')
axis([0 20 0 20])
set(gca,'ytick',0:5:20,'xtick',0:5:20)
set(findall(gcf,'-property','Fontname'),'Fontname','Times')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)

x=linspace(0,20);
y=10-x;
plot(x,y,'k--')
y=20-x;
plot(x,y,'k--')
lh=legend({'Fast','Medium','Slow'},'fontsize',16,'box','off','location','southeast');



