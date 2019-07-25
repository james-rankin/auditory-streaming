clear all
close all

data = readtable('data-all-switches.csv');
% Loads table with columns: trialnum,subject,percept,TMod,dur,ontime,offtime,direction,phase,SwTGroup

% Experimental parameters
NumSubj = 17;                                   % Known number of subjects
TrialLength = 180;                              % Known trial length
TMods = [0 5 10 20];                            % Modulation periods
trials = {1:6,7:12,13:21,22:33};                % Trial indices for each modulation period

colors = {'r','b','g'};

%% Duration histograms A 
bins=0:0.5:TrialLength;

figure(1);clf;hold on

set(gcf,'units','centimeters','position',[6,12,15,14]);

for i = 1:length(TMods)
    
    indF = find(strcmp(['TMod',num2str(TMods(i))],data{:,4})&strcmp('fast',data{:,10}));
    indM = find(strcmp(['TMod',num2str(TMods(i))],data{:,4})&strcmp('medium',data{:,10}));
    indS = find(strcmp(['TMod',num2str(TMods(i))],data{:,4})&strcmp('slow',data{:,10}));
    
    DurF = data{indF,5};
    DurM = data{indM,5};
    DurS = data{indS,5};
    
    f1 = figure(i);
    set(f1,'units','centimeters','position',[3,6,12.5,12]);
    hold on
    % Duration histograms
    [NF,XF]=histnorm(DurF,bins);
    plot(XF,NF,'Color',colors{1},'LineWidth',2);
    [NM,XM]=histnorm(DurM,bins);
    plot(XM,NM,'Color',colors{2},'LineWidth',2);
    [NS,XS]=histnorm(DurS,bins);
    plot(XS,NS,'Color',colors{3},'LineWidth',2);
    
    axis([0 20 0 0.45])
    axis square
    xlabel('Duration (s)')
    ylabel('Probability')
    
    set(gca,'linewidth',1.5)
    set(findall(gcf,'-property','Fontname'),'Fontname','Times')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    legend({'Fast','Medium','Slow'},'fontsize',16,'box','off','location','northeast')

end


%% Duration boxplot B
MeanDur = zeros(length(TMods),NumSubj);
MeanDurI = zeros(length(TMods),NumSubj);
MeanDurS = zeros(length(TMods),NumSubj);
SwType = {};

for i = 1:length(TMods)
    for j = 1:NumSubj
        substr=['s',num2str(j)];
        
        ind = find(strcmp(substr,data{:,2})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
        
        MeanDur(i,j) = mean(data{[ind],5});
    
        SwType{j} = data{ind(1),10}{1};
    end
end

MeanDurF = MeanDur(:,strcmp('fast',SwType));
MeanDurM = MeanDur(:,strcmp('medium',SwType));
MeanDurS = MeanDur(:,strcmp('slow',SwType));

figure(5)
hold on
boxplot(MeanDurF',TMods,'orientation','horizontal','colors',colors{1})
boxplot(MeanDurM',TMods,'orientation','horizontal','colors',colors{2})
boxplot(MeanDurS',TMods,'orientation','horizontal','colors',colors{3})
figure(5)
axis square
axis([0 30 0.5 4.5])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
xlabel('Duration (s)','Interpreter', 'Latex','FontSize',24)
ylabel('$T_{\rm mod}$','Interpreter', 'Latex','FontSize',24)
set(gca,'fontsize',20,'fontname','Times')



%% Switches per cycle C
overlap=0.2;
n = 5;                      % number of TMod cycles for binning
type = {'fast','medium','slow'};
NumPerType = [sum(strcmp('fast',SwType)),sum(strcmp('medium',SwType)),sum(strcmp('slow',SwType))];
SubRange = 1:17;

for i = 2:length(TMods)
   
    binwidth = n*TMods(i);
    binS = 0:overlap*binwidth:TrialLength-binwidth;
    A = [];
    
    figure
    hold on
    set(gcf,'units','centimeters','position',[3,6,12.5,12]);
    
    for m = 1:length(type)
        
        if strcmp(type(m),'fast')
            jrange = SubRange(strcmp('fast',SwType));
        elseif strcmp(type(m),'medium')
            jrange = SubRange(strcmp('medium',SwType));
        elseif strcmp(type(m),'slow')
            jrange = SubRange(strcmp('slow',SwType));
        end
        
        for j = jrange
            
            for k = trials{i} %size(SwTimesCell,3)
                
               ind = find(strcmp(['t',num2str(k)],data{:,1})&strcmp(['s',num2str(j)],data{:,2}));
            
            SwTimes = data{ind,7};

                for p = 1:length(binS)
                    
                    A = [A,(sum(SwTimes>binS(p)&SwTimes<binS(p)+binwidth)/(2*n))];
                    
                end
                
            end
            
        end
        
        A = A(A>0);
        [X,Y] = histnorm(reshape(A,size(A,1)*size(A,2),1),0:0.1:5);
        plot(Y,X,colors{m},'LineWidth',2)
        
        
    end
    
    xlabel('Switches per cycle')
    ylabel('Probability density')
    axis square
    axis([0 5 0 4])
    set(gca,'linewidth',1.5)
    set(findall(gcf,'-property','Fontname'),'Fontname','Times')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    lh=legend({'Fast','Medium','Slow'},'fontsize',16,'box','off','location','east');
   
end