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

type = 'ALL';

if strcmp(type,'FAST')
    jrange = [1,3,4,5,11,14];
elseif strcmp(type,'MEDIUM')
    jrange = [2,8,10,12,13,16,17];
elseif strcmp(type,'SLOW')
    jrange = [6,7,9,15];
else
    jrange = 1:NumSubj;
end

for i = 2:length(TMods)
    
    temp = linspace(0,2*pi,Tcycle(i)+1);
    Phase = [repmat(temp(1:end-1),1,max(t)/TMods(i)-1),temp];
    
    DataForMapI = zeros(Tcycle(i)*6+1,length(bins));
    DataForMapS = zeros(Tcycle(i)*6+1,length(bins));
    
    for j = jrange
        
        for k = trials{i}
            
            indI = find(strcmp(['t',num2str(k)],data{:,1})&strcmp(['s',num2str(j)],data{:,2})&strcmp('int',data{:,3})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
            indS = find(strcmp(['t',num2str(k)],data{:,1})&strcmp(['s',num2str(j)],data{:,2})&strcmp('seg',data{:,3})&strcmp(['TMod',num2str(TMods(i))],data{:,4}));
            
            A = [data{indI,6},data{indI,7}];
            Int = zeros(length(t),1);
            for nI = 1:size(A,1)
                Int = Int + (t>=A(nI,1)&t<A(nI,2))';
            end
            
            B = [data{indS,6},data{indS,7}];
            Seg = zeros(length(t),1);
            for nS = 1:size(B,1)
                Seg = Seg + (t>=B(nS,1)&t<B(nS,2))';
            end
            
            S2IIndex = x(Int(2:end)>Int(1:end-1));
            S2I = Phase(Int(2:end)>Int(1:end-1));
            I2SIndex = x(Seg(2:end)>Seg(1:end-1));
            I2S = Phase(Seg(2:end)>Seg(1:end-1));
            
            for n = 1:length(bins)-1
                
                S2Itemp = S2IIndex.*(S2I>=bins(n)&S2I<bins(n+1));
                I2Stemp = I2SIndex.*(I2S>=bins(n)&I2S<bins(n+1));
                S2Itemp = S2Itemp(S2Itemp~=0 & S2Itemp>round(Tcycle(i)*(3)) & S2Itemp<length(Int)-round(Tcycle(i)*(3)));
                I2Stemp = I2Stemp(I2Stemp~=0 & I2Stemp>round(Tcycle(i)*(3)) & I2Stemp<length(Seg)-round(Tcycle(i)*(3)));
                
                for mI = 1:length(S2Itemp)
                    % Normalized by number of TMod cycles
                    DataForMapI(:,n) = DataForMapI(:,n) + Int(S2Itemp(mI)-round(Tcycle(i)*(3)):S2Itemp(mI)+round(Tcycle(i)*(3)))/(180/TMods(i)*length(trials{i})*length(jrange));
                end
                for mS = 1:length(I2Stemp)
                    % Normalized by number of TMod cycles
                    DataForMapS(:,n) = DataForMapS(:,n) + Seg(I2Stemp(mS)-round(Tcycle(i)*(3)):I2Stemp(mS)+round(Tcycle(i)*(3)))/(180/TMods(i)*length(trials{i})*length(jrange));
                end
                
            end
            
        end
        
    end
    
    cax = 0.0362;%max(max([DataForMapI,DataForMapS]))
    
    figure
    surf(bins,linspace(-3,3,size(DataForMapI,1)),DataForMapI)
    view(2)
    shading interp
    axis tight
    set(gcf, 'Position', [100 500 520 370])
    box on
    colormap jet
    colorbar
    caxis([0 cax]);
    set(gca,'linewidth',1.5,'fontsize',18,'fontname','Times')
    set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
    set(gca,'XTickLabels',{'0', '$\pi/2$','$\pi$', '$3\pi/2$', '$2\pi$'})
    set(gca,'TickLabelInterpreter','latex')
    
    figure
    surf(bins,linspace(-3,3,size(DataForMapS,1)),DataForMapS)
    view(2)
    shading interp
    axis tight
    set(gcf, 'Position', [100 500 520 370])
    colormap jet
    colorbar
    caxis([0 cax]);
    set(gca,'linewidth',1.5,'fontsize',18,'fontname','Times')
    set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
    set(gca,'XTickLabels',{'0', '$\pi/2$','$\pi$', '$3\pi/2$', '$2\pi$' })
    set(gca,'TickLabelInterpreter','latex')
    
end
