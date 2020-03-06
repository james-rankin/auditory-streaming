function [y,varargout] = SDERun(tint,f,Sarg,u0,taua,taunm,taux,taud,synd,gamma,InVals,yn_in)
%% 
% Code to accompany the paper Byrne, Rinzel and Rankin (2019)
% Entrainment of stream segregation in a dynamic environment
% Contact: james.rankin@gmail.com
% If you use or adapt this code acknowledge us by citing our paper
%
% Euler-Muryama scheme to timestep model
dim=length(u0);
dimn=3;


y=zeros(length(tint),dim);
yn=zeros(length(tint),dimn);
dt=tint(2)-tint(1);

% disp(['dt=',num2str(dt)]);
N=length(tint);

y(1,:)=u0(:);
yn(1,:)=zeros(dimn,1);
doloadnoise=~isempty(yn_in);

for i=2:N
    ti=tint(i-1);
    ui=y(i-1,:);
   	uin=yn(i-1,:);
    Fn=-1/taux*uin;
    SargVals=Sarg(ti,ui);
    Fu=[-ui(1)+f(SargVals(1)+InVals(i-1,1)+uin(1)),...
        -ui(2)+f(SargVals(2)+InVals(i-1,2)+uin(2)),...
        -ui(3)+f(SargVals(3)+InVals(i-1,3)+uin(3)),...
        (-ui(4)+ui(1))/taua,...
        (-ui(5)+ui(2))/taua,...
        (-ui(6)+ui(3))/taua,...
        (-ui(7)+ui(1))/taunm,...
        (-ui(8)+ui(2))/taunm,...
        (-ui(9)+ui(3))/taunm,...
        (-ui(10)+(1-synd*ui(1)))/taud,...
        (-ui(11)+(1-synd*ui(2)))/taud,...   
        (-ui(12)+(1-synd*ui(3)))/taud];
    y(i,:)= ui + dt*(Fu);
    if doloadnoise  
        yn(i,:) = yn_in(i,:);
    else
     yn(i,:) = uin + dt *Fn + gamma*sqrt(dt)*sqrt(2/taux)*randn(1,dimn);
    end
end

if nargout==2;
    varargout{1}=yn;
end
