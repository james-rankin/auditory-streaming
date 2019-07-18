function [NumberOfSwitches,First,Durations,Durationsu,Durationsv,PropIntegrated,...
    SwitchTimes,SwDirections]=...
    ProcSwTimesT0(trsc,t0,MinWindow,uFilt,vFilt,ifplot)

t0idx=find(trsc>t0,1);
First=[0,sign(uFilt(t0idx)-vFilt(t0idx))];
SwDirections=diff(sign(uFilt-vFilt));
SwIdx=find(diff(sign(uFilt-vFilt))~=0);
SwIdx(trsc(SwIdx)<trsc(t0idx))=[];
RemoveIdx=find(diff(trsc(SwIdx))<MinWindow);
ToBeRemovedIdx=[];
for i=1:length(RemoveIdx)
    ToBeRemovedIdx=[ToBeRemovedIdx;RemoveIdx(i);RemoveIdx(i)+1];
end
SwIdx(ToBeRemovedIdx)=[];
SwitchTimes=trsc(SwIdx);
SwDirections=sign(SwDirections(SwIdx));
NumberOfSwitches=numel(SwIdx);
Durations=[];
Durationsv=[];
Durationsu=[];
PropIntegrated=numel(find(uFilt>vFilt))/numel(trsc);
if numel(SwIdx)>=1
    firstIdx=t0idx;
    First=[trsc(SwIdx(1)),sign(uFilt(firstIdx)-vFilt(firstIdx))];
end
if numel(SwIdx)>=2
    Durations=diff([0,trsc(SwIdx)]);
    Durations=Durations(:);
    for i=1:length(SwIdx)
        if ifplot
        subplot(3,1,1);
        plot([trsc(SwIdx(i)),trsc(SwIdx(i))],[-1 2],'k-','linewidth',2);
        end
        if i<length(SwIdx)
        if SwDirections((i))>0
            if ifplot
            text(trsc(SwIdx(i)),0.5,'In')
            end
            Durationsu=[Durationsu;trsc(SwIdx(i+1))-trsc(SwIdx(i))];
        else
            if ifplot
            text(trsc(SwIdx(i)),0.5,'Sg')
            end
            Durationsv=[Durationsv;trsc(SwIdx(i+1))-trsc(SwIdx(i))];
        end
        end
        
    end
    
end