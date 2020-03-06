% Supporting information for manuscript ``Neuromechanistic Model of
% Auditory Bistability'' by James Rankin, Elyse Sussman and John Rinzel
% Contact: james.rankin@nyu.edu

% The data file IndividualDurationsByTrial.mat contains the raw durations
% from the experiments reported in the manuscript

% This file illustrates how to extract individual durations from the data
% structure

% Note subject 1 was excluded from most of the analysis, see the MS for
% details


data=load('IndividualDurationsByTrial.mat');
% Data are index by 
% - the condition number i (=1..8)
% - the subject number j (=1..16)
% - the repetition number k (=1..3)

% The frequency difference values for the 8 conditions are in 
DFvals=data.DFvals;

% For a given trial i,j,k we can extract the integrated (One stream) and
% segregated (Two streams) durations from cell arrays as follows
i=4;j=2;k=1; % DF=5, subject 2, repetition 1
disp(['Frequency difference: ',num2str(DFvals(i)),' st']);
disp(['Subject number: ',num2str(j)]);
disp(['Repetition number: ',num2str(k)]);
DurationsInt=data.DurationsStructOne{i,j,k}
DurationsSeg=data.DurationsStructTwo{i,j,k}
% Empty cells mean no durations after the first were registed on that trial


% We can extract the first duration as follows
DurationFst=data.DurationsStructFst{i,j,k}
% The first element in the 1x2 vector is the duration and the second
% percept type 1=Integrated/One stream 2=Segregated/Two streams