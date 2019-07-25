# auditory-streaming
Neuromechanistic modelling of auditory streaming

Code to accompany the paper Byrne, Rinzel and Rankin (2019)
Auditory streaming and bistability paradigm extended to a dynamic environment
Contact: james.rankin@gmail.com

If you use or adapt this code acknowledge us by citing our paper (presently under review)

To produce Fig 1B:\
RunModelPF\
or\
RunModelPF('single')

To produce quick  (less 'trials') version of Fig 1C--D:\
RunModelPF('pfp') \
or the version as in the paper\
RunModelPF('pfp',320)

If you don't have a Parallel Computing Toolbox license replace two instances of 'parfor'
with 'for' in RunModelPF.m
