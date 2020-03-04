function opts = SetAnalysisOptions(varargin)
%CNMF Options
opts.K = 50;              %10                                  %Number of factors
opts.L = 13;              %100                                 %Length (timebins) of each factor exemplar
opts.maxiter =300;        %100                                 %Maximum # iterations to run for seqNMf
opts.lambda = 0.0005;     %0.0005                              %Regularization parameter
opts.lambdaL1H = 1;       %0                                   %L1 sparsity parameter; Increase to make H's more sparse
opts.lambdaOrthoH = 1;    %0                                   %||HSH^T||_1,i~=j; Encourages events-based factorizations
opts.tolerance = 0;       %0                                   %Stop fitting if error reaches said value
opts.shift = 0;    

%Additional Options
opts.showPlot = 0; 
opts.movie_flag = 0; 
opts.block = 1; 
opts.RemoveEmptyWs = 1; 

%Directory Navigation Options
opts.base = 'Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018';
opts.bucket = '/jukebox/buschman/'; 
opts.data_file_name = '/AnalyzedData_MesomappingManuscript_5_2019/DFF_Data/AllData_binned_SmallMask_3x_2minTraining.mat';



end