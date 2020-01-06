load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\AnalyzedData_MesomappingManuscript_5_2019\TrainRepitoires\TrainingFit_CompareDiscoveryMethods\block_32.mat')

%remove any empty rows from nmf
empty_elems = arrayfun(@(s) all(structfun(@isempty,s)), nmf);
nmf(empty_elems==1)=[];

%Get the values used
xvals = [1:1:25,50:25:300];

figure; hold on; 
%Plot the nmf results
plot(xvals,[nmf.ExpVar_all]*100,'linewidth',2,'color','r');

%Plot the PCA results
plot(1:300,spca.ExpVar_all(1:300)*100,'linewidth',2,'color','b');

%Plot the line of the 
plot(1:300,repmat(seq.ExpVar_all,1,300)*100,'k','linewidth',2,'linestyle','--')

xlim([0 150])

xlabel('Number of Components/Factors')
ylabel('PEV')

legend('nmf','pca','cnmf (22motifs)','Location','southeast')
title('Comparing Motif CNMF to 1D nmf and pca')