function stats = BehavioralStateClassification_Nway(avg_h_classification)

stats = [];
for cur_round = 3
rng('default')
%balance by grabbing random samples from each state
num_samples = cellfun(@(x) size(x,1),avg_h_classification,'UniformOutput',0);
num_samples = min([num_samples{:}]);

n_test = floor(cur_round/10*num_samples);
n_train = num_samples-n_test-13;

classification_data = cell(1,size(avg_h_classification,2));
for i = 1:size(avg_h_classification,2)
   temp = avg_h_classification{i};
   %get linearly spaced subsamples from throughout the recording
   temp = temp(floor(linspace(1,size(temp,1),num_samples)),:);   
   %leave a 1 sec frame gap between final two samples
   temp = cat(1,temp(1:n_train,:),temp((size(temp,1)-n_test)+1:end,:));
   classification_data{i} = temp;
end
    
labels = ones(n_test+n_train,1).*(1:size(classification_data,2));
labels = labels(:);
% classification_data = cat(1,classification_data{:});

n_class =numel(unique(labels(:)));

targetsVector = cell(1,size(comparisons,1));
outputsVector = cell(1,size(comparisons,1));
auc = zeros(numel(unique(labels(:))));
cm = zeros(numel(unique(labels(:))));
for i = 1:n_class
   fprintf('\nworking on comparision %d',i)
   idx = (1:n_class)~=i;
   
   fact = n_class-1;
   class1_train = classification_data{i}(1:(floor(n_train/fact)*fact),:);
   class1_test = classification_data{i}(end-(floor(n_test/fact)*fact)+1:end,:);
   class2_train = cellfun(@(x) x(randperm(n_train,floor(n_train/fact)),:),classification_data(idx),'UniformOutput',0);
   class2_test = cellfun(@(x) x(randperm(n_test,floor(n_test/fact))+(n_train),:),classification_data(idx),'UniformOutput',0);
      
   class2_train = cat(1,class2_train{:});
   class2_test = cat(1,class2_test{:});
   
   temp_data = cat(1,class1_train,class1_test,class2_train,class2_test); 
   temp_resp = cat(1,ones(size(class1_train,1),1), ones(size(class1_test,1),1), ones(size(class2_train,1),1)*2, ones(size(class2_test,1),1)*2);
   
   %split into training and test seperateing by 1 sec;
   cvp.training = logical(cat(1,ones(size(class1_train,1),1), zeros(size(class1_test,1),1), ones(size(class2_train,1),1), zeros(size(class2_test,1),1))); 
   cvp.test = logical(cat(1,zeros(size(class1_train,1),1), ones(size(class1_test,1),1), zeros(size(class2_train,1),1), ones(size(class2_test,1),1))); 

   %set to zero and 1
   [~, Observed, ~] = SVMClassifier_Binary([temp_data,temp_resp],cvp,'holdout',cur_round/10,...
      'nshuf',0,'featureselect','none','optimize',1,'pca',0,'solver',1,'kernel','rbf','numkfold',2);
   auc(comparisons(i,1),comparisons(i,2)) = Observed.AUC;
   
   %get the off diagnol error
   targets = zeros(numel(unique(Observed.CorrectResponse)),numel(Observed.CorrectResponse));
   outputs = zeros(numel(unique(Observed.CorrectResponse)),numel(Observed.CorrectResponse));
   targetsIdx = sub2ind(size(targets), Observed.CorrectResponse', 1:numel(Observed.CorrectResponse));
   outputsIdx = sub2ind(size(outputs), Observed.Predictions', 1:numel(Observed.CorrectResponse));
   targets(targetsIdx) = 1;
   outputs(outputsIdx) = 1;
   [~,temp_cm] = confusion(targets,outputs);
   
   %populate confusion matrix (transpose of temp_cm)
   cm(i,1)=temp_cm(2,1)/sum(temp_cm(2,:));
%    cm(i,2,comparisons(i,1))=temp_cm(1,2)/sum(temp_cm(1,:));
end
cm = round(cm,2);
cm(cm==0)=NaN;



%% Plot the confusion matrix

figure('position',[680, 296, 560, 682]);  
h = heatmap(cm);
% Temporarily change axis units 
originalUnits = h.Units;  % save original units (probaly normalized)
h.Units = 'centimeters';  % any unit that will result in squares
% Get number of rows & columns
sz = size(h.ColorData); 
% Change axis size & position;
originalPos = h.Position; 
% make axes square (not the table cells, just the axes)
h.Position(3:4) = min(h.Position(3:4))*[1,1]; 
if sz(1)>sz(2)
    % make the axis size more narrow and re-center
    h.Position(3) = h.Position(3)*(sz(2)/sz(1)); 
    h.Position(1) = (originalPos(1)+originalPos(3)/2)-(h.Position(3)/2);
else
    % make the axis size shorter and re-center
    h.Position(4) = h.Position(4)*(sz(1)/sz(2));
    h.Position(2) = (originalPos(2)+originalPos(4)/2)-(h.Position(4)/2);
end
% Return axis to original units
h.Units = originalUnits; 
% title({'Motif Activity Discriminates';'Behavioral States'})
xlabel('Correct Behavioral State');
ylabel('Predicted Behavioral State');
set(gca,'FontSize',16,'FontName','Arial')


%%
% Plot the AUC
figure('position',[680, 296, 560, 682]); hold on; 
y = auc(auc~=0);
x = 0.05-rand(1,numel(y))/10;
plot(x,y,'color',[0.5 0.5 0.5],'linestyle','none','marker','.','markersize',20);
line([-0.05,0.05],[nanmean(y),nanmean(y)],'color','k','Linewidth',2);
line([0,0],[nanmean(y)-sem(y),mean(y)+sem(y)],'color','k','Linewidth',2)
line([-1,1],[0.5 0.5],'color',[0.75 0.75 0.75],'Linewidth',2,'Linestyle','--')
ylabel('Classifier AUC');
set(gca,'ylim',[0 1],'xlim',[-0.3 .3],'xtick','')
setFigureDefaults
set(gca,'position',[3 3 4 8.5])
[h,p] = ttest(y,0.5,'Tail','right');
AddSig(h,p,[0 0 0.9 0.9],1,0,1);

%%
stats(cur_round).mean_auc = mean(y);
stats(cur_round).sem_auc = sem(y);
stats(cur_round).targets = targets;
stats(cur_round).outputs = outputs;
stats(cur_round).auc =auc;

end


% % 
% % % Plot the confusion matrix
% % figure('position',[680, 296, 560, 682]); hold on; 
% % plotconfusion(targets,outputs)
% % title({'Motif Activity Discriminates';'Behavioral States'},'FontName','Arial','FontSize',16,'FontWeight','normal')
% % xlabel('Correct Behavioral State');
% % ylabel('Predicted Behavioral State');
% % set(gca,'xticklabelrotation',0)
% % setFigureDefaults;
% % set(gca,'position',[3 3 8.5 8.5])
% 
% % Convert this data to a [numClasses x 6] matrix
% targetsVector = cat(1,targetsVector{:});
% outputsVector = cat(1,outputsVector{:});
% targets = zeros(numel(unique(labels)),numel(targetsVector));
% outputs = zeros(numel(unique(labels)),numel(targetsVector));
% targetsIdx = sub2ind(size(targets), targetsVector', 1:numel(targetsVector));
% outputsIdx = sub2ind(size(outputs), outputsVector', 1:numel(targetsVector));
% targets(targetsIdx) = 1;
% outputs(outputsIdx) = 1;


%    %convert output to original labels
%    temp = NaN(numel(Observed.CorrectResponse),1);
%    temp(Observed.CorrectResponse==1)=comparisons(i,1);
%    temp(isnan(temp))=comparisons(i,2);
%    targetsVector{i} = temp;
%    
%    temp = NaN(numel(Observed.Predictions),1);
%    temp(Observed.Predictions==1)=comparisons(i,1);
%    temp(isnan(temp))=comparisons(i,2);
%    outputsVector{i} = temp;      









