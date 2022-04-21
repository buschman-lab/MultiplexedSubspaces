function [stPC, rho] = PopulationCoupling(x,type)
%Camden MacDowell

%computes the spike-triggered coupling or pearsons coupling of each neuron in the neuron x time
%matrix x.

%stPC follows "Diverse coupling of neurons to populations in
%sensory cortex" 2015 - Nature


[s,~] = size(x);
stPC = NaN(s,1);
rho = NaN(s,1);
for i = 1:s
    neu_idx = ~ismember(1:s,i);
    if strcmp(type,'rho')        
        rho(i,:) = nanmean(corr(x(i,:)',x(neu_idx,:)','type','pearson'));
    else strcmp(type,'stPC')
        popFR = sum( x(neu_idx,:)-nanmean(x(neu_idx,:),2) );
        stPC(i,:) = (1/sum(x(i,:))) * trapz(x(i,:) .* popFR); 
    end
end


end