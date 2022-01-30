function [theta,sigstat] = AngleBetweenWeights(x,y,sig_flag)
%Camden - timeless
%designed for use within CommunicationSubspace analysis
%compute the angle between neuron weightings x and y
%x and y can be vectors or matrices. 
%all permutations computed on x

%%
%normalize to max value
[x,y] = normToMax(x,y);
%get angle
theta = getAngle(x,y); 

switch sig_flag %how do you want to compute significance
    case 'bootstrap'
        error('not yet supported')
        sigstat = []; %confidence interval
    case 'perm_weight' %test if smaller than permutated distribution
        rng('default')
        n_perm = 1000;
        theta_perm = NaN(1,n_perm);
        for i = 1:n_perm
            %randomly permute the neural weighings (maintain any structure across subspaces if matrix
            idx = randperm(size(x,1),size(x,1));
            x_perm = x(idx,:);
            x_perm = normToMax(x_perm);
            theta_perm(i) = getAngle(x_perm,y);
        end
        %tailed significance        
        sigstat(1,1) = sum([theta_perm,theta]<=theta)/numel([theta_perm,theta]);
        sigstat(1,2) = sum([theta_perm,theta]>=theta)/numel([theta_perm,theta]);
    case 'none'        
        theta = getAngle(x,y);
        sigstat = [];
    otherwise
        error('unknown significance test option');
end
%%
    

end %function end

function theta = getAngle(x,y)
%vectors 
if (size(x,2)+size(y,2))==2
    theta = acosd( dot(x,y) / ( norm(x)*norm(y) ) );
    %adjust for random sign attach of CCA
    if theta>90
       theta = 180-theta; 
    end
else
    %subspace only considers between pi/2 so no sign correction needed
    theta = rad2deg(subspace(x,y));
end
end

function [x,y] = normToMax(x,y)
if nargin <2; y=[]; end
for i = 1:size(x,2)
   x(:,i) = x(:,i)/max(x(:,i));
end
if ~isempty(y)
    for i = 1:size(y,2)
       y(:,i) = y(:,i)/max(y(:,i));
    end    
else
    y=[];
end
end