function rho = InternalXValidation(x,y,lambda)
%Camden - timeless
%returns the correlation between betas for cross validations of the same
%subspace between x and y
%input x and y as a 2D matix for timepoitn cross validation or 3D tensor
%for trialwise
if nargin <3; lambda = []; end
ndim = 10; %number of dimensions to compute (def:10)
nxval = 10; %number of cross valdiations (def:10)
rng('default')

if size(x,3)>1 %trialwise
    %10 fold cross validation        
    cvp = cvpartition(size(x,3),'kfold',nxval);

    % Cross validated Betas
    B = cell(1,nxval); V = cell(1,nxval);
    for i = 1:nxval
        [~,B{i},V{i}] = RRR_simple(x(:,:,cvp.training(i)),y(:,:,cvp.training(i)),lambda);              
    end    
else
    %10 fold cross validation        
    cvp = cvpartition(size(x,1),'kfold',nxval);

    % Cross validated Betas
    B = cell(1,nxval); V = cell(1,nxval);
    for i = 1:nxval
        [~,B{i},V{i}] = RRR_simple(x(cvp.training(i),:),y(cvp.training(i),:),lambda);              
    end
end

xx = nchoosek(1:ndim,2);
rho = cell(1,ndim);
for cur_d = 1:ndim
    rho{cur_d} = arrayfun(@(a,b) corr(B{a}(:,cur_d),B{b}(:,cur_d)),xx(:,1),xx(:,2),'UniformOutput',1);
end
rho = cat(2,rho{:});
rho = abs(rho);
        

end %function 

%reason for doing this, instead of split haves is the number of
%observations x predictors... basically, this is more fair

% x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';     
% xx = reshape(xx,[size(xx,1),size(xx,2)*size(xx,3)])';  
% y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';    
% yy = reshape(yy,[size(yy,1),size(yy,2)*size(yy,3)])';    
% x_full = reshape(x_full,[size(x_full,1),size(x_full,2)*size(x_full,3)])';





















