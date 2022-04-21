function [ndim,y] = SVCA(x,npc,cutoff,tpoint)
%Camden MacDowell - timeless
%wrapper for stringer et., science 2019
if nargin <2; npc = floor(size(x,1)/2); end
if nargin <3; cutoff = []; end
if nargin <4; tpoint = []; end

rng('default');
npc = min(npc,floor(size(x,1)/2));
ntrain = randperm(size(x,1),floor(size(x,1)/2));
ntest = find(ismember(1:size(x,1),ntrain)==0);

tidx = randperm(size(x,3),floor(size(x,3)/2));

t = ones(size(x,2),size(x,3));
t(:,tidx)=0;
t = t(:);
itrain = find(t==1);
itest = find(t==0);

xx = reshape(x,[size(x,1),size(x,2)*size(x,3)]);
xx = xx-nanmean(xx,2);
[sneur, varneur, ~, ~, s1, s2]=runSVCA(xx, npc, ntrain, ntest, itrain, itest);

if ~isempty(tpoint) %just use a specific timepoint within the trial
   s1 = reshape(s1,size(s1,1),size(x,2),size(s1,2)/size(x,2));
   ss1 = s1(:,tpoint,:);
   ss1 = reshape(ss1,[size(ss1,1),size(ss1,2)*size(ss1,3)]);
   ss2 = s2(:,tpoint,:);
   ss2 = reshape(ss2,[size(ss2,1),size(ss2,2)*size(ss2,3)]);   
   sneur = sum(ss1 .* ss2, 2);
   varneur = sum(ss1.^2 + ss2.^2,2)/2;
end

if ~isempty(cutoff) %set a threshold in reliable variance
    y = cumsum(sneur./varneur);
    y = y/max(y); 
    ndim = find(y>=cutoff,1,'first');
else %get AUC
    y = sneur./varneur;
    %stop at first negative
    y(find(y<0,1,'first'):end)=[];
    y = cumsum(y)/max(cumsum(y));
    y = [0, y', 1];
    ndim = 1-trapz(y/numel(y));
end

y = cumsum(sneur./varneur);

end %function end


function [sneur, varneur, u, v, s1, s2] = runSVCA(Ff, npc, ntrain, ntest, itrain, itest)
%from stringer et., science 2019
cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
[u,~,v] = svdecon(cov);
u = u(:,1:npc);
v = v(:,1:npc);
s1 = u' * Ff(ntrain,itest);
s2 = v' * Ff(ntest,itest);
sneur = sum(s1 .* s2, 2);
varneur = sum(s1.^2 + s2.^2,2)/2;

end

% INPUTS: 
%     Ff (neurons x timepts)
%     npc (number of PCs to compute)
%     ntrain (one half of neurons)
%     ntest (other half of neurons)
%     itrain (one half of timepts)
%     itest (other half of timepts)
% OUTPUTS:
%     sneur (shared variance of each covariance component)
%     vneur (total variance of each covariance component)
%     u (left eigenvectors of covariance matrix btw ntrain and ntest on
%        itrain timepts)
%     v (right eigenvectors of covariance matrix btw ntrain and ntest on
%        itrain timepts)