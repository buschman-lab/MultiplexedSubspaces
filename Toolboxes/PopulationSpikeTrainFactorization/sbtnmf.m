% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <http://www.gnu.org/licenses/>.

% This function is based on the function sNM3F_basic by Ioannis Delis and
% Bastien Berret.
%--- Contact the authors
%  Ioannis Delis (id2286@columbia.edu)
%  Bastien Berret (bastien.berret@u-psud.fr) 
%  Arno Onken (arno.onken@iit.it)

function [Acal,Wi,Wb,vaf,err] = sbtnmf(X,n_tm,n_sm,n_e,...
                                       Wi,Wb,max_iter,err_tol,display_iter,n_stops)
% Space-by-time non-negative matrix factorization.
% Input arguments:
%  X            - Input matrix
%                 (size #cells x #time_points_per_sample * #samples)
%                 composed of horizontally concatenated data samples
%  n_tm         - Number of temporal (or row) modules
%  n_sm         - Number of spatial (or column) modules
%  n_e          - Number of samples (i.e. #trials over all stimuli)
% Optional input arguments:
%  Wi           - n_tm temporal modules (not optimized if Wi and Wb specified)
%  Wb           - n_sm spatial modules (not optimized if Wi and Wb specified)
%  max_iter     - Maximum iteration before stopping is enforced
%  err_tol      - Tolerance on the recontruction error changes
%  display_iter - Display or not some information at each iteration
%  n_stops      - Number of steps for the stopping criterion
% Output:
%  Acal         - Activation coefficients
%  Wi           - n_tm temporal modules
%  Wb           - n_sm spatial modules
%  vaf          - Variance accounted for
%  err          - Total reconstruction error (across iterations)

% Check whether we should optimize modules
optimize_modules = nargin<6 || isempty(Wi) || isempty(Wb);
% Set default parameters
if nargin<7 || isempty(max_iter)
    max_iter=1000; % Maximum iteration before stopping is enforced
end
if nargin<8 || isempty(err_tol)
    err_tol=1e-10; % Tolerance on the recontruction error changes
end
if nargin<9 || isempty(display_iter)
    display_iter=false; % Display or not some information at each iteration
end
if nargin<10 || isempty(n_stops)
    n_stops=5; % Number of steps for the stopping criterion
end

Mb = X';

if any(Mb(:)<0)
	error('sbtnmf: Please check your input matrix... It should only contain non-negative entries');  
end
if rem(size(Mb,1),n_e)==0
	T=size(Mb,1)/n_e; % number of temporal dimensions (e.g. time frames)
else
   error(['sbtnmf: Please check your input matrix... It should be [N, T*n_e] with ' ...
          'N=number of spatial units, n_e=number of samples and T=number of time steps']);
end
N=size(Mb,2); % number of spatial dimensions (e.g. neurons)

% Get the block transpose of Mb
%%
Mi=block_transpose(Mb,'r',n_e);
%%
% Start space-by-time NMF algorithm
% Initialization of arrays - YOU CAN EDIT THE INITIAL GUESS
if optimize_modules
    Wb=rand(n_sm,N);
    Wi=rand(T,n_tm);
%     disp(['sbtnmf: Starting to extract ' int2str(n_tm) ' temporal and ' ...
%            int2str(n_sm) ' spatial modules...']);
end
A=rand(n_tm,n_sm*n_e);  

% Error container
err_all=NaN(1,max_iter);

% Initialize some variables
count=0;
it=0;

% Main iterative loop of the algorithm  
% N.B.: we decompose the data as Mi=Wi*A*Wb (for each sample)
while count<n_stops && it<max_iter
    it=it+1; 
    % Clean the Output, by reordering or rescaling etc.
    % May improve the robustness & convergence
    if optimize_modules
        %%
        % Norm rows and colums to one
        [A,Wi,Wb] = normalize_output(A,Wi,Wb,n_sm,n_tm,n_e);
        
        %%
        % Then, order the elements
        [A,Wi,Wb] = order_output(A,Wi,Wb,n_tm,n_sm,n_e);
    
        %%
        % Update Wb 
        % We approximate Mb=Mi^{\prime}=Cb^{\prime}*Wb=Cbb*Wb
        Cb=Wi*A;
        Cbb=block_transpose(Cb,'c',n_e);
        %%
        numer=Cbb'*Mb;

        denom=(Cbb'*Cbb)*Wb;
        Wb= Wb.*numer./(denom+eps(numer));

%%
        % Update Wi 
        Atil=block_transpose(A,'c',n_e);
        Ci=Atil*Wb;
        
        %%
        Cii=block_transpose(Ci,'r',n_e);
        numer=Mi*Cii';
        
        %%

        denom=Wi*(Cii*Cii');
        Wi= Wi.*numer./(denom+eps(numer));
  
    else
        %%
        Atil=block_transpose(A,'c',n_e);
        %%
    end

    % Update A for all samples s=1..S and compute the error err(it)
    err_all(it)=0;
    for s=1:n_e
    	denom=(Wi'*Wi)*Atil(n_tm*(s-1)+1:n_tm*s,:)*(Wb*Wb');
        numer=Wi'*Mi(:,N*(s-1)+1:N*s)*Wb';
    	A(:,n_sm*(s-1)+1:n_sm*s)=A(:,n_sm*(s-1)+1:n_sm*s).*(numer./(denom+eps(numer)));
    	err_all(it)=err_all(it)+norm(Mi(:,N*(s-1)+1:N*s)-Wi*A(:,n_sm*(s-1)+1:n_sm*s)*Wb,'fro')^2;
    end
    
	if display_iter
    	disp(['Iteration: ' num2str(it) '; error: ' num2str( err_all(it))]);
	end
                 
    % Implement convergence criterion
    if it>2
        if abs(err_all(it)-err_all(it-1))<err_tol
            count=count+1;
            if display_iter
                disp(['stop counter #' num2str(count)]);
            end
        elseif err_all(it)>err_all(it-1)
            break;
        else
            count=0;
        end
    end
end

SST=0;
for s=1:n_e
	SST=SST+norm(Mi(:,N*(s-1)+1:N*s)-mean(mean(Mi(:,N*(s-1)+1:N*s))),'fro')^2;
end

err=err_all(it);
vaf=1-err/SST;

if optimize_modules
    % Clean the Output, by reordering or rescaling etc.
    % Norm rows and colums to one
    %%
    [A,Wi,Wb] = normalize_output(A,Wi,Wb,n_sm,n_tm,n_e);
    %%
    % Then, order the elements
    [A,Wi,Wb] = order_output(A,Wi,Wb,n_tm,n_sm,n_e);
end

% Build the full Acal cell array of activations
Acal=zeros(n_tm,n_sm,n_e);
for s=1:n_e
	Acal(:,:,s)=A(:,n_sm*(s-1)+1:s*n_sm);
end

end % sbtnmf


function Mat_out=block_transpose(Mat_in,type,n_e)

if strcmpi(type,'r') % in rows, Mat_in is [T*n_e,N]
    [r,c]=size(Mat_in);
    D=r/n_e; % equal to T
    Mat_out=zeros(D,c*n_e); % We are going to get a [T,N*n_e] matrix
    for s=1:n_e
        Mat_out(:,c*(s-1)+1:c*s)=Mat_in(D*(s-1)+1:D*s,:);
    end
elseif strcmpi(type,'c') % in columns, Mat_in is [T,N*n_e]
    [r,c]=size(Mat_in);
    D=c/n_e; % equal to N
    Mat_out=zeros(r*n_e,D); % We are going to get a [T*n_e,N] matrix
    for s=1:n_e
    	Mat_out(r*(s-1)+1:r*s,:)=Mat_in(:,D*(s-1)+1:D*s); 
    end
else
    error('sbtnmf:block_transpose: Unknown option. It is either "r" or "c".')
end

end % block_transpose


function [A,Wi,Wb] = normalize_output(A,Wi,Wb,n_sm,n_tm,n_e)

% Row-wise normalization
sP=zeros(1,n_tm);
for i=1:n_tm
   sP(i)=norm(Wi(:,i));
   Wi(:,i)=Wi(:,i)./sP(i);
end
% Column-wise normalization
sN=zeros(1,n_sm);
for j=1:n_sm
   sN(j)=norm(Wb(j,:));
   Wb(j,:)=Wb(j,:)./sN(j);
end
% Element-wise normalization of A to make the error unchanged
for s=1:n_e
    for i=1:n_tm
       for j=1:n_sm  
        A(i,n_sm*(s-1)+j)=A(i,n_sm*(s-1)+j).*(sP(i)*sN(j));
       end
    end
end

end % normalize_output


function [Anew,Wi,Wb] = order_output(A,Wi,Wb,n_tm,n_sm,n_e)

Anew=zeros(size(A));

% Order elements by the occurence of the maximal element
% in columns for Wi
iP=zeros(1,n_tm);
for i=1:n_tm
    [~,iP(i)]=max(Wi(:,i));
end
[~,idxP] = sort(iP,'ascend');

% in rows for Wb
iN=zeros(1,n_sm);
for j=1:n_sm
    [~,iN(j)]=max(Wb(j,:));
end
[~,idxN] = sort(iN,'ascend');

Wi=Wi(:,idxP);
Wb=Wb(idxN,:);

for s=1:n_e
    for j=1:n_sm
        Anew(:,n_sm*(s-1)+j)=A(:,n_sm*(s-1)+idxN(j));
    end
end
Anew(:,:)=Anew(idxP,:);

end % order_output
