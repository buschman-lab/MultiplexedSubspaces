% Copyright (C) 2016 Arno Onken
%
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

function [W,H,err] = nmf(X,k,...
                         W,b,max_iter,max_tol,max_err,replicates)
% Calculates the nonnegative matrix factorization X = W*H, where X is a
% [n, m] matrix, W is a [n, k] matrix and H is a [k, m] matrix. If W is
% specified then W is fixed and not optimized (for test sets). The
% optimization terminates when the change in H and W is not greater than
% max_tol, after max_iter iterations or when the norm of the error matrix
% is smaller than max_err. The multiplicative update rules are used for
% the optimization (Lee and Seung, 2001).
% Input arguments:
%  X          - Input matrix to be factorized (size n x m)
%  k          - Number of modules
% Optional input arguments:
%  W          - First matrix factor (size n x k); if specified, then W is
%               not optimized
%  b          - Divergence index; Set b=2 for Frobenius norm (default),
%               b=1 for KL-divergence or b=0 for Itakura-Saito divergence
%  max_iter   - Maximum number of iterations after which to stop
%  max_tol    - Maximum tolerance below which to stop
%  max_err    - Maximum error below which to stop
%  replicates - Number of optimization repetitions for avoiding local
%               minima
% Output:
%  W          - First matrix factor (size n x k)
%  H          - Second matrix factor (size k x m)
%  err        - Final reconstruction error

optimize_W = nargin<3 || isempty(W); % If false, do not optimize W
if ~optimize_W && (size(W,1)~=size(X,1) || size(W,2)~=k)
    error('nmf: X, W and k must have compatible sizes.');
end
if nargin<4 || isempty(b)
    % Default: minimize Frobenius norm
    b = 2;
end
if nargin<5 || isempty(max_iter)
    max_iter = 1000;
end
if nargin<6 || isempty(max_tol)
    max_tol = 1e-3;
end
if nargin<7 || isempty(max_err)
    max_err = 0;
end
if nargin<8 || isempty(replicates)
    replicates = 1;
end

n = size(X,1);

err_best = inf;
for j = 1:replicates
    if optimize_W
        W = rand(n,k);
    end
    H = pinv(W)*X;
    H(H<0) = 0;
    i = 0;
    err = inf;
    tol = inf;

    while tol>max_tol && i<max_iter && err>max_err
        H_old = H;

        % Apply multiplicative update rules
        WH = W * H + eps;
        H = H .* (W'*(WH.^(b-2).*X)) ./ (W'*WH.^(b-1) + eps);
        WH = W * H + eps;
        if optimize_W
            W_old = W;
            W = W .* ((WH.^(b-2).*X)*H') ./ (WH.^(b-1)*H' + eps);
            WH = W * H + eps;
        end
        if b==0
            % Minimize Itakura-Saito divergence
            err = sum(X(:) ./ (WH(:)+eps) - log(X(:) ./ (WH(:)+eps)+eps) - 1);
        elseif b==1
            % Minimize KL divergence
            err = sum(X(:) .* log(X(:) ./ (WH(:)+eps)+eps) - X(:) + WH(:));
        else
            % Minimize Frobenius norm
            err = norm(X - WH,'fro');
        end

        i = i+1;
        tol = sum(abs(H(:)-H_old(:)));
        if optimize_W
            tol = tol + sum(abs(W(:)-W_old(:)));
        end
    end
    if err<err_best
        err_best = err;
        W_best = W;
        H_best = H;
    end
end
err = err_best;
W = W_best;
H = H_best;

end

