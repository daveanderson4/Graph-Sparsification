function pi = ucs(U, k, l, minc)
%Unweighted Column Selection algorithm
%   deterministic, unweighted algorithm for use in sparse PCA CUR
%   David G. Anderson
%   2014
%
%   inputs:
%       U       orthogonal matrix to select cols
%               typically is V_k of rank k TSVD
%       k       rank of approximation (k<<n)
%       l       max rank of C and R (k<l<<n)
%       minc    if true, finds column to minimize objective
%               if false, finds first column to make objective negative
%               defaults true
%
%   outputs:    
%       pi      index set |pi|=l s.t.
%               sigma_min(U(:,id))>=sqrt(l)-sqrt(k))/sqrt(k)
%
%   notes:
%       roots are found with brent's method for speed and stability

if nargin < 4, minc = true; end

VERBOSE = false;
n = size(U,1);
options.maxiter = 1000;
tol = 1e-6;

% fix inverse trace
T_hat = (k*(n+(l+1)/2-k)+sqrt(k*l*(n-(l-1)/2)*(n+(l+1)/2-k)))/(l-k);
F     = @(t) (1-k/t)*l/(n-(l-1)/2-k+t)-k/t;
T     = T_hat*(1+F(T_hat));
clear F

% initialize vars
A = sparse(k, k);
I = [];
lambda = -Inf;
lambda_hat = -Inf;

mindx = 1;

lbond3 = -n/T;

% main loop to build C and R
for r = 0:l-1
    % step 1: compute lambda so tr(A-lambda*I)^-1=T
    s = eig(A);
    %mins = s(2); % obsolete % first is always zero for sparsifiers
    mins = max(0, min(s));
    
    % inputs for root finder
    fun = @(x) sum(1./(s-x)) - T;
    a   = max(lambda_hat, mins - 1); % extra room for first iteration
    b   = mins - (1e-12);
        
    if fun(a)>0 & fun(b)>0
        [lambda_new, num_steps] = brent_zero(lambda,b,tol,fun,options);
        if lambda_new-lambda_hat < -1e-6
            error('decreasing bound');
        else
            lambda = lambda_new;
        end
    else
        [lambda, num_steps] = brent_zero(a,b,tol,fun,options);
    end

    % step 2: solve for lambda_hat
    c = n-r+sum((1-s)./(s-lambda));
    lh_update = @(lh) (lh-lambda)*c- ...
        sum((1-s)./((s-lambda).*(s-lh)))/ ...
        sum((1)./((s-lambda).*(s-lh)));
    l0 = .5*(lambda+s(k));
    [lambda_hat,num_steps] = brent_zero(lambda,mins-(1e-12),tol,lh_update,options);
    
    % step 3: select i and update
    trinv  = sum(1./(s-lambda));
    trhinv = sum(1./(s-lambda_hat));
    IC = setdiff([1:n],I);
    if minc, adj_max = -1; ix_min = 0; end % any neg. num. to initialize
    
    % continue where last iteration left off
    lenic = length(IC);                         
    mindx= min(mindx,lenic);                        
    nn = mod([1:lenic] + mindx - 2, lenic) + 1;        
    for ii = 1:lenic                                
        ix = IC(nn(ii));                           
        ui  = U(ix, :)'; % transpose of U/V
        d   = 1+ui'*((A-lambda_hat*eye(k))\ui);
        adj = ui'*((A-lambda_hat*eye(k))\((A-lambda_hat*eye(k))\ui))/d;
        
        % testing
        pdrptdp = trhinv-adj - trinv;
        if minc
            if adj > adj_max
                adj_max = adj;
                ix_min  = ix;
            end
        else
            if trhinv-adj <= trinv
                A = A + ui * ui';
                I = [I, ix];
                mindx = ii;         
                break;
            end
        end
    end
    if minc
        ui = U(ix_min, :)';
        A  = A + ui * ui';
        I  = [I, ix_min];
    end
    
    % calculate new lower bound in steps
    lbond=((r-k)*T_hat-k*(n+(r+1)/2-k))/ ...
        (T_hat*(n-(r-1)/2-k+T_hat)+(r-k)*T_hat-k*(n+(r+1)/2-k));
    lbond2=2*k-T-r/2-n;
    lbond2=lbond2+sqrt((n+r/2+T-2*k)+4*T*(-r+r*k/(2*T)+k*n/T+k-k^2/T));
    lbond2=lbond2/(-2*T);
    lbond3=[lbond3;sum((1-lbond3-k/T)./(n-r+(1-lbond3)*T-k))-k/T];
    if (VERBOSE)
        fprintf('%d %12.12f %12.12f %12.12f %12.12f %12.12f %12.12f\n',r,lambda_hat,mins,pdrptdp,lbond,lbond2,lbond3(end));
    end
end

pi = I;

end

