function [b, num_steps] = brent_zero(a, b, t, f, options)
% Brent's method for root finding

% INPUTS:
% a = lower bound on root
% b = upper bound on root
% t = tolerance
% f = implicit function containing root in (a,b)
% options
%   .maxiter = max iterations

% Example:
%{
f = @(x) (x+3)*(x-1)^2;
a = -4;
b = 4/3;
[b, num_steps] = brent_zero(a, b, '', f)
%}

% Initialize parameters
c = a;
fa = f(a);
fb = f(b);
fc = f(c);
flag = true;
num_steps = 0;
maxiter = options.maxiter;

% Check parameters
if (fa*fb >= 0) || isnan(fa) || isnan(fb)
    error('function must have differing sign on interval points');
end
if isempty(t)
    t = eps;
end

% Iterate
while abs(fb) > t && a ~= b && num_steps < maxiter
    if (fa ~= fc) && (fb ~= fc)
        dab = fa-fb;
        dac = fa-fc;
        dbc = fb-fc;
        s =     a*fb*fc / (dab*dac);
        s = s + b*fa*fc / ((-dab)*dbc);
        s = s + c*fa*fb / ((-dac)*(-dbc));
    else
        s = b - fb * (b-a) / (fb-fa);
    end
    
    if ((s<(3*a+b)/4) || s>b || ( flag && (abs(s-b) >= abs(b-c)/2)) || (~flag && (abs(s-b) >= abs(c-d)/2)) || ( flag && (abs(b-c) < t)) || (~flag && (abs(c-d) < t)))
        
        s = (a+b)/2;
        flag = true;
    else
        flag = false;
    end
    
    fs = f(s);
    d = c;
    c = b;
    fc = f(c);
    if fa*fs < 0
        b = s;
        fb = f(b);
    else
        a = s;
        fa = f(a);
    end
    if abs(fa) < abs(fb)
        temp = b;
        b = a;
        a = temp;
        fa = f(a);
        fb = f(b);
    end
    
    num_steps = num_steps+1;
    if num_steps == maxiter
        fprintf('max iterations reached\n');
    end
end
        
        
           
        
        