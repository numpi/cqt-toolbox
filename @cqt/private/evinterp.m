function [fm, fp] = evinterp(f, tol, varargin)
%EVINTERP Evaluation / interpolation scheme
%
% This implementation assumes that f maps the real axis into
% itself, so that if real input is given then real output is
% expected. 

if ~exist('tol', 'var')
    tol = cqtoption('threshold');
end

% Number of arguments of the function f
n = length(varargin) / 2;

realflag = true;
for j = 1 : n
    realflag = realflag && isreal(varargin{j});
end

am = varargin(1:2:end);
ap = varargin(2:2:end);

nm = arrayfun(@(i) length(am{i}), 1 : n);
np = arrayfun(@(i) length(ap{i}), 1 : n);

% Estimate the length of the result, so that it is longer then all the
% inputs
km = max(nm) + 1;
kp = max(np) + 1;

converged = false;

while ~converged
    % Build the polynomials of degree (km + kp - 1)
    p = cell(1, n); fp = cell(1, n);
    xi = ifft([ zeros(1, km-1), 1, zeros(1, kp-1) ]);
    for j = 1 : n
        % Make the arrays am, ap long enough
        am{j}(km) = 0; ap{j}(kp) = 0;
        p{j} = [ am{j}(end:-1:1), ap{j}(2:end) ];
        
        fp{j} = ifft(p{j}) ./ xi;
    end
    
    % Evaluation phase -- the function needs to support vectorized form
    % evaluation for this to work. 
    ev = f( fp{:} ) .* xi;
    
    % Interpolation
    ff = fft(ev);
    
    % Check if the decay is strong enough
    if norm(ff(1:ceil(km/2)), 1) + norm(ff(end-ceil(kp/2)+1:end), 1) < ...
            tol
        converged = true;
    else
        km = 2 * km;
        kp = 2 * kp;
    end
end

fm = cln(ff(km:-1:1));
fp = cln(ff(km:end));

if realflag
    fm = real(fm);
    fp = real(fp);
end

end

