function [n, p] = symbol(T, format)
%SYMBOL Returns the symbol of the Toeplitz part in the CQT matrix T.
%
%   [N, P] = SYMBOL(T) returns two vectors containing the negative and
%       positive coefficients of the symbol of the Toeplitz part of T,
%       respectively.
%
%   F = SYMBOL(T, 'fun') returns a function handle that evaluates the
%   symbol at any point z. 

if ~exist('format', 'var')
    format = 'vector';
end

switch format
    case 'vector'
        n = T.n;
        p = T.p;
    case 'fun'
        n = @(z) evaluate_symbol(T.n, T.p, z);
    otherwise
        error('Unsupported format required for SYMBOL');
end

    function r = evaluate_symbol(n, p, z)
        if isempty(n)
            r = 0.0;
        else
            r = polyval(p(end:-1:1), z) + ...
                polyval([ n(end:-1:2) 0 ], 1 ./ z);
        end
    end

end

