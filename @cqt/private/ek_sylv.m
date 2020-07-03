function [Xu, Xv, VA, VB] = ek_sylv(A, B, u, v, k, tol, debug, nrm_type)
%EK_SYLV Approximate the solution of a Sylvester equation AX + XB + U*V' = 0.
%
% [XU,XV] = EK_SYLV(A, B, U, V, K) approximates the solution of the
%     Sylvester equation in the factored form XU * XV'.
%
% [XU, VA] = EK_SYLV(A, B, U, V, K, TOL, DEBUG) also returns the bases VA
%     and VB, and the optional parameters TOL and DEBUG control the
%     stopping criterion and the debugging during the iteration.
%
% The tolerance TOL can also be specified as a function TOL(R, N) that
% takes as input the residual norm and the norm of the solution (R and N,
% respectively), and returns true if the solution can be accepted.

if ~exist('debug', 'var')
    debug = false;
end

if ~exist('tol', 'var')
    tol = 1e-8;
end

if ~exist('nrm_type', 'var')
    nrm_type = 2;
end

nrmA = norm(A, nrm_type);
nrmB = norm(B, nrm_type);

% Dimension of the space
sa = size(u, 2);
sb = size(v, 2);

bsa = sa;
bsb = sb;

% tol can be function tol(r, n) that is given the residual and the norm, or
% a scalar. In the latter case, we turn it into a function
if isfloat(tol)
    tol_eps = tol;
    tol = @(r, n) r < tol_eps * n;
end

it=1;
while max(sa-2*bsa, sb-2*bsb) < k
    
    if ~exist('VA', 'var')
        [VA, KA, HA, param_A] = ek_krylov(A, u);
        [VB, KB, HB, param_B] = ek_krylov(B', v);
    else    
        [VA, KA, HA, param_A] = ek_krylov(VA, KA, HA, param_A);
        [VB, KB, HB, param_B] = ek_krylov(VB, KB, HB, param_B);
    end
    
    sa = size(VA, 2);
    sb = size(VB, 2);
    
    % Compute the solution and residual of the projected Lyapunov equation
    As = HA / KA(1:end-bsa,:);
    Bs = HB / KB(1:end-bsb,:);
    Cs = zeros(size(As, 1), size(Bs, 1));

    if ~exist('Cprojected', 'var')
        Cprojected = ( VA(1:size(u,1),1:size(u,2))' * u ) * ( v' * VB(1:size(v,1),1:size(v,2)) );
    end
    
    Cs(1:size(u,2), 1:size(v,2)) = Cprojected;
    
    [Y, res] = lyap_galerkin(As, Bs, Cs, bsa, bsb);
    
    % You might want to enable this for debugging purposes
    if debug
        Xu = VA(:,1:size(Y,1)) * Y;
        Xv = VB(:,1:size(Y,2));
        X = cqt([], [], Xu, Xv);
        C = cqt([], [], u, v);
        fprintf('%d %e\n', 2*it-1, norm(A*X + X*B + C, 2));
        fprintf('%d %e\n', 2*it, norm(A*X + X*B + C, 2));
    end
    
    if tol(res, norm(Y, nrm_type)) % res < norm(Y) * tol
        break
    end
    it=it+1;
end

[UU,SS,VV] = svd(Y);

% rk = sum(diag(SS) > SS(1,1) * tol / max(nrmA, nrmB));
rk = sum(arrayfun(@(s) tol(s, SS(1,1) / max(nrmA, nrmB)), diag(SS)) == false);

Xu = VA(:,1:size(Y,1)) * UU(:,1:rk) * sqrt(SS(1:rk,1:rk));
Xv = VB(:,1:size(Y,2)) * VV(:,1:rk) * sqrt(SS(1:rk,1:rk));

end

function [Y, res] = lyap_galerkin(varargin)
%LYAP_GALERKIN Solve the reduced Lyapunov equation and check the residual.
%
% Y = LYAP_GALERKIN(HA, HB, C, bsa, bsb) solves the (projected) Lyapunov
%        equation HA1 * Y + Y * HB1' = C1, where HA1 is obtained shrinking
%        HA by BSA columns and rows, HB1 is obtained shrinking HB1 by BSB
%        columns and rows, and C1 is obtained cutting C accordingly.
%
% [Y, RES] = LYAP_GALERKIN(HA, HB, C, bsa, bsb) does the same computation
%        but also computes the residual of the unprojected equation,
%        assuming that HA, HB, and C have been projected using a
%        Krylov-type method and that the action of A on the basis excluding
%        the last BSA (resp. BSB) columns is contained in the full basis.

if length(varargin) == 3
    HA = varargin{1};
    C = varargin{2};
    bsa = varargin{3};
    is_lyapunov = true;
else
    HA = varargin{1};
    HB = varargin{2};
    C = varargin{3};
    bsa = varargin{4};
    bsb = varargin{5};
    is_lyapunov = false;
end

% Consider the projected matrices at the previous step, which is needed to
% check the Galerkin condition
HA1 = HA(1 : end - bsa, :);

if ~is_lyapunov
    HB1 = HB(1 : end - bsb, :);
end

% Compute the solution of the Lyapunov equation (word of warning: please
% check the sign of C in the implementation of SylvKrylov).
if ~is_lyapunov
    Y = lyap(HA1, HB1', C(1:end-bsa,1:end-bsb));
else
    Y = lyap(HA1, C(1:end-bsa,1:end-bsa));
end

% Check the residual
if ~is_lyapunov
    res = max(norm(HA(end-bsa+1:end, :) * Y), ...
        norm(Y * HB(end-bsb+1 : end, :)'));
else
    res = max(norm(HA(end-bsa+1:end, :) * Y), ...
        norm(Y * HA(end-bsa+1 : end, :)'));
end

end
