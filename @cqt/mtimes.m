function [ T ] = mtimes(T1, T2)
%MTIMES Multiply two CQT matrices T1 and T2.
%
% T = MTIMES(T1, T2) computes the CQT matrix T obtained multiplying two CQT
%     matrices T1 and T2. If T2 is a dense finite matrix it is considered
%     as a semiinfinite matrix with only the leading top-left corner
%     different from zero.

if isa(T1,'cqt') && isa(T2, 'cqt')
    if T1.sz(2) == T2.sz(1)
        if min(T1.sz(1), T2.sz(2)) == 0
            T = cqt([], [], [], [], T1.sz(1), T2.sz(2));
            return;
        end
        
        if T1.sz(2) == inf
            [cm, cp, cU, cV] = qt_mult(T1.n, T1.p, T1.U, T1.V, ...
                T2.n, T2.p, T2.U, T2.V);
            cm = cm(1:min(length(cm), T1.sz(1)));
            cp = cp(1:min(length(cp), T2.sz(2)));
            cU = cU(1:min(size(cU,1), T1.sz(1)), :);
            cV = cV(1:min(size(cV,1), T2.sz(2)), :);
            T = cqt(cm, cp, cU, cV, [], [], T1.sz(1), T2.sz(2));
        else
            if T1.sz(1) == inf
                rw = max(size(T1.U, 1), length(T1.n) + T1.sz(2));
            else
                rw = T1.sz(1);
            end
            
            if T2.sz(2) == inf
                cl = max(size(T2.V,1), length(T2.p) + T2.sz(1));
            else
                cl = T2.sz(2);
            end
            
            [cm, cp, cU, cV, cW, cZ] = fqt_mult2(T1.n, T1.p, T1.U, T1.V,...
                T1.W, T1.Z, T2.n, T2.p, T2.U, T2.V, T2.W, T2.Z, rw, ...
                T1.sz(2), cl);
            if min(T1.sz(1), T2.sz(2)) == inf
                T = cqt(cm, cp, cU, cV, cW(end:-1:1,end:-1:1), ...
                    cZ(end:-1:1,end:-1:1), rw, cl);
                T = cqt([], [], full(T), [], T1.sz(1), T2.sz(2));
            elseif max(T1.sz(1), T2.sz(2)) == inf
                T = cqt(cm, cp, cU, cV, [], [], T1.sz(1), T2.sz(2));
            else
                T = cqt(cm, cp, cU, cV, cW(end:-1:1,end:-1:1), ...
                    cZ(end:-1:1,end:-1:1), T1.sz(1), T2.sz(2));
            end
            
        end
        
        if ~isempty(T1.c) || ~isempty(T2.c)
            [na, pa] = symbol(T1);
            sa = sum([ na, pa(2:end) ]);
            
            c2 = T2.c; c1 = T1.c;
            T2.c = zeros(1,0); T1.c = zeros(1,0);
            
            T.c = formatted_sum(correction(cqt(c1) * T2), sum(c1) * c2);
            T = extend(T, formatted_sum(T.c, sa * c2));
            
            Te = cumsum(na(end:-1:2)).';
            Te = Te(end:-1:1);
            [U1,V1] = correction(T1);
            T = T + cqt([], [], formatted_sum(-Te, U1 * sum(V1, 1)'), c2');
        end
        
    else
        error('Incompatible inner dimensions');
    end
elseif  isa(T1, 'cqt') && isscalar(T2)
    T = cqt(T1.n * T2, T1.p * T2, T1.U* T2, T1.V, ...
        T1.W(end:-1:1,end:-1:1) * T2, T1.Z(end:-1:1,end:-1:1), ...
        T1.sz(1), T1.sz(2));
    if ~isempty(T1.c)
        T = extend(T, T1.c * T2);
    end
elseif isscalar(T1) && isa(T2, 'cqt')
    T = cqt(T2.n * T1, T2.p * T1, T2.U * T1, ...
        T2.V, T2.W(end:-1:1,end:-1:1) * T1, ...
        T2.Z(end:-1:1,end:-1:1), T2.sz(1), T2.sz(2));
    if ~isempty(T2.c)
        T = extend(T, T2.c * T1);
    end
elseif isa(T1,'cqt') && T1.sz(1) == inf && ~isa(T2, 'cqt')
    error([ 'Incompatible types multiplication. \nIf you want' ...
        ' to multiply a cqt matrix T with a finite matrix of' ...
        ' %s A you can use T * cqt(A) ' ], class(T2));
elseif isa(T1,'cqt') && T1.sz(1) ~= inf && ~isa(T2, 'cqt')
    T = full(T1 * cqt([], [], T2, [], size(T2, 1), size(T2, 2)));
elseif isa(T2,'cqt') && T2.sz(1) == inf && ~isa(T1, 'cqt')
    error([ 'Incompatible types multiplication. \nIf you want to' ...
        'multiply a finite matrix of %s A with a cqt matrix you can' ...
        'use cqt(A) * T' ], class(T1));
elseif isa(T2,'cqt') && T2.sz(1) ~= inf && ~isa(T1, 'cqt')
    T = full(cqt([], [], T1, [], size(T1, 1), size(T1, 2)) * T2);
end

