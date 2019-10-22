function C = qt_mldivide(A, B)

% Handle the triangular case
if length(A.p) == 1 || length(A.n) == 1
    C = inv(A) * B;
    return;
end

[U, L, E] = ul(A);

c = limit(E)';
UE = E.U; VE = E.V;

Uinv = inv(U);
Linv = inv(L);

UB = Uinv * B;

if ~isempty(c)
    if size(VE, 1) < length(c)
        VE = [ VE ; zeros(length(c) - size(VE, 1), size(VE, 2)) ];
    elseif size(VE, 1) > length(c)
        c(size(VE, 1)) = 0;
        
    end
    
    Vh = [ VE, c(:) ];
else
    Vh = VE;
end

[ku1, ku2] = size(UE);
[kv1, kv2] = size(VE);

UE0 = Uinv  * cqt(UE); 
VE0 = Linv' * cqt(Vh);

if ~isempty(c)
    kv2 = kv2 + 1;
    gamma = sum(Uinv.p);
end

UE0 = slicemat(UE0, { 1:ku1, 1:ku2 });
VE0 = slicemat(VE0, { 1:kv1, 1:kv2 });

if size(UE0, 1) < kv1
    UE0 = [ UE0 ; zeros(kv1 - size(UE0, 1), size(UE0, 2)) ];
end

if ~isempty(c)
    S = eye(kv2) + VE0' * [ UE0(1:kv1,:), gamma * ones(kv1, 1) ];
else
    S = eye(kv2) + VE0' * UE0(1:kv1,:);
end

VE0 = VE0 / S';

if ~isempty(c)
    C = cqt('extended', 1, 1, -UE0, VE0(:,1:end-1), -gamma * VE0(:,end)');
else
    C = cqt(1, 1, -UE0, VE0);
end

C = Linv * (C * UB);

end
