function [l, u] = reciprocal_cr(vm, vp)
%RECIPROCAL_CR 

accurate = false;

kleft = 1;
kright = 1;

vm = reshape(vm, 1, length(vm));
vp = reshape(vp, 1, length(vp));

threshold = eps * ( norm(vm) + norm(vp) );

while ~accurate
    [~,~,l,u] = spectral_cr([ vm , zeros(1, kleft) ], ...
        [ vp, zeros(1, kright) ]);
    
    if abs(l(end)) < threshold && abs(u(end)) < threshold
        accurate = true;
    end
    
    % Increase the padding to gain more accuracy
    kleft  = 2 * kleft;
    kright = 2 * kright;
end

l = cln(l).';
u = cln(u).';

end

