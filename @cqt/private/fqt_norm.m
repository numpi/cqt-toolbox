function nrm = fqt_norm(am, ap, u, v, w, z)
%FQT_NORM 

nrm = 0.0;

if ~isempty(am)
    nrm = nrm + norm([ am, ap(2:end) ], 1);
end

[~,ru] = qr(u, 0);
[~,rv] = qr(v, 0);
[~,rw] = qr(w, 0);
[~,rz] = qr(z, 0);

nrm = nrm + max([ norm(ru*rv'), norm(rw*rz') ]);

end

