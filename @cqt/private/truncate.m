function v=truncate(u,n)
% If the length of u is greater than n then truncate u
v=u;
if (length(u)>n)
    v=v(1:n);
end
