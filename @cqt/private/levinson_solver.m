function x = levinson_solver(c,r,b)
%LEVINSON Solve a Toeplitz linear system.

if size(b, 2) > 1
    x = zeros(length(c), size(b, 2));
    for i = 1 : size(x, 2)
        x(:,i) = levinson_solver(c, r, b(:,i));
    end
    return;
end

c = reshape(c, length(c), 1);
r = reshape(r, length(r), 1);

s = c(2:end);
r = r(2:end);
rho0 = c(1);

n=length(s);
n=n+1;

x=zeros(n,1);
y=x;
z=x;
ga=rho0;
eta=-r(1)/ga;
phi=-s(1)/ga;

x(1)=b(1)/ga;
y(1)=-r(1)/ga;
z(1)=-s(1)/ga;

for k=1:n-1
    ga=(1-eta*phi)*ga;
    alpha = ( b(k+1)-s(1:k)'*flipud(x(1:k)))/ga;
    x(1:k)=x(1:k)+flipud(y(1:k))*alpha;
    x(k+1)=alpha;
    if (k<n-1)
        eta   = (-r(k+1)-r(1:k)'*flipud(y(1:k)))/ga;
        phi   = (-s(k+1)-s(1:k)'*flipud(z(1:k)))/ga;
        yold=y;
        y(1:k)=y(1:k)+flipud(z(1:k))*eta;
        y(k+1)=eta;
        z(1:k)=z(1:k)+flipud(yold(1:k))*phi;
        z(k+1)=phi;
    end
end
