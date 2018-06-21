function [Uf,Vf]=lanczos_hankel_product(a,b)
% Function that find a good low-rank approximation of the product
% between the two anti-triangular Hankel matrices represented with the
% vectors a and b, by means of the Lanczos bidiagonalization method applied with
% starting vector e1 (first vector of the canonical basis)

if size(a,1)>1
    a=a.';
end
if size(b,1)>1
    b=b.';
end
if(sum(abs(a)>eps)==0 || sum(abs(b)>eps)==0)
    Uf=0; Vf=0;
    return
end
m=length(a); n=length(b);
h=min(m,n);
k_max=50;
alpha=zeros(k_max,1); betha=zeros(k_max-1,1);
U=zeros(m,k_max);
V=zeros(n,k_max);
omega=1;

U(:,1)=toepmult_fft(a(end:-1:1),a(end),m,h,b(1:h).'); % we choose as starting poi the first column of H(a)*H(b)
U(:,1)=U(end:-1:1,1);
%U(:,1)=ones(m,1);
%U(:,1)=eye(m,1);
%U(:,1)=randn(m,1);
U(:,1)=U(:,1)/norm(U(:,1));

V(1:h,1)=toepmult_fft(conj(a(h:-1:1)),conj(a(h:end)),h,m,U(:,1));
V(1:h,1)=V(h:-1:1,1);
V(:,1)=toepmult_fft(conj(b(end:-1:1)),conj(b(end)),n,h,V(1:h,1));
V(:,1)=V(end:-1:1,1);
alpha(1)=norm(V(:,1));
%tol=eps*alpha(1);
%tol = eps*alpha(1)*m;
tol=eps*sqrt([1:m]*abs(a.').^2)*sqrt([1:n]*abs(b.').^2)*10;

V(:,1)=V(:,1)/alpha(1);

it=2;
while   1 %&& omega>tol
    omega_old=omega;
    if it>k_max
        U=[U,zeros(m,k_max)];
        V=[V,zeros(n,k_max)];
        alpha=[alpha;zeros(k_max,1)];
        betha=[betha;zeros(k_max-1,1)];
        k_max=2*k_max;
    end
    
    U(1:h,it)=toepmult_fft(b(h:-1:1),b(h:end),h,n,V(:,it-1));
    
    U(1:h,it)=U(h:-1:1,it);
    
    U(:,it)=toepmult_fft(a(end:-1:1),a(end),m,h,U(1:h,it));
    
    U(:,it)=U(end:-1:1,it)-alpha(it-1)*U(:,it-1);
    
    U(:,it)=U(:,it)-U(:,1:it-1)*(U(:,1:it-1)'*U(:,it)); % re-orthogonalization step
    U(:,it)=U(:,it)-U(:,1:it-1)*(U(:,1:it-1)'*U(:,it));
    betha(it-1)=norm(U(:,it));
    
    U(:,it)=U(:,it)/betha(it-1);
    
    V(1:h,it)=toepmult_fft(conj(a(h:-1:1)),conj(a(h:end)),h,m,U(:,it));
    V(1:h,it)=V(h:-1:1,it);
    
    V(:,it)=toepmult_fft(conj(b(end:-1:1)),conj(b(end)),n,h,V(1:h,it));
    V(:,it)=V(end:-1:1,it)-betha(it-1)*V(:,it-1);
    
    V(:,it)=V(:,it)-V(:,1:it-1)*(V(:,1:it-1).'*V(:,it)); % re-orthogonalization step
    V(:,it)=V(:,it)-V(:,1:it-1)*(V(:,1:it-1).'*V(:,it));
    alpha(it)=norm(V(:,it));
    V(:,it)=V(:,it)/alpha(it);
    %  [it, betha(it-1),alpha(it), tol]
    if(betha(it-1)<tol && alpha(it)<tol)
        break
    end
    
    it=it+1;
    
    fprintf('Iteration %d, beta = %e\n', it, alpha(it));
    
    
    
end

it=it-1;

[Uf,S,Vf]=svd(diag(alpha(1:it))+diag(betha(1:it-1),-1));

Uf=U(:,1:it)*Uf*S;
Vf=V(:,1:it)*Vf;
