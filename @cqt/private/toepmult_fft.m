function u = toepmult_fft(am, ap, m, n, v)
%TOEPMULT_FFT Fast multiplication of a CQT matrix times a vector. 
%
% U = TOEPMULT_FFT(AM, AP, M, N, V) computes the matrix-matrix product M*V
%     where M is a Toeplitz matrix defined by the vectors AM and AP and has
%     size M x N, where N = length(V). 
%

  realflag=isreal(am)*isreal(ap)*isreal(v);
  n1=length(am);n2=length(ap);
  
  if size(v,1) == n
    N= max(m+min(n2,n),n+min(m,n1));
    N=2^ceil(log(N)/log(2));
    a=zeros(1,N);
    mn1=min(n1,m);
    mn2=min(n,n2);
    a(1:mn1)=am(1:mn1); a(end:-1:end-mn2+2)=ap(2:mn2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w=zeros(N,size(v,2));w(1:n,:)=v;
    wf=fft(w);af=fft(a);
    u=ifft(diag(sparse(af))*wf);
    u=u(1:m,:);
    if(realflag == 1)
    	u=real(u);
    end
  else
    disp('in toepmult_fft size mismatch')
    sizeV=size(v)
    sizeA=[m,n]
    return
  end

