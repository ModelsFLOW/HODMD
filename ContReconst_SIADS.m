function ContReconst=ContReconst(t,t0,u,deltas,omegas)
[N,M]=size(u);
vv=zeros(M,1);
for m=1:M
 vv(m)=exp((deltas(m)+i*omegas(m))*(t-t0));   
end
ContReconst=u*vv;
