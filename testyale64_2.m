load('Yale_64x64.mat');
M=fea(1:165,:)';
svd(M)';

% seed=rng;
rng(seed);
k=3;
% save('seedforyale22.mat','seed');
n=4096;m=size(M,2);

tol=1e-6;
maxtime=60;
ne=zeros(3,10);bb=ne;bpp=bb;two=bb;

for i=1:10

X0=rand(n,k);Y0=rand(k,m);

[X,W,iter,HIS]=nmf(M,K,'MAX_TIME',maxtime,'W_INIT',X0,'H_INIT',Y0,'TOL',tol);
[X,W,iter,elapse,hisne]=NeNMF(M,K,'MAX_TIME',maxtime,'W_INIT',X0,'H_INIT',Y0,'TOL',tol);
[X,W,iter,elapse,hisbb]=NMF_QRPBB(M,K,'MAX_TIME',maxtime,'W_INIT',X0,'H_INIT',Y0,'TOL',tol);
[X,W,histtime,histg,histf]=nmf_twophase( M,X0,Y0 );
% 
ne(1,i)=hisne.cpus(end);ne(2,i)=hisne.prjg(end);ne(3,i)=hisne.objf(end)*1e-7;
bb(1,i)=hisbb.t(end);bb(2,i)=hisbb.p(end);bb(3,i)=hisbb.f(end)*1e-7;
bpp(1,i)=HIS(end,3);bpp(2,i)=HIS(end,9);bpp(3,i)=HIS(end,6)*1e-7;
two(1,i)=histtime(end);two(2,i)=histg(end);two(3,i)=histf(end)*1e-7;
end
%%
format long
[mean(ne,2),min(ne,[],2),max(ne,[],2)]
[mean(bb,2),min(bb,[],2),max(bb,[],2)]
[mean(bpp,2),min(bpp,[],2),max(bpp,[],2)]
[mean(two,2),min(two,[],2),max(two,[],2)]
% save('yale11.mat','ne','bb','bpp','two')









