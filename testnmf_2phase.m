%%
n=1000;m=100;
% rng(1);
seed=rng;
rng(seed);
K=8;
X00=rand(n,K);Y00=rand(K,m);M=max(0,X00*Y00+randn(n,m)*1e-1);
% M=rand(n,m);
svd(M)'
k=K;
X0=rand(n,k);Y0=rand(k,m);
%%
tol=1e-6;
maxtime=30;
[X,W,iter,HIS]=nmf(M,K,'MAX_TIME',maxtime,'W_INIT',X0,'H_INIT',Y0,'TOL',tol);
[X,W,iter,elapse,hisne]=NeNMF(M,K,'MAX_TIME',maxtime,'W_INIT',X0,'H_INIT',Y0,'TOL',tol);
[X,W,iter,elapse,hisbb]=NMF_QRPBB(M,K,'MAX_TIME',maxtime,'W_INIT',X0,'H_INIT',Y0,'TOL',tol);
[X,W,histtime,histg,histf]=nmf_twophase( M,X0,Y0 );
figure;semilogy(hisne.cpus,hisne.prjg,':b',hisbb.t,hisbb.p,'--r',HIS(:,3),HIS(:,9),'-.g',histtime,histg,'-m');
legend('NeNMF','QRPBB','ANLS-BPP','2-STAGE');xlim([0,maxtime]);xlabel('cpu time(s)');ylabel('E');
% figure;semilogy(hisne.cpus,hisne.objf,hisbb.t,hisbb.f,HIS(:,3),HIS(:,6),histtime,histf);
hisne.objf(end),hisbb.f(end),histf(end)