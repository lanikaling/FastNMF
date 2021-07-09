Data=AgMorData(:,2:end);
rr=AgMorData(:,1);
fileNum=size(Data,2);
% Data=Data(10:2990,:);rr=rr(10:2990);
figure;plot(rr,Data(:,1:end));
figure;plot(rr,Data+repmat((0:35)*2,2999,1),'k');xlabel('r');ylabel('G')
%%
K=3;
[U,S,V]=svd(Data);U=U(:,1:K);S=S(1:K,1:K);V=V(:,1:K);
GBc=repmat(min(Data),size(Data,1),1);
%%
randX=rand(size(U));randW=rand(size(V'));

maxtime=10;
[X,W,iter,HIS]=nmf(Data-GBc,K,'MAX_TIME',maxtime,'W_INIT',randX,'H_INIT',randW);
[X,W,iter,elapse,hisne]=NeNMF(Data-GBc,K,'MAX_TIME',maxtime,'W_INIT',randX,'H_INIT',randW);
[X,W,iter,elapse,hisbb]=NMF_QRPBB(Data-GBc,K,'MAX_TIME',maxtime,'W_INIT',randX,'H_INIT',randW);
[X,W,histtime,histg,histf]=nmf_twophase( Data-GBc,randX,randW );
figure;semilogy(hisne.cpus,hisne.prjg,':b',hisbb.t,hisbb.p,'--r',HIS(:,3),HIS(:,9),'-.g',histtime,histg,'-m');
legend('NeNMF','QRPBB','ANLS-BPP','2-STAGE');xlim([0,maxtime]);xlabel('cpu time(s)');ylabel('E');
% figure;semilogy(hisne.cpus,hisne.objf,hisbb.t,hisbb.f,HIS(:,3),HIS(:,6),histtime,histf);


