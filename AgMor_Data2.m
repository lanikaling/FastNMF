Data=AgMorData(:,2:end);
rr=AgMorData(:,1);
fileNum=size(Data,2);
% Data=Data(10:2990,:);rr=rr(10:2990);
% figure;plot(rr,Data(:,1:end));
% figure;plot(rr,Data+repmat((0:35)*2,2999,1),'k');xlabel('r');ylabel('G')
%%
K=3;
[U,S,V]=svd(Data);U=U(:,1:K);S=S(1:K,1:K);V=V(:,1:K);
GBc=repmat(min(Data),size(Data,1),1);
%%

ne=zeros(3,10);bb=ne;bpp=bb;two=bb;

maxtime=20;
for i=1:10
randX=rand(size(U));randW=rand(size(V'));

[X,W,iter,HIS]=nmf(Data-GBc,K,'MAX_TIME',maxtime,'W_INIT',randX,'H_INIT',randW);
[X,W,iter,elapse,hisne]=NeNMF(Data-GBc,K,'MAX_TIME',maxtime,'W_INIT',randX,'H_INIT',randW);
[X,W,iter,elapse,hisbb]=NMF_QRPBB(Data-GBc,K,'MAX_TIME',maxtime,'W_INIT',randX,'H_INIT',randW);
[X,W,histtime,histg,histf]=nmf_twophase( Data-GBc,randX,randW );

ne(1,i)=hisne.cpus(end);ne(2,i)=hisne.prjg(end);ne(3,i)=hisne.objf(end);
bb(1,i)=hisbb.t(end);bb(2,i)=hisbb.p(end);bb(3,i)=hisbb.f(end);
bpp(1,i)=HIS(end,3);bpp(2,i)=HIS(end,9);bpp(3,i)=HIS(end,6);
two(1,i)=histtime(end);two(2,i)=histg(end);two(3,i)=histf(end);
end
%%
format long
[mean(ne,2),min(ne,[],2),max(ne,[],2)]
[mean(bb,2),min(bb,[],2),max(bb,[],2)]
[mean(bpp,2),min(bpp,[],2),max(bpp,[],2)]
[mean(two,2),min(two,[],2),max(two,[],2)]



