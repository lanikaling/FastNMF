% Yale 64
%The face image can be displayed in matlab with the following command lines:
%===========================================
%%
load('Yale_64x64.mat');
% faceW = 64; 
% faceH = 64; 
% numPerLine = 11; 
% ShowLine = 2; 
% 
% Y = zeros(faceH*ShowLine,faceW*numPerLine); 
% for i=0:ShowLine-1 
%   	for j=0:numPerLine-1 
%     	Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(fea(i*numPerLine+j+1,:),[faceH,faceW]); 
%   	end 
% end 
% 
% imagesc(Y);colormap(gray);
%===========================================
%%
M=fea(1:44,:)';
svd(M)';

% seed=rng;
rng(seed);
k=3;
X0=rand(4096,k);Y0=rand(k,size(M,2));
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
