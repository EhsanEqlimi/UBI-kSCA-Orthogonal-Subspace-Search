clc;clear;close all;

m=64;
n=4;
k=2;
%% Test On randn sources
T= 2000;%c*30; %1000; % The Number of data poins a.k.a the number of samples (time points)
% zero sources (non-strict saprsity)
AMode=0; % k-SCA Condition for A is satisfied
Iter_Num_A=1000;% The number of iteration to generate a good mixing matrix
Rank_Th_A=0.1; % to generate a good mixing matrix
DualMode=0; %Twice subpspace clusteing
BAS_Set=[];NMSE_Set=[];Norm2Error_Set=[];Sub_Error_Set=[];Error_Set=[];
%% SCA Mixing (Simulted Sources)
%%'kSCA';%'kSCA';% Uniform Sparse Component Mixing  %'MSCA'=>Multiple Sparse Component Mixing, 'PermKSCA'
MixingMode='kSCA';% {'kSCANoisy'%'PermkSCA'%kSCA';%'MSCA';%OnOffGauss};
Orth=0; %if Orth=1 A is orth
c=nchoosek(n,k); % The number of all possible r-dim subspaces
Actualc=c;
N=zeros(1,c);
N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mode
for j=1:k
    for i=1:nchoosek(n,j)
        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% N showes the number of each subspace in MSCA mode
    end
end
Sigma=0;
[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Sigma,Iter_Num_A,Rank_Th_A,AMode,MixingMode);
ColorLab={'r','b','g','m'};

figure,
for i=1:n
%     Temp=[zeros(1,(i-1)*L*OL) Cosine zeros(1,n*L-L-(i-1)*L*OL)] ;%figure,plot(S);
%     S(i,:)=Temp(:,1:RealL);
    subplot(n,1,i)
    plot(1:size(S,2),S(i,:),ColorLab{i})
end
GroundTruthLables=Labels;
alpha=10;
r = 0; affine = false; outlier = true; rho = 1.5;
[missrate,C,grps,Clus] = SSC(X,r,affine,alpha,outlier,rho,GroundTruthLables);
Sub=FnSubspaceCalc(X,grps',1,Actualc);
Subb=squeeze(Sub)
figure,stem(grps)
title('SCA-ADMM-Alpha-10')

figure,stem(GroundTruthLables)
title('GroundTruthLables')
[Sub,OrthSub]=FnSubspaceCalc(X,grps',k,n);
Subb=squeeze(Sub);
% Ahat=FnMixingCalc(OrthSub,2,6,4)
[Ahat,MinEV]=FnMixingCalc_Ehsan(OrthSub,k,c,n)
[MixingIdentificatioError, MixingVectorerror,NMSE,NMSSum,AhatNew,Norm2Error,Error] = FnMixingIdentificationError(A,Ahat);
MixingVectorerror

%%K-means
 %%K-means
 [LabelsKmenas,KmeansClusters]=kmeans(X',Actualc);
 figure,stem(LabelsKmenas)
 title('k-means')
 Missrate_Kmeans = Misclassification(LabelsKmenas,GroundTruthLables)