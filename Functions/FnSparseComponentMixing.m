% Sparse Componnet Mixing Function: FnSparseComponentMixing
%% Inputs:
% m: The number of sensors a.k.a the ambient dimension
% n: The number of sources
% k: The number of acitive source in each time point (nonzero elements
% in each col of S) or the dimension of each subspace
% N: The cardinal of each subspace in kSCA mode
% Nk: the cardinal of  each subspace in MSCA mode
% MixingMode: Sparse Component Mixing Mode:'kSCA'=> Uniform Sparse Component
% Mixing,'MSCA'=> Multiple Sparse Component Mixing, 'OnOffGauss' (Feb 2020
% was added by EE)
%% Outputs:
% X: The mixtures / mixed data/observed data/data points : m*T
% S: The source matrix/The sparse componenets/The weights in linaer
% combination of the columns of A a basis of the subspaces:n*T
% A: The mixing matrix/ the basis of subspaces/ The dictionary: m*n
% OrthA: The orthonormalized A: each column of OrthA is obtaind of
% orth(A(:,P(i,:)) which P(i,:) is the inices of selected column of A:
% Labels: ground-Truth for subspace clustering:1*T,Label(i) is a number
% between 1 and the total number of the subspaces
% SubspaceInds: each row of P shows the indices of the columns of A that
% has generated each subspace.P is c*k in kSCA mode which c is the total
% number of subspace and is a cell in MSCA mode.
% SubspaceNum: The total number of the subspaces.
% *************************************************************************
% Ehsan Eqlimi, 26 Bahman 1394
% Edit 1:Add a parameter to control the variance of gaussian noise over the
% inactive(nonzeros) sources: 9 Xordad 1394 and Add PermkSCA and kSCANoisy
%Edit2:  Ksca ans KscaNoisy mode were modified using
%(randn(n-k,1)-->randn(n-k,N(i))). OneOffGauss mode was added (based on
%Movahedi's paper).
% *************************************************************************
function [X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Var,IterNum,RankTh,AMode,MixingMode)

if AMode
    BestKa=0;
    Iter=0;
    while BestKa<RankTh && Iter<IterNum
        A=FnColNormalizer(randn(m,n));
        [R,Ka]=k_rank_EE(A);
        if Ka>BestKa
            BestKa=Ka;
            BestA=A;
        end
        Iter=Iter+1;
    end
    A=BestA;
else
    A=FnColNormalizer(randn(m,n));
end
% A=randn(m,n);
% A=FnColNormalizer(A);
P=nchoosek(1:n,k);
c=length(P);
OrthA=[];
% S=[];
X=[];
Labels=[];
S=[];
% % % A=[0.5525 0.3919 0.5707 0.3934 0.6904;0.6863 -0.6066 0.5166 0.8634
% -0.6007;0.4730 0.6917 -0.6383 -0.3158 0.4032];%Whashizava example
switch MixingMode
    case 'kSCA' %
        SubspaceInds=P;
        for i=1:c
            %             a=orth(A(:,P(i,:)));
            a=(A(:,P(i,:)));
            OrthA=[OrthA a];
            s=randn(k,N(i));
            sWithZeros=zeros(n,N(i));
            sWithZeros(P(i,:),:)=s;
            IdxActive = zeros(n,1);
            IdxActive(P(i,:))=1;
            IdxActive=logical(IdxActive);
            sWithZeros(~IdxActive,:)=Var*randn(n-k,N(i));
            S=[S sWithZeros];
            x=a*s;
            X=[X x];
            %             S=[S s];
            Labels=[Labels i*ones(1,N(i))];
        end
        
    case 'PermkSCA'
        SubspaceInds=P;
        S=zeros(n,T);
        for i=1:T
            IdxActive = zeros(n,1);
            %         IdxRandPerm = randperm(k);
            Labels(i) = randperm(size(P,1),1);
            IdxRandPerm=P(Labels(i),:);
            IdxActive(IdxRandPerm) = 1;
            IdxActive = logical(IdxActive);
            S(IdxActive,i) = randn(k,1);
            
            S(~IdxActive,i) = Var*randn(n-k,1);
        end
     
        
        X=A*S;
        
    case 'kSCANoisy'
        SubspaceInds=P;
        for i=1:c
            %             a=orth(A(:,P(i,:)));
            % %             a=(A(:,P(i,:)));
            % %             OrthA=[OrthA a];
            s=randn(k,N(i));
            sWithZeros=zeros(n,N(i));
            sWithZeros(P(i,:),:)=s;
            IdxActive = zeros(n,1);
            IdxActive(P(i,:))=1;
            IdxActive=logical(IdxActive);
            sWithZeros(~IdxActive,:)=Var*randn(n-k,N(i));
            S=[S sWithZeros];
            %             x=a*s;
            %             X=[X x];
            %             S=[S s];
            Labels=[Labels i*ones(1,N(i))];
        end
        %         for j=1:n
        %             S(j,:)=S(j,:)-mean(S(j,:));
        %         end
        
        X=A*S;
        
        
    case 'MSCA'
        Labels=0;
        for j=1:k
            Added(j)=max(Labels);
            kj=j;
            Pj=nchoosek(1:n,kj);
            SubspaceInds{j}=Pj;
            cj=length(Pj);
            for i=1:cj
                a=orth(A(:,Pj(i,:)));
                OrthA=[OrthA a];
                s=randn(kj,Nk(i,j));
                sWithZeros=zeros(n,Nk(i,j));
                sWithZeros(Pj(i,:),:)=s;
                
                S=[S sWithZeros];
                x=a*s;
                X=[X x];
                %                 S=[S s];
                Labels=[Labels (i+Added(j))*ones(1,Nk(i,j))];
                
            end
            
        end
        Labels=Labels(2:end);
        
    case 'OnOffGauss' %Equation 8 of Movahedi's paper
        %         %generating the sources
        %         p=k/n;
        %         z=rand(n,T);
        %         S=(z>p).*randn(n,T)*std_noise+(z<=p).*(randn(n,T)*std_source);
        %         for j=1:n
        %             S(j,:)=S(j,:)-mean(S(j,:));
        %         end
        
        %         mu=[0 ;0];
        %         sigma=[Var];
        %         p=[k/n;1-(k/n)];
        %         gm = gmdistribution(mu,sigma,p)
        %         Pr=k/n;
        %         Sigma_off=Var;
        %         idxNZ = zeros(n,1);
        %         rndPerm = randperm(n);
        %         idxNZ(rndPerm(1:k)) = 1;
        %         idxNZ = logical(idxNZ);
        %
        %         s = zeros(n,1);
        %         % if nargin < 4
        %         s(idxNZ) = randn(k,1);
        %         % else
        %         %     s(idxNZ) = fixedActiveValue;
        %         % end
        %         s(~idxNZ) = Var*randn(n-k,1);
end
SubspaceNum=max(Labels);