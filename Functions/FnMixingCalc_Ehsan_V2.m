function [Ahat,MinEV]=FnMixingCalc_Ehsan_V2(w,k,c,n)
g=nchoosek(n-1,k-1);
% Cx=nchoosek(1:c,g);%Slow
Cx=VChooseK(1:c,g);%Fast (c code)
% AllPerms=nchoosek(1:size(w,2),size(w,1));
Cg=size(Cx,1);
% f=nchoosek(n,g);
% Cg=nchoosek(c,g);
for i=1:Cg
    if numel(size(w))>=3
        wbar=w(:,:,Cx(i,:));
        wbar=reshape(wbar, [size(wbar,1),size(wbar,2)*size(wbar,3)]);
    else
        wbar=w(:,Cx(i,:));
    end
    R=wbar*wbar';
    [V,D]= eig(R);
    [Val,Ix] =sort(diag(D));
    MinEV(i)=abs(Val(1))/abs(Val(end));
    MinEV(i)=abs(Val(1));
    MinInd(i)=Ix(1);
    VV(:,i)=V(:,Ix(1));
end
[Val,Ix] =sort( MinEV);
f=n;
% temp=VV(:,Ix(1:Cg));
Ahat=VV(:,Ix(1:f));
