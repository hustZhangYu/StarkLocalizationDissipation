function Lk = BoundaryLossOperator(L,N,gamma,site)
%BOUNDARYLOSSOPERATOR Historical boundary-loss matrix in descending sectors.
% Basis order is the N-particle sector, then N-1, ..., 1, then vacuum.
% This preserves the occupation-basis construction used in L_N1.m.
validateattributes(L, {'numeric'}, {'scalar','integer','positive'});
validateattributes(N, {'numeric'}, {'scalar','integer','nonnegative','<=',L});
validateattributes(site, {'numeric'}, {'scalar','integer','positive','<=',L});

basis=[];
for k1=1:N
    N1=N-k1+1;
    key=nchoosek(1:L,N1);
    basis1=zeros(size(key,1),L);
    for i=1:size(key,1)
        for j=1:N1
            basis1(i,key(i,j))=1;
        end
    end
    basis=[basis;basis1];
end

basis=[basis;zeros(1,L)];
dim=size(basis,1);
Lk=zeros(dim,dim);
bz=zeros(1,L);
bz(1,site)=1;

for i=1:dim-1
    for j=i+1:dim
        if all(basis(i,:)-basis(j,:)==bz)
            Lk(j,i)=gamma;
        end
    end
end
end

