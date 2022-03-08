clear all;
clc;

L=5;J=-1/2;F=0.1;gamma=1;

data=zeros(1,L);
for d=1:L
    
    Na=d;

    H=0;
    for N=1:Na
        H1=Ha(L,N,J,F);
        H=blkdiag(H1,H);
    end


    Lk1=La(L,Na,gamma);
    Lk2=Lb(L,Na,gamma);

    dim=size(H,1);

    I=eye(dim);
    L1=-1i*(kron(I,H)-kron(H.',I));

    Lk=Lk1;
    L2=kron(conj(Lk),Lk)-(kron(I,Lk'*Lk))/2-kron((Lk'*Lk).',I)/2;
    L1=L1+L2;

%     Lk=Lk2;
%     L2=kron(conj(Lk),Lk)-(kron(I,Lk'*Lk))/2-kron((Lk'*Lk).',I)/2;
%     L1=L1+L2;

    U=expm(L1*2*pi);
    [Ev,E]=eig(U,'vector');

    plot(E,'-*')
    hold on;

    E1=sort(abs(E));
    k=E1(length(E1)-1);
    k1=-log(k)/(2*pi);
    data(1,d)=k1;
end
figure()
plot(linspace(1,L,L)/L,data,'*')

%blkdiag 直和函数

function H=Ha(L,N,J,F)
% 这个函数的主要目的是给出固定粒子数和固定大小的哈密顿量
% 要得到特定的Liouvillian 算符，还需要对全粒子空间做直和

key=nchoosek(1:L,N);
basis=zeros(size(key,1),L);
for i=1:size(key,1)
    for j=1:N
        basis(i,key(i,j))=1;     
    end
end

% 接下来构造哈密顿量
dim=size(key,1);
H=zeros(dim,dim);


% hopping term

for k=1:L-1
    bz=zeros(1,L);
    bz(1,k)=1;bz(1,k+1)=-1;            % 这里主要思路是，前后两组基矢相减 ,（10）-（01）会出现（1，-1）其余都是0
    for i=1:dim-1
        for j=i+1:size(key,1)
            if basis(i,:)-basis(j,:)==bz
               H(i,j)=J;
               H(j,i)=J;
            end
        end
    end
end

% potential

V=linspace(1,L,L);
for i=1:dim
   H(i,i)=F*sum(basis(i,:).*V); 
end


end

function Lk=La(L,N,gamma)
% 不同尺寸下面的边界左耗散算符
% N 是最大粒子数

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
bz(1,1)=1;

for i=1:dim-1
    for j=i+1:dim
        if i==21 && j==26
            basis(i,:)-basis(j,:)
        end
        if basis(i,:)-basis(j,:)==bz
           Lk(j,i)=gamma;
        end
    end
end

end

function Lk=Lb(L,N,gamma)
% 不同尺寸下面的边界左耗散算符
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
bz(1,L)=1;

for i=1:dim-1
    for j=i+1:dim
        if i==21 && j==26
            basis(i,:)-basis(j,:)
        end
        if basis(i,:)-basis(j,:)==bz
           Lk(j,i)=gamma;
        end
    end
end

end

