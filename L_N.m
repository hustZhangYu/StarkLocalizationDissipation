% 一个小插曲，关于无相互作用系统当中，Liouvillan gap随粒子数密度的关系
% 我认为这可能和费米子的pauil不相容原理有关。
% 第一阶段进攻命令： 计算不同粒子数下 Liouvillan gap的变化
 % 二进制数 从00000-11111

 
 L=5; % 体系大小
 dim=2^L;
m=['%0',num2str(L),'d'];
F=1;
%  bin2dec(x) 查看基矢粒子占据情况
 
 % Hamlitonian  粒子数守恒的子空间
 
 
 
 %------------------------------------------------------------------------
 % 构造哈密顿量与耗散算符
 
 J=-1/2;
 % 哈密顿量hopping term构造
 H=zeros(dim,dim);
 
 for x=0:dim-2
     for i=1:L-1
         if tl(x,i,L)==1 && tl(x,i+1,L)==0

            a=dec2bin(x);
            a1=str2num(a);
            m=['%0',num2str(L),'d'];
            a=num2str(a1,m);
            
            a(i)='0';a(i+1)='1';
            x1=bin2dec(a);
            H(x+1,x1+1)=J; 
         end
     end     
 end
 
 H=H+H';
 

 % 哈密顿量势能构造
 for x=0:dim-1
    for i=1:L
%        H(x+1,x+1)=H(x+1,x+1)+F*(L-i+1)*tl(x,i,L);
    H(x+1,x+1)=H(x+1,x+1)+F*(i)*tl(x,i,L);
    end
 end
 [Ev,E]=eig(H,'vector')
 

 
 
%  surf(H)
 
 %  耗散算符
 % 这里我们构造的是位于两端的耗散算符 (湮灭)
 C1=zeros(dim,dim);
 CL=zeros(dim,dim);
 gamma=1; %代表了耗散强度
 for x=0:dim-1
     if tl(x,1,L)==1 
        a=dec2bin(x);
        a1=str2num(a);
        m=['%0',num2str(L),'d'];
        a=num2str(a1,m);
        a(1)='0';
        x1=bin2dec(a);
        C1(x+1,x1+1)=gamma; 
     end   
 end
 
  for x=0:dim-1
     if tl(x,L,L)==1 
        a=dec2bin(x);
        a1=str2num(a);
        m=['%0',num2str(L),'d'];
        a=num2str(a1,m);
        a(L)='0';
        x1=bin2dec(a);
        CL(x+1,x1+1)=gamma; 
     end   
  end
 
%--------------------------------------------------------------------------- 
% 构造Liouvillian operator
  I=eye(dim);
  L1=-1i*(kron(I,H)-kron(H.',I));
  
  Lk=C1;
  L2=kron(conj(Lk),Lk)-(kron(I,Lk'*Lk))/2-kron((Lk'*Lk).',I)/2;
  L1=L1+L2;
  
   Lk=CL;
  L2=kron(conj(Lk),Lk)-(kron(I,Lk'*Lk))/2-kron((Lk'*Lk).',I)/2;
   L1=L1+L2;
  
  U=expm(L1*2*pi);
  [Ev,E]=eig(U,'vector');
  
  plot(E,'.')
  
  E1=sort(abs(E));
  k=E1(length(E1)-1);
  k1=-log(k)/(2*pi)
  
 
%  检验哈密顿量构造是否正确
%   [Ev1,E1]=eig(H,'vector');
%   figure()
%   plot(E1)
% mesh(Ev1.*conj(Ev1))
  
%  figure()
%  subplot(1,2,1)
%  surf(C1)
%  subplot(1,2,2)
%  surf(CL)
 

 function kf=tl(x,n,L)
 % 该函数给定一个十进制的数值
 % 返回值是 该数值二进制表达下第n个位置的数（从左往右数）
 
 a=dec2bin(x);
 a1=str2num(a);
 m=['%0',num2str(L),'d'];
 a2=num2str(a1,m);
 kf=str2num(a2(n)); 
 
 end
