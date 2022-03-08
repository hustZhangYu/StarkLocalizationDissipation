% һ��С�������������໥����ϵͳ���У�Liouvillan gap���������ܶȵĹ�ϵ
% ����Ϊ����ܺͷ����ӵ�pauil������ԭ���йء�
% ��һ�׶ν������ ���㲻ͬ�������� Liouvillan gap�ı仯
 % �������� ��00000-11111

 
 L=5; % ��ϵ��С
 dim=2^L;
m=['%0',num2str(L),'d'];
F=1;
%  bin2dec(x) �鿴��ʸ����ռ�����
 
 % Hamlitonian  �������غ���ӿռ�
 
 
 
 %------------------------------------------------------------------------
 % ������ܶ������ɢ���
 
 J=-1/2;
 % ���ܶ���hopping term����
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
 

 % ���ܶ������ܹ���
 for x=0:dim-1
    for i=1:L
%        H(x+1,x+1)=H(x+1,x+1)+F*(L-i+1)*tl(x,i,L);
    H(x+1,x+1)=H(x+1,x+1)+F*(i)*tl(x,i,L);
    end
 end
 [Ev,E]=eig(H,'vector')
 

 
 
%  surf(H)
 
 %  ��ɢ���
 % �������ǹ������λ�����˵ĺ�ɢ��� (����)
 C1=zeros(dim,dim);
 CL=zeros(dim,dim);
 gamma=1; %�����˺�ɢǿ��
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
% ����Liouvillian operator
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
  
 
%  ������ܶ��������Ƿ���ȷ
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
 % �ú�������һ��ʮ���Ƶ���ֵ
 % ����ֵ�� ����ֵ�����Ʊ���µ�n��λ�õ�����������������
 
 a=dec2bin(x);
 a1=str2num(a);
 m=['%0',num2str(L),'d'];
 a2=num2str(a1,m);
 kf=str2num(a2(n)); 
 
 end
