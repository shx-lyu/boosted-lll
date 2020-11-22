function [B,T,nloop] = LLLboost(B,num_routes)
% function: boosted LLL algorithm (num_routes=1,3,9)
% input: lattice basis B
% output: reduced basis B
% author: Shanxiang Lyu, shanxianglyu@gmail.com
% ref:Boosted KZ and LLL Algorithms. IEEE Trans. Signal Process. 65(18): 4784-4796 (2017)

if nargin == 1
    num_routes =1; 
end

[Q,R]=qr(B);

[m,n]=size(B);
T=eye(n);
delta=0.99;%0.99
nloop=0;
i=2;
while i<=n

   [r1,t1,r2,t2]=PSIC(R(:,1:i),T(:,1:i),num_routes);%num_routes=1,3,9
   nloop=nloop+(i);
   
    R(:,i)=r1;T(:,i)=t1;%PSIC reduction
    
    cdelta=round(R(i-1,i)/R(i-1,i-1));
    if delta*norm(R(i-1,i-1))^2>abs(R(i,i))^2+abs(R(i-1,i)-cdelta*R(i-1,i-1))^2  %Lovasz fails
        
        R(:,i)=r2;T(:,i)=t2;%the shortest one that has <1/2 on one
        %return to a valid value and swap
        V=[0 1;1 0];   
        u=V(:,1);
          %Givens matrix
           tempsum=sqrt((R(i-1,i-1)*u(1)+R(i-1,i)*u(2))^2+(R(i,i)*u(2))^2);
           alpha=(R(i-1,i-1)*u(1)+R(i-1,i)*u(2))/tempsum;
           beta=-(R(i,i)*u(2))/tempsum;
           
           R(:,i-1:i)=R(:,i-1:i)*V;%SWAP
           T(:,i-1:i)=T(:,i-1:i)*V;
           G=eye(m);
           G(i-1:i,i-1:i)=[alpha,-beta;beta,alpha];
  
           R=G*R;%RESTORE
           Q=Q*(G)';
           i=max(i-1,2);
    else
        i=i+1;
    end
    
    
end
B=Q*R;
end

function [r1,t1,r2,t2]=PSIC(R,T,num_routes)
%size reduction one column; 1 is non-greedy and second is greedy
% r: a column of the R matrix;
% t: corresponding Unitary transform matrix
if nargin == 2
    num_routes=1;
end

%use the flag for the PSIC mod; 1 3 9
if num_routes==1
    M=0;
elseif num_routes==3
    M=[0,-1,1];
else
    M=[0 0 0 -1 -1 -1 1 1 1;
       0 -1 1 0 -1 1 0 -1 1];
end
flag2=size(M,1);

%estab coe matrix and reduced matrix
[m,n]=size(R);
mt=size(T,1);
%m=n;%added in Feb20-2017, to settle over complete

Rbig=zeros(m,num_routes+1);
Tbig=zeros(mt,num_routes+1);%modify from m to n, Feb24
val=zeros(1,num_routes+1);
ind=zeros(1,num_routes+1);


%Push the original inside, at the last position
Rbig(:,num_routes+1)=R(:,n);
Tbig(:,num_routes+1)=T(:,n);
val(num_routes+1)=norm(R(:,n));
ind(num_routes+1)=round(R(n-1,n)/R(n-1,n-1));

r0=R(:,n);t0=T(:,n);
for iter=1:num_routes
      for k=n-1:-1:1 
          q=round(R(k,n)/R(k,k));
          if n-k<=flag2
              q=q+M(n-k,iter);%shited value based on choices
          end
          if q~=0  %do a unimodular for B, BS is not changed
              R(:,n)=R(:,n)-q*R(:,k);
              T(:,n)=T(:,n)-q*T(:,k);
          end    
      end
  Rbig(:,iter)=R(:,n);%store the new values in a table
  Tbig(:,iter)=T(:,n);
  val(iter)=norm(R(:,n));
  ind(iter)=round(R(n-1,n)/R(n-1,n-1));
  %RESTORE THE ORIGINAL VALUE
  R(:,n)=r0;T(:,n)=t0;
end

%SHORTEST SET
pos1=find(val==min(val));
r1=Rbig(:,pos1(1));
t1=Tbig(:,pos1(1));

indSet=find(ind==0);%find the one that has 1/2 in the first position
R3=Rbig(:,indSet);
T3=Tbig(:,indSet);
val3=val(indSet);

pos3=find(val3==min(val3));
if  numel(pos3)==0 %this cannot happen
    r2=r1;
    t2=t1;
else
    r2=R3(:,pos3(1));
    t2=T3(:,pos3(1));
end
end
