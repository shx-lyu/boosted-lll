function [B] = LLLdeep(B)
% function: deep LLL algorithm for lattice reduction
% input: lattice basis B
% output: reduced basis B
% author: Shanxiang Lyu, shanxianglyu@gmail.com
% ref: S-E¡°Lattice basis reduction: Improved practical algorithms and solving subset sum problems,¡± Math. Program., vol. 66, pp. 181¨C199, 1994.

[Q,R]=qr(B);
[~,n]=size(B);
delta=0.99;
k=2;
flag=0;
while k<=n

    for l=k-1:-1:1  %size reduction
     q=round(R(l,k)/R(l,l));
          if q~=0 
              R(:,k)=R(:,k)-q*R(:,l);
          end    
    end

    B=Q*R;%get the current size reduced version
    C=norm(B(:,k))^2;i=1;
    while i < k
     if C>delta*abs(R(i,i))^2
         C=C-R(i,k)^2;
         i=i+1;
         flag=0;
     else %otherwise, swap
         tmp=B(:,k);
         B(:,k)=[];
         if i>=2
             B0=B(:,1:i-1);
             B1=B(:,i:end);
             B=[B0,tmp,B1];
         else
             B1=B(:,i:end);
             B=[tmp,B1];
         end
         %updage QR
         [Q,R]=qr(B);
         k=max(i,2);
         flag=1;
         break;%goto outer-size-reduction
     end
    end

    if flag==0
    k=k+1;
    end
end
B=Q*R;
end