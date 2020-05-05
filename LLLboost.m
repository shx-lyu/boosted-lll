function [B] = LLLboost(B)
% function: boosted LLL algorithm (L=1)
% input: lattice basis B
% output: reduced basis B
% author: Shanxiang Lyu, shanxianglyu@gmail.com
% ref:Boosted KZ and LLL Algorithms. IEEE Trans. Signal Process. 65(18): 4784-4796 (2017)


[Q,R]=qr(B);
[m,n]=size(B);
delta=0.99;
i=2;
V=[0 1;1 0];   
while i<=n

     Ri=R(:,i);
     for k=i-1:-1:1 
         q=round(R(k,i)/R(k,k));
              if q~=0  
                  R(:,i)=R(:,i)-q*R(:,k);
              end    
     end
     
    if norm(R(:,i))>norm(Ri) && abs(Ri(i-1))<.5 %if the non-reduced Ri is shorter && effective
        R(:,i)=Ri;
    end
     
    if delta*norm(R(i-1,i-1))^2>abs(R(i,i))^2+abs(R(i-1,i))^2  %Lovasz fails
        %Givens rotation by G0
           tempsum=sqrt(R(i-1,i)^2+R(i,i)^2);
           alpha=R(i-1,i)/tempsum;
           beta=R(i,i)/tempsum;
           R(:,i-1:i)=R(:,i-1:i)*V;%SWAP
           G0=[alpha,beta;-beta,alpha];
           %Restore
           R(i-1:i,1:n)=G0*R(i-1:i,1:n);
           Q(1:m,i-1:i)=Q(1:m,i-1:i)*G0';

           i=max(i-1,2);
    else
        i=i+1;
    end   
end
B=Q*R;
end