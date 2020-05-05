function [B] = LLL(B)
% function: LLL algorithm for lattice reduction
% input: lattice basis B
% output: reduced basis B
% author: Shanxiang Lyu, shanxianglyu@gmail.com
% ref: Lenstra, A. K.; Lenstra, H. W.; and Lovasz, L. "Factoring Polynomials with Rational Coefficients." Math. Ann. 261, 515-534, 1982.

[Q,R]=qr(B);
[m,n]=size(B);
delta=0.99;
i=2;
V=[0 1;1 0];   
while i<=n

     for k=i-1:-1:1 
         q=round(R(k,i)/R(k,k));
              if q~=0 
                  R(:,i)=R(:,i)-q*R(:,k);
              end    
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