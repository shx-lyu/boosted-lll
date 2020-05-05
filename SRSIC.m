function [B]=SRSIC(B)
% function: SR-SIC algorithm for lattice reduction
% input: lattice basis B
% output: reduced basis B
% author: Shanxiang Lyu, shanxianglyu@gmail.com

n=size(B,2);
eta=0;
while eta<1
    [~,IND]=sort(diag(B'*B));
    B=B(:,IND);%basis vectors are sorted from short to long
   %------an SIC algorithm----%
    x=zeros(n-1,1);
    [Q,R]=qr(B);
    for k=n-1:-1:1 
     x(k)=round(R(k,n)/R(k,k));
          if x(k)~=0 
              R(:,n)=R(:,n)-x(k)*R(:,k);
          end    
    end

    if norm(R(:,n))^2<norm(B(:,n))^2
       B(:,n)=Q*R(:,n);
       eta=0;
    else
       eta=eta+1;
    end
end
end