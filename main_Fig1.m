% MAIN function 1: compare 'LLL','Boosted LLL','Deep LLL','SR-SIC' for
% random bases
% author: Shanxiang Lyu, shanxianglyu@gmail.com

clc;
clear all;
linestyles = cellstr(char('-','--','-.','--','-'));
SetColors=lines(10); 
Markers=['o','x','+','*','s'];
rho_range=[0.1:0.1:0.6];
FIT=[];
ALGORITHMS=[1,2,3,4]; %indicate the four types of algorithms
algorithmbox={'LLL','Boosted LLL','Deep LLL','SR-SIC'};
metricbox={'Orthogonality defect','Basis length','Shortest Vector','Average reduction time/s'};
M=1;%M indicates the desired metric
 
 for rho=rho_range 
     
    n=20;%dimension
    %Correlated matrix PHI
    PHI=zeros(n,n);
    for k=1:n
        PHI(1,k)=rho^(k-1);%first row
    end
    PHI(1,1)=0;
    c1=PHI(1,:);
    for k=1:n %upper triangular
        PHI(k,k:n)= c1(1:n-k+1);
    end
    PHI=PHI+PHI'+eye(n);%upper+lower+eye(n)

    for monte=1:500 %the number of Monte-Carlo runs
       B=PHI*randn(n,n);%correlated n-dim random matrix
        for j=ALGORITHMS
            switch j
               case 1 %normal lll
                  expression = 'BR=LLL(B);';
                case 2 % boosted LLL
                    expression = 'BR=LLLboost(B);';
                case 3 %LLL deep
                    expression = 'BR=LLLdeep(B);';
                case 4 %SR-SIC
                    expression = 'BR=SRSIC(B);';
            end
           tic
           eval(expression);
           ttoc=toc;

           switch M
               case 1
                   fit(monte,j)=abs(prod(sqrt(abs(diag(BR'*BR))))/sqrt(abs(det(BR'*BR))));
               case 2
                   fit(monte,j)=max(diag(BR'*BR).^.5);
               case 3
                   fit(monte,j)=min(diag(BR'*BR).^.5);
               case 4
                   fit(monte,j)=ttoc;
           end

        end
    end
        
    FIT=[FIT,mean(fit,1)'];
end

figure(1)
    for j=ALGORITHMS
semilogy(rho_range,FIT(j,:),[linestyles{j} Markers(j)],'Color',SetColors(j,:),'Linewidth',1.5);
        hold on;
        grid on;
    end
hold off;
legend(algorithmbox(ALGORITHMS));
xlabel('Correlation coefficient $\rho$','Interpreter','latex');
ylabel(metricbox(M),'Interpreter','latex');

