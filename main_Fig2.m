% MAIN function 1: compare 'LLL','Boosted LLL','Deep LLL','SR-SIC' for
% bases in Integer-Forcing
% author: Shanxiang Lyu, shanxianglyu@gmail.com

clc;
clear all;
linestyles = cellstr(char('-','--','-.','--','-'));
SetColors=lines(10); 
Markers=['o','x','+','*','s'];
SNR_range=[0:1:10];
FIT=[];
ALGORITHMS=[1,2,3,4]; %indicate the four types of algorithms
algorithmbox={'LLL','Boosted LLL','Deep LLL','SR-SIC'};
metricbox={'Ergodic rate $R_E$/bpcu','Basis length','Shortest Vector','Average reduction time/s'};
M=2;%M indicates the desired metric

 for SNR_val=SNR_range 
     
    n=20;%dimension
    P=10^(SNR_val/10);%the SNR value
    
        for monte=1:500 %the number of Monte-Carlo runs
        
           H=randn(n,n);
           [U,S,V]=svd(H'*H+(1/P)*eye(n));
           B=S^(-0.5)*V';%basis in integer-forcing
              
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
                       fit(monte,j)=n*max(0,.5*log2(P/max(diag(BR'*BR))));
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
semilogy(SNR_range,FIT(j,:),[linestyles{j} Markers(j)],'Color',SetColors(j,:),'Linewidth',1.5);
        hold on;
        grid on;
    end
hold off;
legend(algorithmbox(ALGORITHMS));
xlabel('SNR/dB','Interpreter','latex');
ylabel(metricbox(M),'Interpreter','latex');

