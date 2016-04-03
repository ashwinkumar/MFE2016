%Empirical Method HW2 Problem 2

clear all
clc

%Problem 1 part 1.1
filename ='datahwk2_problem1.xlsx';
data = xlsread(filename);
[m,n]=size(data);

meanReturn = mean(data(:,:));
meanReturn = transpose(meanReturn)

%part 1.2
bsMeanSTD = NaN(n,1);
% 10000 bootstraps samples
nboot=10000;
for i=1:n
bootstat=bootstrp(nboot,@mean,data(:,i));
bsMeanSTD(i)=std(bootstat);
end;
bsMeanSTD

fundNames={'CSTCVAH' 'CSTEMNH' 'CSTEVDH' 'CSTDISH' 'CSTRARH' 'CSTFIAH' 'CSTGLMH' 'CSTMNFH'};

T=table(meanReturn, bsMeanSTD,'RowName', fundNames)

%part 1.3 & 1.4

autocorr_data = zeros(12+1,n);
figure
for iacf =1:n
   subplot(2,4,iacf)
   autocorr(data(:,iacf),12,[],2)
   autocorr_data(:,iacf) = autocorr(data(:,iacf),12,[],2);
   % Ljung-Box test 
   [h(iacf),pValue(iacf),stat(iacf),cValue(iacf)] = lbqtest(data(:,iacf),'lags',12);
end
mat2dataset(autocorr_data,'VarNames',{'CSTCVAH', 'CSTEMNH', 'CSTEVDH', 'CSTDISH', 'CSTRARH', 'CSTFIAH', 'CSTGLMH', 'CSTMNFH'},'ObsNames',{'lag 0','lag 1', 'lag 2', 'lag 3', 'lag 4', 'lag 5','lag 6', 'lag 7', 'lag 8', 'lag 9', 'lag 10','lag 11', 'lag 12'})
decision=transpose(h);
Pval=transpose(pValue);
table(decision,Pval,'RowName', fundNames)
%if ture, data are autocorrelated, if falsee, not autocorrelated

%part 1.5
VWRETD = xlsread('CSRP_VW_problem1.xlsx');
figure
autocorr(VWRETD(:,2),12,[],2)

%part 1.7
% Vector block bootstraps with a block size of 4
blockSize=4;
for j=1:n
bsData = block_bootstrap(data(:,j), nboot, blockSize);
blockbsSTD(j)=std(mean(bsData));
end
blockbsSTD=transpose(blockbsSTD);
T=table(meanReturn, bsMeanSTD,blockbsSTD,'RowName', fundNames)


%problem 2
%part 2.1
factorsData=xlsread('datahwk2_problem2.xlsx');
factorsData=factorsData(85:286,2:7);
factorsData(:,4)=[];
[nrow ncol]=size(factorsData)
factors=zeros(nrow,ncol);
factors(:,1:2)=factorsData(:,4:5);
factors(:,3:5)=factorsData(:,1:3);
data2=data/100;
b_hat = zeros(5,n);
Aeq=ones(1,5);
for j=1:n
b_hat(:,j)=lsqlin(factors,data2(:,j),[],[],Aeq,1);
R_star(:,j)=factors(:,:)*b_hat(:,j);
end
ds=mat2dataset(b_hat,'VarNames',{'CSTCVAH', 'CSTEMNH', 'CSTEVDH', 'CSTDISH', 'CSTRARH', 'CSTFIAH', 'CSTGLMH', 'CSTMNFH'},'ObsNames',{'beta_1', 'beta_2', 'beta_3', 'beta_4', 'beta_5'})
%part 2.2
meanR_star = mean(R_star(:,:));
meanR_star = transpose(meanR_star);
meanR_hat=zeros(n,1);
gamma=zeros(n,1);
for i=1:n
gamma(i,1)=std(data2(:,i))/std(R_star(:,i));
R_hat(:,i)=R_star(:,i)*gamma(i,1);
meanR_hat(i,1)=mean(R_hat(:,i));
end
meanR_hat;
meanR=meanReturn/100,
table(meanR,meanR_star,meanR_hat)
%part 2.3
autocorr_clonedata = zeros(12+1,n);
figure
for ihat =1:n
subplot(2,4,ihat)
autocorr(R_hat(:,ihat),12,[],2)
autocorr_clonedata(:,ihat) = autocorr(R_hat(:,ihat),12,[],2);
% Ljung-Box test 
[h_hat(ihat),pValue_hat(ihat),stat_hat(ihat),cValue_hat(ihat)] = lbqtest(R_hat(:,ihat),'lags',12);
end
mat2dataset(autocorr_clonedata,'VarNames',{'CSTCVAH', 'CSTEMNH', 'CSTEVDH', 'CSTDISH', 'CSTRARH', 'CSTFIAH', 'CSTGLMH', 'CSTMNFH'},'ObsNames',{'lag 0','lag 1', 'lag 2', 'lag 3', 'lag 4', 'lag 5','lag 6', 'lag 7', 'lag 8', 'lag 9', 'lag 10','lag 11', 'lag 12'})
decision_hat=transpose(h_hat);
Pval_hat=transpose(pValue_hat);
table(decision_hat,Pval_hat,'RowName', fundNames)
