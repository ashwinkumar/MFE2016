%% Empirical Method HW5

clear; clc; close all;
% this part is subject to changes
addpath('/Users/akumar/Documents/MATLAB/MFEToolbox/timeseries')
addpath('/Users/akumar/Documents/MATLAB/MFEToolbox/crosssection')
addpath('/Users/akumar/Documents/MATLAB/MFEToolbox/utility')
addpath('/Users/akumar/Acads/MFE 2016/2nd Quarter/Empirical Methods in Finance/HW/HW5')
%%  question 1
prices = xlsread('Fama_bond_prices.xlsx');
% get annually price based on 1 dollar face value
yrprices = prices(1:length(prices(:,1)),2:6)/100;
% generate matrix for log price based on 1 dollar face value
logprices = log(yrprices);
% generate matrix for yields
yields = zeros(length(logprices(:,1)),5);
for n=1:5
    yields(:,n)= -logprices(:,n)/n;
end 
% generate matrix for forward rate 
fwd = zeros(length(logprices(:,1)),4);
for n=1:4
    fwd(:,n)= logprices(:,n)-logprices(:,n+1);
end 
% generate matrix for dependent variables for regression 1
yieldchanges=NaN(length(logprices(:,1)),4);
for n=1:4
    yieldchanges(1:end-n*12,n)= yields(1+n*12:end,1)-yields(1:end-n*12,1);
end 
% gernerate matrix for independent variables
fwdspread=NaN(length(logprices(:,1)),4);
for n=1:4
    fwdspread(:,n)=fwd(:,n)-yields(:,1);
end 

% figure(1)
% plot(fwdspread)
% 
% figure(2)
% plot(yieldchanges)

% compute OLS standard error
B = NaN(4,2);
TSTAT = NaN(4,2);
S2 = NaN(4,1);
SE = NaN(4,2);
SEWHITE = NaN(4,2);
VCV= NaN(2,2);
VCVWHITE = NaN(2,2);
for i=1:4
   [B(i,:),TSTAT(i,:),S2(i,1),VCV,VCVWHITE] = ols(yieldchanges(1:end-i*12,i),fwdspread(1:end-i*12,i));
   SE(i,:)=sqrt(diag(VCV));
   SEWHITE(i,:)=sqrt(diag(VCVWHITE));
end

OLS = table(B(:,2),SE(:,2),TSTAT(:,2));
OLS.Properties.VariableNames = {'Beta' 'OLS_SE' 'T_stat' };
disp(OLS);

OLSWHITE = table(B(:,2),SEWHITE(:,2),TSTAT(:,2));
OLSWHITE.Properties.VariableNames = {'Beta' 'OLS_SEWHITE' 'T_stat' };
disp(OLSWHITE);

% compute HAC standard erros
B_hac = NaN(4,2);
TSTAT_hac = NaN(4,2);
S2_hac = NaN(4,1);
SE_hac = NaN(4,2);
VCVNW_hac = NaN(2,2);

for i=1:4
   [B_hac(i,:),TSTAT_hac(i,:),S2_hac(i,1),VCVNW_hac] = olsnw(yieldchanges(1:end-i*12,i),fwdspread(1:end-i*12,i));
   SE_hac(i,:)=sqrt(diag(VCVNW_hac));
end

HAC = table(B_hac(:,2),SE_hac(:,2),TSTAT_hac(:,2));
HAC.Properties.VariableNames = {'Beta_hac' 'Hac_SE' 'T_stat' };
disp(HAC);

%% question 2
% generate matrix for holding period returns
hpr = NaN(length(logprices(:,1)),4);
for n=1:4
    hpr(1:end-12,n)= logprices(13:end,n)-logprices(1:end-12,n+1);
end 
% generate matrix for dependent variables for regression 2
hprchanges= NaN(length(logprices(:,1)),4);
for n=1:4
    hprchanges(:,n)=hpr(:,n)-yields(:,1);
end 

% figure(1)
% plot(fwdspread)
% 
% figure(3)
% plot(hprchanges)


% compute OLS standard error
B2 = NaN(4,2);
TSTAT2 = NaN(4,2);
S22 = NaN(4,1);
SE22 = NaN(4,2);
SEWHITE2 = NaN(4,2);
VCV2= NaN(2,2);
VCVWHITE2 = NaN(2,2);
for i=1:4
   [B2(i,:),TSTAT2(i,:),S22(i,1),VCV2,VCVWHITE2] = ols(hprchanges(1:end-12,i),fwdspread(1:end-12,i));
   SE2(i,:)=sqrt(diag(VCV2));
   SEWHITE2(i,:)=sqrt(diag(VCVWHITE2));
end

OLS2 = table(B2(:,2),SE2(:,2),TSTAT2(:,2));
OLS2.Properties.VariableNames = {'Gamma' 'OLS_SE' 'T_stat' };
disp(OLS2);

OLSWHITE2 = table(B2(:,2),SEWHITE2(:,2),TSTAT2(:,2));
OLSWHITE2.Properties.VariableNames = {'Gamma' 'OLS_SEWHITE' 'T_stat' };
disp(OLSWHITE2);

% compute HAC standard erros
B2_hac = NaN(4,2);
TSTAT2_hac = NaN(4,2);
S22_hac = NaN(4,1);
SE2_hac = NaN(4,2);
VCVNW2_hac = NaN(2,2);

for i=1:4
   [B2_hac(i,:),TSTAT2_hac(i,:),S22_hac(i,1),VCVNW2_hac] = olsnw(hprchanges(1:end-12,i),fwdspread(1:end-12,i));
   SE2_hac(i,:)=sqrt(diag(VCVNW2_hac));
end

HAC2 = table(B2_hac(:,2),SE2_hac(:,2),TSTAT2_hac(:,2));
HAC2.Properties.VariableNames = {'Gamma_hac' 'Hac_SE' 'T_stat' };
disp(HAC2);
