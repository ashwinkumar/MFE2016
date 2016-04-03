%Empirical Method HW3 
clear all
%% Problem 3: AR(1) model of the yield curve
%% (1). estimate different AR(p) models for the 1-month yield
data = xlsread('fama_bliss_data.xlsx');
for i =1:length(data(:,1))
    data(i, 1) = datenum(num2str(data(i,1)), 'yyyymmdd');
end
dropped_data = data(1:(751-72), :);
y = dropped_data(:,2);

p=20;
aic=zeros(p,1);
bic=zeros(p,1);
for i =1:p
    Mdl=arima(i,0,0);
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,y);
    [aic(i,1), bic(i,1)]=aicbic(logL,i,679);   
end
plot(bic)
% I choose AR(1) model because the BIC value is the lowest.
ToEstMdl = arima(1, 0, 0);
ToEstMdl.Constant = 0;
Y = dropped_data(:, 2);
[EstMdl, EstParamCov, logL, info] = estimate(ToEstMdl, Y );
OurModel = arima('AR',{0.992729},'Variance',0.2,'Constant',0);
[E,V] = infer(OurModel,Y);

ts1 = timeseries(E, dropped_data(:, 1));
plot(ts1)
datetick('x','mmm yyyy')
autocorr(E)
%%
%% (2). Estimate an AR(1) model on monthly data for the 1-month yield
Y = dropped_data(:,2);
X = zeros(length(Y), 3);
X(:, 1) = 1;
X(:, 2) = lagmatrix(Y,1);
[parameters, bint, residuals,rint,stats] = regress(Y,X);
mu = parameters(1)/(1-parameters(2));
phi = parameters(2);
sigma = sqrt(var(residuals, 'omitnan'));
%%
%% (3). calibrate the Vasicek model to get the best ?t for the average yield curve in the sample
largest = 60;
a_bar = zeros;
b_bar = zeros;
a_bar(1) = 0;
b_bar(1) = -1;
delta0 = 0;
delta1 = 1;
% initial guess for lamda_0 and lamda_1 from class notes
lamda_0 = -0.1144; 
lamda_1 = -10.741; 
step = 0.1;
step_size = 100;
l0 = lamda_0-step*step_size:step:lamda_0+ step*step_size;
l1 = lamda_1-step*step_size:step:lamda_1+step*step_size;
minsum = 1000;
l0min = 0;
l1min = 0;
Vasicek_Y = zeros;
Vasicek_Y(1) = mean(Y);


for lambda0 = l0
    for lambda1 = l1
        for i =2:largest
            b_bar(i) = b_bar(i-1)*(phi - sigma*lambda1) - delta1;
            a_bar(i) = a_bar(i-1) - delta0 + b_bar(i-1)*((1-phi)*mu - sigma*lambda0) + 1/2*sigma^2*(b_bar(i-1))^2;
        end
        g = Y;
        z = zeros(length(Y),6);
        k = 1;
        for i = [3,12,24,36,48,60]
            z(:,k) = -(a_bar(i)+b_bar(i)*g)/i;
            k = k+1;
        end
        if sum((mean(dropped_data(:,2:8))-mean([g,z]))'.^2) < minsum
            minsum = sum((mean(dropped_data(:,2:8))-mean([g,z]))'.^2);
            l0min = lambda0;
            l1min = lambda1;
        end
    end
end

%Using the l0min and l1min
for i =2:largest
     b_bar(i) = b_bar(i-1)*(phi - sigma*l1min) - delta1;
     a_bar(i) = a_bar(i-1) - delta0 + b_bar(i-1)*((1-phi)*mu - sigma*l0min) + 1/2*sigma^2*(b_bar(i-1))^2;
end

for i =2:largest
    Vasicek_Y(i) = -1/i*a_bar(i) -1/i*b_bar(i)*Vasicek_Y(1);
end
periods = [1 3 12 24 36 48 60];
TermStructure = zeros;
for i=1:length(periods)
    TermStructure(i) = mean(dropped_data(:, i+1));
end
x = 1:60;
plot(x, Vasicek_Y, periods, TermStructure);
%%