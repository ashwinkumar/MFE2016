%%Homework 1 Derivatives
%%by YiChi Chan, Sally , George , Ashwin
%%
%%Q1. a) Straddle construction using T=4, r=0.02 , h=0.25 , u = exp(rh+
%%0.2*sqrt(h)) d = exp(rh - 0.2*sqrt(h)), delta = 0, S0 = 90, K= 100
clear;
h = 0.25;
T = 4*h;
r = 0.02;
N = ceil(T/h);
u = exp(r*h + 0.2*sqrt(h));
d = exp(r*h -0.2*sqrt(h));
delta = 0;
S0 = 100;
K=90;
discrete_div= zeros(1,N+1);
%[EuropeanCallPrice_1, Delta_Call_1, B_Call_1,] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'EC');
%[EuropeanPutPrice_1, Delta_Put_1, B_Put_1,] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'EP');
[EuropeanStraddle1,EuropeanDelta_Straddle1,EuropeanBond_Straddle1,] =OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'ES');
disp(sprintf('European Straddle Option Price = %f',EuropeanStraddle1(1,1)));
% Straddle
%%Q1. b) Straddle construction using T=40, r=0.02 , h=0.025 , u = exp(rh+
%%0.2*sqrt(h)) d = exp(rh - 0.2*sqrt(h)), delta = 0, S0 = 90, K= 100
clear;
h = 0.025;
T = 40*h;
r= 0.02;
N = ceil(T/h);
u = exp(r*h + 0.2*sqrt(h));
d = exp(r*h -0.2*sqrt(h));
delta = 0;
S0 = 100;
K=90;
discrete_div= zeros(1,N+1);
%[EuropeanCallPrice_2, Delta_Call_2, B_Call_2] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'EC');
%[EuropeanPutPrice_2, Delta_Put_2, B_Put_2] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'EP');
[EuropeanStraddle2,EuropeanDelta_Straddle2,EuropeanBond_Straddle2,] =OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'ES');
disp((sprintf('European Straddle2 Option Price = %f',EuropeanStraddle2(1,1))));

%%Q1. c) Binary Call Option construction using T=4, r=0.02 , h=0.25 , u = exp(rh+
%%0.2*sqrt(h)) d = exp(rh - 0.2*sqrt(h)), delta = 0, S0 = 90, K= 100
h=0.25;
T= 4*h;
r=0.02;
N= ceil(T/h);
u = exp(r*h + 0.2*sqrt(h));
d = exp(r*h -0.2*sqrt(h));
delta = 0;
S0 = 100;
K=90;
discrete_div= zeros(1,N+1);
[EuropeanBinaryCallPrice, Delta_Call_Binary, B_Call_Binary,] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'EB');
disp(sprintf('European Binary Call Option Price = %f',EuropeanBinaryCallPrice(1,1)));
%%Q2. American Option 
clear;
h = 1/365;
T = 250/365;
N = ceil(T/h);
S0 = 10;
K = 10;
r= 0.01;
delta = 0;
discrete_div= zeros(1,N+1);
u = exp(r*h + 0.15*sqrt(h));
d = exp(r*h -0.15*sqrt(h));
[AmericanCallOption, Delta_Call_American, B_Call_American,Exercise_Call] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'AC');
[AmericanPutOption, Delta_Put_American, B_Put_American,Exercise_Put] = OptionPricing(S0,K,u,d,N,r,delta,h,discrete_div,'AP');

disp(sprintf('American Put Option Price = %f',AmericanPutOption(1,1)));
disp(sprintf('American Call Option Price = %f',AmericanCallOption(1,1)));
%%Q3. Dividend payments
clear;
K =10;
r= 0.02;
S0 =10;
h = 1/365;
T = 200/365;
N = ceil(T/h);
u = exp(0.2*sqrt(h));
d = 1/u;
discrete_div= zeros(1,N+1);
DD = 50:50:150;
for i= 1:length(DD)
    discrete_div(DD(i)+1) = 0.05;
end


[AmericanCallOption_DD] = OptionPricing_RiskNeutral(S0,K,u,d,N,r,h,discrete_div,'AC');
[AmericanPutOption_DD] = OptionPricing_RiskNeutral(S0,K,u,d,N,r,h,discrete_div,'AP');

[AmericanStraddleOption_DD] = OptionPricing_RiskNeutral(S0,K,u,d,N,r,h,discrete_div,'AS');
disp(sprintf('American Straddle Option Price = %f',AmericanStraddleOption_DD(1,1)));
disp(sprintf('American Put + Call Option Price = %f',AmericanPutOption_DD(1,1) + AmericanCallOption_DD(1,1)));

%The reason for the straddle to be lower than the combined put and call is in 
%individual put/call, we have option to exercise either one of them.'); 
%But for a straddle they have to be exercised together. Also the payoff of
%a straddle american option is max(c(t) + p(t), mod(S(t)-K)) which is
%lesser than max(c(t),S(t)-K)) + max(p(t), K-S(t))');

%%Q4. Monte Carlo Simulation
clear all;
%part a
S0   = 200;
K     = 220;
T     = 1;
N     = 365;
r     = .02;
h     = T/N;
sigma = 0.2;
n = 100000; % number of Monte Carlo samples

[AsianOptionPrice, AsianOptionStd] = MC_Asian(S0, K, r, sigma, T, N, n);
lower_bound = AsianOptionPrice - 1.96*AsianOptionStd/sqrt(n);
upper_bound = AsianOptionPrice + 1.96*AsianOptionStd/sqrt(n);
disp(sprintf('Asian Option Estimate Price = %f with 95 perc conf (%f,%f)',AsianOptionPrice,lower_bound,upper_bound));

%%Q5. Pricing American Option using Least Square Method
clear;
S0 = 200;
r= 0.1;
sigma= 0.3;
K=220;
T=1;
N=250;
n = 100000;
LSM_AmericanOptionPrice = LSM_American(S0, K, r, sigma, T, N, n,-1);
disp(LSM_AmericanOptionPrice);
%
Paths = [10 100 1000 10000 100000];
price = zeros;
for i = 1: length(Paths)
    price(i) = LSM_American(S0, K, r, sigma, T,N, Paths(i), -1);
end
figure
plot(log10(Paths),price)
xlabel ('Log(Simulations)');
ylabel ('American Option Price');
Steps = [3 10 100 250 1000];
price2 = zeros;
for i = 1: length(Steps)
    price2(i) = LSM_American(S0, K, r, sigma, T, Steps(i), n, -1);
end
figure
plot(log(Steps),price2);
xlabel('Ln(Steps)');
ylabel('American Option Price');
