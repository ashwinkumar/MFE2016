%%Homework 2 Derivatives
%%by YiChi Chan, Sally Shi, George Bonebrigh, Ashwin Kumar
%%
%%Q1. Black-Scholes Closed Form Solution vs Monte-Carlo Simulation
%%a) Simulate 5 paths for Stock Price
clear;
S0 = 100;
Sb= -1;
N_Days = 360;
T  = 1/4;
K  = 100;
r  = 0.05;
sigma = 0.2;
div_yield = 0;
n_sims = 5;
trading_hrs = 8;
N = 12*trading_hrs*N_Days*T;
dt = (5/60)/trading_hrs/N_Days;
[St,]= GetSimulatedGBMStockPrice(S0,Sb,r,sigma,N,dt,n_sims,'N');
figure;
x = linspace(0,T*N_Days,N+1);

plot(x,St(:,:));
title('5 simulated paths using GBM');
xlabel('Time in Days') % x-axis label
ylabel('Stock Price') % y-axis label
%%b Calculate Black Scholes Call Option Price
[BS_Call_Price] = bs_call_price(S0,K,r,T,sigma,div_yield);
fprintf('BS Call Price = %0.3f \n',BS_Call_Price);

%%c Monte Carlo Simulation for Option Price
n_sims = [100,1000,10000,1000000];
P = zeros;
L = zeros;
U = zeros;
for i=1:length(n_sims)
    [P(i),L(i),U(i)]=MC_Option_Price_Direct(S0,K,r,sigma,T,'EC',n_sims(i));
end
figure; hold on;
p1 =plot(n_sims,U);
p2 =plot(n_sims,L); 
p3 = refline(0,BS_Call_Price);
title('MC Simulation: Option Price vs No of simulations');
xlabel('Number of simulations') % x-axis label
ylabel('MC European Call Option Price') % y-axis label
legend([p1,p2,p3],'Upper bound','Lower Bound','BS Price');
fprintf('We can see from the above graph the length of the confidence interval\n decreases as the number of simulations increases. This satisfies with our understanding\n of standard error which is inversely proportional to number of simulations. Also the MC price\n tends towards the BlackScholes price as N increases\n');
%%2 Down and Out Put Option Price by MC Simulation

clear;
S0= 100;
T = 1/4;
K = 95;
N_Days = 360;
Sb = 75;
r = 0.05;
sigma = 0.2;
delta = 0;
n_sims  = 1000;
trading_hrs = 8;
N = 12*trading_hrs*N_Days*T;
dt = (5/60)/trading_hrs/N_Days;
[St, Within_Barrier,]= GetSimulatedGBMStockPrice(S0,Sb,r,sigma,N,dt,n_sims,'D-O-P');

fprintf('Number of simulations that crossed the barrier = %d\n',n_sims-sum(Within_Barrier));
r_array = r*ones(n_sims,1);
[P1,L1,U1]=MC_Option_Price(St(:,N+1),K,r_array(:,1), T,'EP',Within_Barrier,n_sims);
fprintf('Down and Put Option Estimate Price = %f with 95 perc conf (%f,%f)\n',P1,L1,U1);


[BS_Put_Price] = bs_put_price(S0,K,r,T,sigma,delta);
fprintf('In general we expect the down and put option price to be less when\n compared to the Black-Scholes put price.The reason being the option is knocked\n out for certain number of paths (5 out of 1000). In our simulation we are \n getting the down and put to be higher. The reason for this is the number of\n simulations we are using is low. The down and put option price will be lesser\n than the BS Option price when we increase the simulation. The argument for this\n is simular to the argument that Binomial model option price tends to Black Scholes\n (continuos) model when we increase the number of nodes ( can be visualized\n  as number of simulations for Monte-Carlo');

%b & c Out Paths vs Vol
sigma = 0.05:0.05:0.40;
PercentageOutOfPaths = zeros;
P1 =zeros;
L1 = zeros;
U1 = zeros;
BS_PutPrice = zeros;
r_array = r*ones(n_sims,1);
for i=1:length(sigma)
    [St, Within_Barrier,]= GetSimulatedGBMStockPrice(S0,Sb,r,sigma(i),N,dt,n_sims,'D-O-P');
    PercentageOutOfPaths(i) = 1-sum(Within_Barrier)/length(Within_Barrier);
    [P1(i),L1(i),U1(i)]=MC_Option_Price(St(:,N+1),K,r_array(:,1),T,'EP',Within_Barrier,n_sims);
    BS_PutPrice(i) = bs_put_price(S0,K,r,T,sigma(i),delta);
end
figure;
plot(sigma,PercentageOutOfPaths);
title('Percentage of knocked out option price path vs Volatality');
xlabel('sigma') % x-axis label
ylabel('% of paths knocked out crossing barrier') % y-axis label

fprintf('We can see from the above graph the percentage increases as sigma grows.\n This is due to the fact that as sigma becomes larger the stock is highly volatile. This \n creates more no of paths under MC  simulation that crosses the barrier.');
figure; hold on;
p2 = plot(sigma,P1);
p3 = plot(sigma,L1);
p4 = plot(sigma,U1);
p5 = plot(sigma,BS_PutPrice);
legend([p2,p3,p4,p5],'Price','Upper Bound', 'Lower Bound','BS Put Price');
title('Option Price vs Volatality');
xlabel('sigma') % x-axis label
ylabel('Option Price'); % y-axis label
fprintf(' We can see that BS- option price increases as volatality increases. This is consistent\n with our expectation. For the down and put option there are two opposing effects due to Vol.\n As Vol increases the probability for option to be in money increases and at the same time \n the the probablity to be knocked out (crossing the barrier) also increases. Thus we see the graph to \n change trend after sigma crosses 0.3');
%%3 Exotic Options in Complicated Market Structures
clear;
r0 = 0.05;
beta = 0.05;
alpha = 0.6;
sigma11 = 0.1;
sigma12 = 0.2;
sigma21 = 0.3;
S10 = 10;
S20 = 10;
delta = 0.1;
coeff_dt(1)= alpha;
coeff_dt(2) =beta;
coeff_dw(1) = delta;
T = 1;
dt = 1/250;
N = T/dt;
n_sims = 1000;
r = zeros(n_sims,N+1);
randn('seed',0);

for i=1:n_sims
[r(i,:),] = Euler_Discretization(r0, coeff_dt, coeff_dw ,dt,N, 'R');
end
figure;
histogram(r(:,N+1));
title('Histogram of rates simualted from Euler Discretization');
xlabel('Interest rate') % x-axis label
ylabel('Frequency') % y-axis label


%3 b Simulate one trajectory
T = 100;
dt = 1/52;
N = T/dt;
n_sims = 1000;
r = zeros(n_sims,N+1);
randn('seed',0);
[r,]= Euler_Discretization(r0, coeff_dt, coeff_dw ,dt,N, 'R');
figure;
x = linspace(0,T,N+1);
plot(x,r);
title('Simulated path for interest rates using Euler Discretization');
xlabel('Time in Days') % x-axis label
ylabel('Rate') % y-axis label

%3 c 
T= 0.5;
K=10;
dt = 1/250;
N = T/dt;
n_sims = 10000;

coeff_dt2(1)= 1;
coeff_dw2(1) = sigma11;
coeff_dw2(2) = sigma12;

S1 = zeros(n_sims,N+1);
r = zeros(n_sims,N+1);
randn('seed',0);
for i=1:n_sims
    [r(i,:),rnd]= Euler_Discretization(r0, coeff_dt, coeff_dw ,dt,N, 'R');
    y= vertcat(r(i,:),rnd);
    [S1(i,:),] = Euler_Discretization(S10, coeff_dt2, coeff_dw2 ,dt,N, 'S1',y);
 end
Within_Barrier = ones(n_sims,1);

[P1,L1,U1]=MC_Option_Price(S1(:,N+1),K,mean(r,2),T,'EC',Within_Barrier,n_sims);
% check if discounting using r_0 or r_t
fprintf('European Call Option Estimate Price = %f with 95 perc conf (%f,%f)\n',P1,L1,U1);

% 3 d
coeff_dw3(1) = sigma21;
S1 = zeros(n_sims,N+1);
S2 = zeros(n_sims,N+1);
r = zeros(n_sims,N+1);
randn('seed',0);
for i=1:n_sims
    [r(i,:),rnd]= Euler_Discretization(r0, coeff_dt, coeff_dw ,dt,N, 'R');
    y= vertcat(r(i,:),rnd);
    [S1(i,:),rnd2] = Euler_Discretization(S10, coeff_dt2, coeff_dw2 ,dt,N, 'S1',y);
    y1= vertcat(r(i,:),S1(i,:),rnd);
    [S2(i,:),] = Euler_Discretization(S20, coeff_dt, coeff_dw3 ,dt,N, 'S2',y1);
end
Within_Barrier = ones(n_sims,1);
S1max = max(S1,[],2);
S2max = max(S2,[],2);
St_Calc = max(S1max,S2max);

[P1,L1,U1]=MC_Option_Price(St_Calc,K,mean(r,2),T,'EC',Within_Barrier,n_sims);
% check if discounting using r_0 or r_t
fprintf('Option Estimate Price = %f with 95 perc conf (%f,%f)\n',P1,L1,U1);

%%4 Hedging, Large Price Movements, and Transaction Costs
%a) 
clear;
S0 = 100;
Sb = -1;
N_Days = 365;
T  = 30/365;
K  = 100;
sigma = 0.3;
mu = 0.2;
yield = 0;
trading_hrs = 8;
N = 12*trading_hrs*N_Days*T;
r  = 0.05* ones(1,N+1);
dt = (5/60)/trading_hrs/N_Days;
randn('seed',0);
[S_T, Within_Barrier,]= GetSimulatedGBMStockPrice(S0,Sb,mu,sigma,N,dt,1,'P');

figure;
x = linspace(0,(N)/(trading_hrs*12),N+1);
plot(x,S_T);
title('Simulated path for Stock price using growth rate');
xlabel('Time in Days') % x-axis label
ylabel('Price') % y-axis label
%b) 
T_array = 2*T:-(1/(365*trading_hrs*12)):T;
K_array = K*ones(1,N+1);
[BSCallPrice] = bs_call_price(S_T,K,r,T_array,sigma,yield);
V_0 = BSCallPrice(1);
figure;
x = linspace(0,(N)/(trading_hrs*12),N+1);
plot(x,BSCallPrice);
title('BS Price for the above simulated Stock price using growth rate');
xlabel('Time in Days (60- time to expiration)') % x-axis label
ylabel('Price') % y-axis label

%c)
t_cost= 0;
d1 = getd1(S_T,K,r,T_array,sigma);
[V_T,~,Call_Delta_T]= Get_ReplicatingPortfolio(V_0,S_T,r,dt,d1,N,t_cost,'C');
x = linspace(0,(N)/(trading_hrs*12),N+1);
figure;
plot(x, S_T);
title('Simulated path for Stock price ');
xlabel('Time in Days') % x-axis label
ylabel('Stock Price') % y-axis label

figure;hold on;
p1 = plot(x, BSCallPrice);
p2= plot(x,V_T);
title('Simulated BS Call Price and Replicating Portfolio');
xlabel('Time in Days') % x-axis label
ylabel('Price') % y-axis label
legend([p1,p2],'BS Call Price', 'Replicating Portfolio');
figure;
plot(x, Call_Delta_T);
title('Hedge ratio for simulated stock price');
xlabel('Time in Days') % x-axis label
ylabel('Delta of Call') % y-axis label


%d)
t_cost = 0;
jump_down_v = zeros(1,N+1);
jump_down_v(1,N/2) = -0.1;

randn('seed',0);
[S_T_down_jump,~,~]= GetSimulatedGBMStockPrice(S0,Sb,mu,sigma,N,dt,1,'P',jump_down_v);
[BSCallPrice_down] = bs_call_price(S_T_down_jump,K,r,T_array,sigma,yield);
d1_down = getd1(S_T_down_jump,K,r,T_array,sigma);
V_0_down = BSCallPrice_down(1);
[V_T_down,~,Call_Delta_T_down]= Get_ReplicatingPortfolio(V_0_down,S_T_down_jump,r,dt,d1_down,N,t_cost,'C');
figure;
plot(x, S_T_down_jump);
title('Simulated path for Stock price with -ve jump');
xlabel('Time in Days') % x-axis label
ylabel('Stock Price') % y-axis label

figure;hold on;
p1 = plot(x, BSCallPrice_down);
p2= plot(x,V_T_down);
title('Simulated BS Call Price and Replicating Portfolio');
xlabel('Time in Days') % x-axis label
ylabel('Price') % y-axis label
legend([p1,p2],'BS Call Price', 'Replicating Portfolio');
figure;
plot(x, Call_Delta_T_down);
title('Hedge ratio for simulated stock price');
xlabel('Time in Days') % x-axis label
ylabel('Delta of Call') % y-axis label
fprintf('The stock price with the jump is as shown in the graph. We can see the sudden\n jump down in stock price causes the hedging to be ineffective at T=15 days(jump time). The replicating\n portfolio falls more than the call price for that moment (hedging is not continuos). From the next time( immediately after\n the jump is absorbed )the replicating portfolio closely follows the call price (with diff = diff in drop between call and replicating portfolio). If we try to synthesize\n the call by using our replicating portfolio, we will suffer loss = %.3f\n' ,BSCallPrice_down(N/2+1)-V_T_down(N/2+1));

%e)
t_cost=0;
jump_up_v = zeros(1,N+1);
jump_up_v(1,N/2) = 0.1;
randn('seed',0);
[S_T_up_jump,~,~]= GetSimulatedGBMStockPrice(S0,Sb,mu,sigma,N,dt,1,'P',jump_up_v);
[BSCallPrice_up] = bs_call_price(S_T_up_jump,K,r,T_array,sigma,yield);
d1_up = getd1(S_T_up_jump,K,r,T_array,sigma);

V_0_up = BSCallPrice_up(1);
[V_T_up,~,Call_Delta_T_up]= Get_ReplicatingPortfolio(V_0_up,S_T_up_jump,r,dt,d1_up,N,t_cost,'C');
figure;
plot(x, S_T_up_jump);
title('Simulated path for Stock price with -ve jump');
xlabel('Time in Days') % x-axis label
ylabel('Stock Price') % y-axis label

figure;hold on;
p1 = plot(x, BSCallPrice_up);
p2= plot(x,V_T_up);
title('Simulated BS Call Price and Replicating Portfolio');
xlabel('Time in Days') % x-axis label
ylabel('Price') % y-axis label
legend([p1,p2],'BS Call Price', 'Replicating Portfolio');
figure;
plot(x, Call_Delta_T_up);
title('Hedge ratio for simulated stock price');
xlabel('Time in Days') % x-axis label
ylabel('Delta of Call') % y-axis label
fprintf('We can see similar observation. Even though stock\n jumps up, due to limitation of non\n continuos time hedging, we will not able to track the call price \n during the jump.\n The differnce in Call Price is abosrbed slowly and the replicating portfolio falls back again.\n');
%f)
t_cost =20;
[V_T,B_T,Call_Delta_T]= Get_ReplicatingPortfolio(V_0,S_T,r,dt,d1,N,t_cost,'C');
figure;hold on;
p1 = plot(x, BSCallPrice);
p2 = plot(x,V_T);
title('Simulated BS Call Price and Replicating Portfolio');
xlabel('Time in Days') % x-axis label
ylabel('Price') % y-axis label
legend([p1,p2],'BS Call Price', 'Replicating Portfolio');
sprintf('With the 20 basis points transaction cost,the replicating portfolio value diverges from \n the call value. This explanation falls with the fact that\n as time increases due to the volatile stock price movements\n we end up adjusting our portfolio more number of times\n' );
figure;
plot(x, Call_Delta_T);
title('Hedge ratio for simulated stock price');
xlabel('Time in Days') % x-axis label
ylabel('Delta of Call') % y-axis label
%5 Heston model

%a) 
clear;
S0= 100;
v_0 = 0.01;
K= 100;
r= 0.04;
T=0.5;
trading_hrs = 8;
%dt is 5 minute interval
N = 12*trading_hrs*365*T ;
dt = 1/N;
lambda=0;
rho = -0.5;
kappa = 6;
theta = 0.02;
div_yield = 0;
 
sigma= 0.3;
% Assuming mu = 0.2
mu = 0.2;
randn('seed',0);
[S_T,v_T] = CIRProcess(S0,v_0,mu,sigma,rho,kappa, theta,N,dt);
figure;
x = linspace(0,(N)/(trading_hrs*12),N+1);
plot(x,S_T);
title('Stock price simulated using CIR volatality model');
xlabel('Time in Days') % x-axis label
ylabel('Price') % y-axis label
figure;
plot(x,v_T);
title('CIR volatality model');
xlabel('Time in Days') % x-axis label
ylabel('Volatality') % y-axis label
T_array = T:-(1/(365*trading_hrs*12)):0;
BSCall_Price = bs_call_price(S_T(1,1),K,r,T_array(1,1),sqrt(theta),div_yield); 


HestonCallPrice = HestonModel(S_T(1,1),K,v_T(1,1),r,rho,sigma,kappa,lambda,theta,T_array(1,1));
fprintf('BS Call Option Price = %f Heston Call Option Price =%f\n',BSCall_Price,HestonCallPrice);


%b) 
S01 = 70;
S02 = 130;
S = S01:1:S02;
H_Call = zeros;
BS_Call = zeros;
for i=1:length(S)
    H_Call(i) = HestonModel(S(1,i),K,v_T(1,1),r,rho,sigma,kappa,lambda,theta,T);
    BS_Call(i) = bs_call_price(S(1,i),K,r,T,sqrt(theta),div_yield); 
end
figure;hold on;
p1 = plot(S,H_Call);
p2=plot(S,BS_Call);
title('Heston vs BS Call Price');
xlabel('Initial Stock Price') % x-axis label
ylabel('Call Price') % y-axis label
legend([p1,p2],'Heston Option Price','BS Option Price');
figure;
plot(S,H_Call- BS_Call);
title('Heston - BS Call Price vs Stock Price');
xlabel('Initial Stock Price') % x-axis label
ylabel('Diff in Call Price') % y-axis label


sprintf('We can see from the graph, the Heston model for option price converges to BS option pricing model, for in the money call (S0 > K)\n. The reason for this is Heston model incorporates time varying volatility with negative rho. This makes out of money\n call options (low stock price) to have lesser probablity to be in the money when compared\n to in the money (for -ve rho). CIR model for vol being mean reverting.\n'); 
%c)
S01 = 60;
S02 = 160;
S = S01:1:S02;
rho =-0.5;
implied_vol = zeros(1,length(S));
for i=1:length(S)
    H_Call(i) = HestonModel(S(1,i),K,v_T(1,1),r,rho,sigma,kappa,lambda,theta,T);
    implied_vol(i)= getHS_Call_Implied_Vol(S(1,i),K,r,T,div_yield,H_Call(i),sigma);
end

figure;
plot(S,implied_vol);
title('Implied Volatality');
xlabel('Stock Price') % x-axis label
ylabel('Implied volatality') % y-axis label
%d
rho=0.5;
implied_vol = zeros(1,length(S));
for i=1:length(S)
    H_Call(i) = HestonModel(S(1,i),K,v_T(1,1),r,rho,sigma,kappa,lambda,theta,T);
    implied_vol(i)= getHS_Call_Implied_Vol(S(1,i),K,r,T,div_yield,H_Call(i),sigma);
end

figure;
plot(S,implied_vol);
title('Implied Volatality');
xlabel('Stock Price') % x-axis label
ylabel('Implied volatality') % y-axis label

fprintf('The volatality smile is right skewed for negative rho and left skewed for positive rho.\n
In the plot, we used stock price from 60 to 160 to get the detail of the full vol plot.\n
When rho <0, then the Heston model minimum occurs at Stcok price close to K (from below),\n
When rho >0, then the Heston model minimum occurs at Stock price close to K (from above)');

