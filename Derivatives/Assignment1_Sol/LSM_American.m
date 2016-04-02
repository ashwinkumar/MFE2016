function OptionPrice = LSM_American(S0, K, r, sigma, T, N, n, flag)
% Returns the LSM American Option price with the
% following parameters
% S0  -  Initial Stock Price
% K   -  Strike Price
% sigma - volatality
% N   - Steps
% T   - Years
% n   - number of simulations
% flag- Call/Put
    if flag ==1
        rltiplier = 1; %Call Option
    else
        rltiplier = -1; % Put Option
    end
    h = T/(N);
    drift = r - sigma*sigma/2;
    randn ('seed',0);
    for i =1:n
        %rng default
        %randn ('seed',0);
        S(:,i) = S0*[1;cumprod(exp(drift*h+sigma*sqrt(h)*randn([N 1])))];
    end
    CashFlow = zeros(N, n);
    Stop = ones(1, n)*(N+1);
    % decide the cash flow in last period
    for i=1:n
        CashFlow(N+1, i) = max(rltiplier*(S(N+1, i) - K), 0);
    end
    % find the option in the money in i period
    for i = N:-1:2
        X = zeros;
        Y = zeros;
        count = 1;
        location = zeros;
        for j=1:n
            if rltiplier*(S(i, j) - K) >0
                X(count) = S(i, j);
                Y(count) = CashFlow(Stop(j), j)* exp(-drift*h*(Stop(j)-i));
                location(count) = j;
                count = count + 1;
            end
        end
        count = count - 1;
        % run regression
        if count >0
            M = [ones(1, length(Y)) ;X ; X.^2];
            b = lsqlin(transpose(M), Y, [], []);
            E_Y = transpose(M) * b;
        end
        % decide to early excercise or not
        for j=1:count
            if max(rltiplier*(S(i, location(j)) - K), 0) > E_Y(j)
                CashFlow(i, location(j)) = max(rltiplier*(S(i, location(j)) - K), 0);
                CashFlow(Stop(location(j)), location(j)) = 0;
                Stop(location(j)) = i;
            end
        end
    end
    OptionPrice = 0;
    for i = 1:n
        OptionPrice = OptionPrice + CashFlow(Stop(i), i) * exp(-r*h*(Stop(i)-1));
    end
    OptionPrice = OptionPrice / n;
end