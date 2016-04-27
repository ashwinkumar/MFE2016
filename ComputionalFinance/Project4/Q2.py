from pandas_datareader import data as web
from pandas_datareader.data import Options
import numpy as np
import BinomialModel as bm
from datetime import datetime
from scipy.optimize import fsolve
import BlackScholesOptionPrice as bp

### Question 2

def historical_volatility(sym, days):
    "Return the annualized stddev of daily log returns of `sym`."
    try:
        quotes = web.DataReader(sym, 'yahoo')['Adj Close'][-days:]
    except Exception:
        print('Error getting data for symbol %s' %format(sym))
        return None
    #logreturns = np.log(quotes / quotes.shift(1))
    returns = quotes/quotes.shift(1)-1
    return np.sqrt(252*returns.var())

def last_price(sym):
    "Return the last traded price of `sym`."
    try:
        quotes = web.DataReader(sym, 'yahoo')['Close'][-1:]
    except Exception:
        print('Error getting data for symbol %s' % format(sym))
        return None
    return quotes.values[0]


def EuropeanCallPayOff(St,X):
    array_zeros = np.zeros(len(St))
    return np.maximum(St-X,array_zeros)


no_hist_days = 60*21
ticker = "GOOG"
src = "google"

r = 0.02
S0 = last_price(ticker)
sigma = historical_volatility(ticker,no_hist_days)
X = (int)(1.1*S0) - (int)(1.1*S0)%10


expiry = datetime(2017, 1, 17)
'''
options_ticker = Options(ticker, src)

options_data = options_ticker.get_call_data(expiry=expiry)

#options_data = options_ticker.get_all_data()
options_data = options_data.loc[(X, 2017, 'call'), :]


print(options_data.iloc[0:5:, 0:5])


#curr_date =  (time.strftime("%m/%d/%Y"))
#goog = web.DataReader("^GOOG", data_source="yahoo", start=curr_date, end=curr_date)
'''

curr_date = datetime.now()
t = ((expiry- curr_date).days)/365

n= 100
dt = t/n


u,d,p = bm.getBinomialParameters(r,sigma,dt,4)
bin_model = bm.BinomialModel(S0,r, t)

bin_model.n= n
if bin_model.verify(u, d,p):
    call_price = bin_model.getOptionPrice(EuropeanCallPayOff,X,"E")
bin_model.resetParams()


print('Q2. a) Using the Method-4 the call option price [S0,X,r,sigma,t] = [%f,%f,%f,%f,%f] = %f' %(S0,X,r,sigma,t,call_price))

call_price_data = 29.30
implied_vol =fsolve(lambda x: bp.EuropeanCallPrice(S0,X,r,x,t)- call_price_data, sigma)*100

print('Q2. b) The implied volatality of the option price = %f = %f %%' %(call_price_data,implied_vol))
