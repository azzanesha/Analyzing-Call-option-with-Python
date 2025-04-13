#!/usr/bin/env python
# coding: utf-8

# ## Group Number: 16
# ### Topic: Bull Spread
# #### Submitted By: Md Ali Azzan Esha - 234935, Md Shahin Sheikh - 234918

# # Time Series Analysis of Fedex Corp (FDX)

# In[2]:


import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import yfinance as yf


# In[115]:


FDXdata= yf.download("FDX", start= "2018-01-01", end= "2018-12-31")


# In[116]:


FDXdata.head(5)


# In[117]:


data1= pd.DataFrame(FDXdata["Adj Close"]).dropna()
data1.head(5)


# In[119]:


plt.figure(figsize=(15,5))
plt.plot(data1["Adj Close"], "b", linewidth= 3.0)
plt.title ("Adjusted Closing Price")
plt.ylabel("Price")
plt.xlabel("Year")
plt.grid(True)


# In[135]:


data1["Logarithmic Returns"]= np.log(data1["Adj Close"]/data1["Adj Close"].shift(1))
data1.tail(3)


# In[120]:


plt.figure(figsize=(15,5))
plt.plot(data1["Logarithmic Returns"], "blue", linewidth= 3.0)
plt.title ("Logarithmic Returns of FDX")
plt.ylabel("Price")
plt.xlabel("Year")
plt.grid(True)


# In[121]:


D_mean= np.mean(data1['Logarithmic Returns'])
Y_mean= D_mean*252 


# In[122]:


D_sigma= np.std(data1['Logarithmic Returns'])
Y_sigma= D_sigma*252


# In[123]:


print("Daily Mean: ", round (D_mean,3))
print("Yearly Mean: ", round (Y_mean,3)) 
print("Daily Sigma: ", round (D_sigma,3))
print("Yearly Sigma: ", round (Y_sigma,3))


# In[124]:


number_bin = 30

Xlabel= np.arange(-0.15,0.15,0.001)
plt.hist(data1["Logarithmic Returns"], bins= number_bin)
plt.plot(Xlabel, norm.pdf(Xlabel,D_mean, D_sigma))
plt.grid(True)


# ### Note: As FDX logarithmic return perfectly belongs to normal distribution graph, so we can say that it follows Geometric Brownian Motion (GBM).

# # Bull Spread

# In[164]:


r     = 0.0377               #riskfree rate of USA 
S     = 248.62               #underlying price of FDX (collected from NASDAQ)
K_Low = 260                  #strike price of FDX (we have forecasted this vale)_Lower strike price
T     = 252/365              #we are assuming 252 trading days in a year
sigma = 0.2120               #volatility has considered based on historical data 


# In[165]:


def blackscholes(r, S, K_Low, T, sigma, type= "Call Option"):
    d1= (np.log(S/K_Low) + (r + sigma**2/2)*T)/(sigma*np.sqrt(T))
    d2= d1 - sigma*np.sqrt(T)
    try:
        if type == "Call Option":
            price = S*norm.cdf(d1, 0, 1) - K_Low*np.exp(-r*T)*norm.cdf(d2, 0, 1)
        elif type == "Put Option":
            price = K_Low*np.exp(-r*T)*norm.cdf(-d2, 0, 1) - S*norm.cdf(-d1, 0, 1)
        return price
    except:
        print("Recheck option parameters")


# In[166]:


print("Call Option Price with Lower Strike Price: ", round(blackscholes(r, S, K_Low, T, sigma, type= "Call Option"), 2))


# In[167]:


print("Put Option Price is with Lower Strike Price: ", round(blackscholes(r, S, K_Low, T, sigma, type= "Put Option"), 2))


# In[180]:


r      = 0.0377               #riskfree rate of USA 
S      = 248.62               #underlying price of FDX (collected from NASDAQ)
K_high = 265                  #strike price of FDX (we have forecasted this vale)_Higher Strike Price
T      = 252/365              #we are assuming 252 trading days in a year
sigma  = 0.2720               #volatility has considered based on historical data 


# In[181]:


def blackscholes(r, S, K_high, T, sigma, type= "Call Option"):
    d1= (np.log(S/K_high) + (r + sigma**2/2)*T)/(sigma*np.sqrt(T))
    d2= d1 - sigma*np.sqrt(T)
    try:
        if type == "Call Option":
            price = S*norm.cdf(d1, 0, 1) - K_high*np.exp(-r*T)*norm.cdf(d2, 0, 1)
        elif type == "Put Option":
            price = K_high*np.exp(-r*T)*norm.cdf(-d2, 0, 1) - S*norm.cdf(-d1, 0, 1)
        return price
    except:
        print("Recheck option parameters")


# In[182]:


print("Call Option Price is with Higher Strike Price: ", round(blackscholes(r, S, K_high, T, sigma, type= "Call Option"), 2))


# In[183]:


print("Put Option Price is with Higher Strike Price: ", round(blackscholes(r, S, K_high, T, sigma, type= "Put Option"), 2))


# In[186]:


call_spread= 18.33 - 15.36
print ("Bull Call Spread: ", round (call_spread, 2))


# In[187]:


Put_Spread= 27.9 - 20.06
print ("Bull Put Spread: ",round (Put_Spread, 2))

