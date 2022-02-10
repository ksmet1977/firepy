# -*- coding: utf-8 -*-
"""
FIRE calculator
---------------

Created on Tue Feb  8 12:56:05 2022

@author: u0032318
"""
import matplotlib.pyplot as plt
import numpy as np
import copy 

__all__ = ['fire_calculator']

def fire_calculator_(prefire_salary = 40000, prefire_salary_increase = 1, 
                    prefire_extra_income = 0,
                    prefire_expenses = 16000, 
                    retirement_income_net = 1550*12, retirement_income_index_rate = None, retirement_expenses = None, retirement_age = 67,
                    fire_expenses = 20000, monthly = False,
                    capital_start = 450000, capital_allocation = [60,30,10,0], capital_returns = [5,1,0,5],
                    capital_extra = 330000, capital_extra_age = None, real_estate_extra = 200000,
                    tax_rates = [0.35, 0.12, 0, 1.32+0.22],
                    real_estate_value = 200000, real_estate_return = 0,
                    withdrawal_rate = 3.5, inflation = 3.5,
                    age = 45, end_age = 85, fire_age = 46,
                    plot = True, ax = None):
    
    # Convert input to yearly values:
    if monthly:
        prefire_salary = prefire_salary*12
        prefire_expenses = prefire_expenses*12
        fire_expenses = fire_expenses*12
        retirement_income_net = retirement_income_net*12
        
    # Number of years to run simulation for:
    years = end_age - age 
    
    # renormalize allocation to 100 % (= 1):
    capital_allocation = np.array(capital_allocation)
    capital_allocation = capital_allocation/capital_allocation.sum()
    
    # allocate start capital over stocks, bonds, cash, etf:
    portfolio = {'years': [0], 
               'stocks': [capital_start*capital_allocation[0]],
               'bonds': [capital_start*capital_allocation[1]],
               'cash': [capital_start*capital_allocation[2]],
               'etf': [capital_start*capital_allocation[3]],
               'real_estate' : [real_estate_value],
               'total_without_real_estate' : [capital_start],
               'total_with_real_estate' : [capital_start + real_estate_value]}

    if plot:
        if ax is None: 
            fig, ax = plt.subplots(1,1)
        ax.set_xlabel('Age (years)')
        ax.set_ylabel('Total capital (euro)')
    
    # Simulate yearly returns:
    bankrupt_without_re = False
    bankrupt_with_re = False
    fire_without_re = False
    fire_with_re = False
    fire_number = (fire_expenses)/(withdrawal_rate/100)
    if retirement_income_index_rate is None: retirement_income_index_rate = inflation
    if retirement_expenses is None: retirement_expenses = fire_expenses
    for yr in range(1,years+1):
        
        # yearly savings:
        if age + yr < fire_age:
            savings = prefire_salary*(1+prefire_salary_increase/100)**yr - prefire_expenses*(1+inflation/100)**yr
            expenses_infl_corr = prefire_expenses*(1+inflation/100)**yr
            income_net_infl_corr = prefire_salary*(1+prefire_salary_increase/100)**yr
        else:
            savings = 0.0 - fire_expenses*(1+inflation/100)**yr
            expenses_infl_corr = fire_expenses*(1+inflation/100)**yr
            income_net_infl_corr = 0 
            
        # add retirement income from state:
        if age + yr >= retirement_age:
            savings = retirement_income_net*(1+inflation/100)**yr - retirement_expenses*(1+inflation/100)**yr
            income_net_infl_corr = retirement_income_net*(1+inflation/100)**yr

        
        print('Age = {:1.0f} / savings = {:1.0f} / prefire_expenses = {:1.0f} / ret. income = {:1.0f} / inflation = {:1.2f}'.format(age + yr, savings, expenses_infl_corr, income_net_infl_corr, (1+inflation/100)**yr))
           
        # add additional capital & real_estate at specific age:
        if capital_extra_age is not None:
            if age + yr == capital_extra_age:
                capital_extra_corr = []
                for i in range(len(capital_allocation)):
                    capital_extra_corr.append(capital_extra*capital_allocation[i]*(1+capital_returns[i]/100)**yr)
                real_estate_extra_infl_corr = real_estate_extra*(1+real_estate_return/100)**yr
            else:
                capital_extra_corr = [0]*len(capital_allocation)
                real_estate_extra_infl_corr = 0
        else:
            capital_extra_corr = [0]*len(capital_allocation)
            real_estate_extra_infl_corr = 0
            
        # calculate portfolio increases for yr:
        portfolio['years'] += [yr]
        portfolio['stocks'] += [(portfolio['stocks'][-1]*(1+capital_returns[0]/100) + savings*capital_allocation[0]*(1-np.sign(savings)*tax_rates[0]/100)) + capital_extra_corr[0]]
        portfolio['bonds'] += [(portfolio['bonds'][-1]*(1+capital_returns[1]/100) + savings*capital_allocation[1]*(1-np.sign(savings)*tax_rates[1]/100)) + capital_extra_corr[1]]
        portfolio['cash'] += [(portfolio['cash'][-1]*(1+capital_returns[2]/100) + savings*capital_allocation[2]*(1-np.sign(savings)*tax_rates[2]/100)) + capital_extra_corr[2]]
        portfolio['etf'] += [(portfolio['cash'][-1]*(1+capital_returns[3]/100) + savings*capital_allocation[3]*(1-np.sign(savings)*tax_rates[3]/100)) + capital_extra_corr[3]]
        portfolio['real_estate'] += [(portfolio['real_estate'][-1]*(1+real_estate_return/100)) + real_estate_extra_infl_corr ]
        portfolio['total_without_real_estate'] += [(portfolio['stocks'][-1] + portfolio['bonds'][-1] + portfolio['cash'][-1] +  portfolio['etf'][-1])]
        portfolio['total_with_real_estate'] += [(portfolio['total_without_real_estate'][-1] + portfolio['real_estate'][-1])]
        
        if plot:
            label_without = 'total_without_real_estate' if yr == 1 else '_no_legend_'
            label_with = 'total_with_real_estate' if yr == 1 else '_no_legend_'
            ax.plot(age + yr, portfolio['total_without_real_estate'][-1],'bo', label = label_without)
            ax.plot(age + yr, portfolio['total_with_real_estate'][-1],'r.', label = label_with) 
            if yr == 1: ax.hlines(0, xmin=age, xmax = end_age, color = 'k', linestyle = '-')
    
        if (portfolio['total_with_real_estate'][-1] < 0) & (bankrupt_with_re == False): 
            print('Bankrupt at age = {:1.0f} years (with real_estate)'.format(age+yr))
            bankrupt_with_re = True
        if (portfolio['total_without_real_estate'][-1] < 0) & (bankrupt_without_re == False): 
            print('Bankrupt at age = {:1.0f} years (without real_estate)'.format(age+yr))
            bankrupt_without_re = True
            
        if (portfolio['total_with_real_estate'][-1] >= fire_number) & (fire_with_re == False): 
            print('Fire at age = {:1.0f} years (with real_estate)'.format(age+yr))
            fire_with_re = True
        if (portfolio['total_without_real_estate'][-1] >= fire_number) & (fire_without_re == False): 
            print('Fire at age = {:1.0f} years (without real_estate)'.format(age+yr))
            fire_without_re = True
            
    if plot:
        ax.hlines(fire_number,xmin=age, xmax = end_age, color = 'r', linestyle = '--', label = 'SWR = {:1.2f}% -> fire at {:1.0f} euro'.format(withdrawal_rate,fire_number))
        ax.legend(loc='upper left')

    return portfolio

class Portfolio:
    def __init__(self, year = 0, aged = False,
                 capital = 0, capital_allocation = [60,30,10,0,0], capital_returns = [5,1,0,5,0], capital_tax = [0.35,0.12,0,1.32+0.22,0],
                 real_estate = 0, real_estate_returns = 0, real_estate_tax = 0,
                 ):
        self.aged = aged
        self.year = year
        self.capital = capital 
        self.capital_allocation = capital_allocation

        self.stocks = capital * capital_allocation[0]/100
        self.bonds = capital * capital_allocation[1]/100
        self.cash =  capital * capital_allocation[2]/100
        self.etfs = capital * capital_allocation[3]/100
        self.crypto = capital * capital_allocation[4]/100
        self.real_estate = real_estate
        
        self.total_without_real_estate = self.capital
        self.total_with_real_estate = self.capital + self.real_estate
        
        self.capital_returns = capital_returns
        self.capital_tax = capital_tax
        self.real_estate_returns = real_estate_returns
        self.real_estate_tax = real_estate_tax
    
    def __repr__(self):
        return """year = {:1.0f} (aged portfolio = {})
capital = {:1.1f} euro

stocks ({:1.1f}%) = {:1.1f} euro
bonds  ({:1.1f}%) = {:1.1f} euro
cash   ({:1.1f}%) = {:1.1f} euro
etfs   ({:1.1f}%) = {:1.1f} euro
crypto ({:1.1f}%) = {:1.1f} euro
real estate = {:1.1f} euro

Total (without real estate) = {:1.1f} euro
Total (with estate estate)  = {:1.1f} euro
""".format(self.year,self.aged, self.capital, 
self.capital_allocation[0], self.stocks,
self.capital_allocation[1], self.bonds,
self.capital_allocation[2], self.cash,
self.capital_allocation[3], self.etfs,
self.capital_allocation[4], self.crypto,
self.real_estate,
self.total_without_real_estate, self.total_with_real_estate)
    
    def __str__(self):
        return self.__repr__()
        
    def __add__(self, new):
        capital = self.capital + new.capital
        stocks = self.stocks + new.stocks
        bonds = self.bonds + new.bonds
        cash = self.cash + new.cash 
        etfs = self.etfs + new.etfs
        crypto = self.crypto + new.crypto
        real_estate = self.real_estate + new.real_estate
        allocation = list(100*np.array([stocks,bonds,cash,etfs,crypto])/capital)
        
        if self.year >= new.year:
            year = self.year 
            capital_returns = self.capital_returns 
            real_estate_returns = self.real_estate_returns
            capital_tax = self.capital_tax 
            real_estate_tax = self.real_estate_tax
        else:
            year = new.year 
            capital_returns = new.capital_returns 
            real_estate_returns = new.real_estate_returns
            capital_tax = new.capital_tax 
            real_estate_tax = new.real_estate_tax
            
        
        return Portfolio(year = year, capital = capital,
                         capital_allocation = allocation, real_estate = real_estate,
                         capital_returns = capital_returns, capital_tax = capital_tax,
                         real_estate_returns = real_estate_returns, real_estate_tax = real_estate_tax)
    
    def __sub__(self, new):
        capital = self.capital - new.capital
        stocks = self.stocks - new.stocks
        bonds = self.bonds - new.bonds
        cash = self.cash - new.cash 
        etfs = self.etfs - new.etfs
        crypto = self.crypto - new.crypto
        real_estate = self.real_estate - new.real_estate
        allocation = list(100*np.array([stocks,bonds,cash,etfs,crypto])/capital)
        
        if self.year >= new.year:
            year = self.year 
            capital_returns = self.capital_returns 
            real_estate_returns = self.real_estate_returns
            capital_tax = self.capital_tax 
            real_estate_tax = self.real_estate_tax
        else:
            year = new.year 
            capital_returns = new.capital_returns 
            real_estate_returns = new.real_estate_returns
            capital_tax = new.capital_tax 
            real_estate_tax = new.real_estate_tax
        
        return Portfolio(year = year, capital = capital,
                         capital_allocation = allocation, real_estate = real_estate,
                         capital_returns = capital_returns, capital_tax = capital_tax,
                         real_estate_returns = real_estate_returns, real_estate_tax = real_estate_tax)
    
    def scale_asset_value(self, scale = 1):
        self.capital *= scale
        self.stocks *= scale
        self.bonds *= scale
        self.cash *= scale
        self.etfs *= scale
        self.crypto *= scale
        self.real_estate *= scale
        self.total_without_real_estate *= scale
        self.total_with_real_estate *= scale
        
    def get_scaled(self, scale = 1):

        return Portfolio(year = self.year, capital = self.capital*scale,
                         capital_allocation = self.capital_allocation, real_estate = self.real_estate,
                         capital_returns = self.capital_returns, capital_tax = self.capital_tax,
                         real_estate_returns = self.real_estate_returns, real_estate_tax = self.real_estate_tax)

    
    def get_aged(self, years, capital_returns = None, real_estate_returns = None, 
                 capital_allocation = None, rebalance = True):
        if capital_returns is None: capital_returns = self.capital_returns
        if capital_allocation is None: capital_allocation = self.capital_allocation
        if real_estate_returns is None: real_estate_returns = self.real_estate_returns

        stocks = self.stocks * (1+capital_returns[0]/100)**years
        bonds = self.bonds * (1+capital_returns[1]/100)**years
        cash = self.cash * (1+capital_returns[2]/100)**years 
        etfs = self.etfs * (1+capital_returns[3]/100)**years
        crypto = self.crypto * (1+capital_returns[4]/100)**years
        real_estate = self.real_estate * (1+real_estate_returns/100)**years
        capital = stocks + bonds + cash + etfs + crypto
    
        if rebalance == False: 
            # calculate new allocation:
            capital_allocation = list(100*np.array([stocks,bonds,cash,etfs,crypto])/capital)
        else:
            # get cost (tax) of redistributing:
            c_v = np.array([stocks,bonds,cash,etfs,crypto]) # value
            c_a = c_v/capital # allo
            capital_alloc = (np.array(self.capital_allocation)/100) # desired allocation
            c_d = c_a - capital_alloc # difference between current and desired allocation
            cost = 0 
            for i in range(len(c_a)):
                if np.abs(c_d[i]) > 0: # sell / buy part of asset
                    tmp = (c_v[i]/c_a[i]) * np.abs(c_d[i]) 
                    cost = cost + tmp*(self.capital_tax[i]/100) # add to tax cost

            capital = capital - cost # update capital to distribute
            c_v = capital * capital_alloc # allocate capital to assets
            stocks,bonds,cash,etfs,crypto = list(c_v) # split
            capital_allocation = list(100*np.array([stocks,bonds,cash,etfs,crypto])/capital)
            
        return Portfolio(year = self.year + years, capital = capital, aged = True,
                         capital_allocation = capital_allocation, real_estate = real_estate)

        
    def buy(self, capital=0, capital_allocation = None, capital_tax = None,
            real_estate = 0, real_estate_tax = None):
        if capital_allocation is None: capital_allocation = self.capital_allocation
        if capital_tax is None: capital_tax = self.capital_tax
        if real_estate_tax is None: real_estate_tax = self.real_estate_tax
        
        if not isinstance(capital,(float,int)):
            capital = np.array(capital).sum()

        # get required buying ratio's from capital to ensure portfolio keeps allocations:
        br = np.array([(c_a/100)/(1-c_t/100) for (c_a,c_t) in zip(capital_allocation, capital_tax)])
        br = br/br.sum()*(1-np.array(capital_tax)/100) # = buying ratio of capital with tax already subtracted for each asset class
        stocks, bonds, cash, etfs, crypto = capital * br 
        capital = stocks + bonds + cash + etfs + crypto # new capital after taxes
        
        self.stocks += stocks
        self.bonds += bonds
        self.cash += cash
        self.etfs += etfs
        self.crypto += crypto 
        self.capital += capital
        self.capital_allocation = 100*np.array([self.stocks, self.bonds, self.cash, self.etfs, self.crypto])/self.capital
        
        real_estate = real_estate * (1 - real_estate_tax/100)
        self.real_estate += real_estate
        
        self.total_without_real_estate = self.capital
        self.total_with_real_estate = self.capital + self.real_estate
        
    def sell(self, capital=0, capital_allocation = None, capital_tax = None,
             real_estate = 0, real_estate_tax = None):
        if capital_allocation is None: capital_allocation = self.capital_allocation
        if capital_tax is None: capital_tax = self.capital_tax
        if real_estate_tax is None: real_estate_tax = self.real_estate_tax
        
        if not isinstance(capital,(float,int)):
            capital = np.array(capital).sum()
            
        # get required selling ratio's from capital to ensure portfolio keeps allocations:
        sr = np.array([(c_a/100)/(1+c_t/100) for (c_a,c_t) in zip(capital_allocation, capital_tax)])
        sr = sr/sr.sum()*(1+np.array(capital_tax)/100) # = selling ratio of capital with tax already subtracted for each asset class

        stocks, bonds, cash, etfs, crypto = capital * sr 
        capital = stocks + bonds + cash + etfs + crypto # new capital after taxes
        
        self.stocks -= stocks
        self.bonds -= bonds
        self.cash -= cash
        self.etfs -= etfs
        self.crypto -= crypto 
        self.capital -= capital
        self.capital_allocation = 100*np.array([self.stocks, self.bonds, self.cash, self.etfs, self.crypto])/self.capital
        
        real_estate = real_estate * (1 + real_estate_tax/100)
        self.real_estate -= real_estate
        
        self.total_without_real_estate = self.capital
        self.total_with_real_estate = self.capital + self.real_estate
        
# rebalance = True
# p1 = Portfolio(year=0,capital = 100000)
# for yr in range(22):
#     # p1.buy(capital = 1000-172.9)
#     p1.buy(capital = 1000+445.173)
#     p1 = p1.get_aged(1,rebalance = rebalance)
#     p1.sell(capital = 10000)
#     p1 = p1.get_aged(1,rebalance = rebalance)
# # p1.sell(capital = 40000)
# # p1 = p1.get_aged(1,rebalance = rebalance)
# # p1.sell(capital = 40000)
# # p1 = p1.get_aged(1,rebalance = rebalance)
# # p1.sell(capital = 22443.7)
# # p1 = p1.get_aged(1,rebalance = rebalance)
# print(p1)
# raise Exception('')        

def _get_portfolio_for_year(year = 0, inflation = 3.5,
                          salary = 0, salary_increase = 1, 
                          extra_income = 0, extra_income_increase = 0,
                          expenses = 0, expenses_increase = 0,
                          
                          portfolio = None,
                          capital = 0, 
                          capital_allocation = [60, 30, 10, 0, 0], 
                          capital_returns = [5,1,0,5,0],
                          capital_tax = [0.35, 0.12, 0, 1.32 + 0.22,0],
                          real_estate = 0, 
                          real_estate_returns = 0, 
                          real_estate_tax = 0,
                          access_usufruct = False,
                          rebalance = True):
    
    capital_year = capital if (((year == 0) & (portfolio is None))  | access_usufruct) else 0 # only on startup (or when usufruct is accessed) is capital (pre-existing) not bought   
    real_estate_year = real_estate if (((year == 0) & (portfolio is None)) | access_usufruct) else 0 # only on startup (or when usufruct is accessed) is real_estate (pre-existing) not bought   

    # calulate (capital) investment for year (>0):
    money_in =  salary*(1 + salary_increase/100)**year + extra_income*(1 + extra_income_increase/100)**year  
    money_out = (expenses*(1 + expenses_increase/100)**year)*(1 + inflation/100)**year # correct expenses for inflation (reduced buying power)
    balance = money_in - money_out # saving or expenses
    buy = balance > 0 # if balance positive buy, else sell    
    
    balance = np.abs(balance)
    
    # renormalize allocation to 100 % (= 1):
    capital_allocation = np.array(capital_allocation)
    capital_allocation = list(100*capital_allocation/capital_allocation.sum())
    
    # create new portfolio for year:
    prtfl = Portfolio(year = year, aged = False, 
                      capital = capital_year,
                      capital_allocation = capital_allocation,
                      capital_returns = capital_returns,
                      capital_tax = capital_tax,
                      real_estate = real_estate_year,
                      real_estate_returns = real_estate_returns,
                      real_estate_tax = real_estate_tax)

    # add existing portfolio to new one:
    if portfolio is not None:
        prtfl = prtfl + portfolio
        
    # Update portfolio to value at end of year:
    prtfl = prtfl.get_aged(1, capital_returns = capital_returns, 
                            real_estate_returns = real_estate_returns, 
                            capital_allocation = capital_allocation, 
                            rebalance = rebalance)  
    
    
    
    # Quick check to see if real estate needs to be sold, sell if necessary:
    if (buy == False) & (prtfl.total_without_real_estate - balance < 0):
        if real_estate == 0: 
            real_estate_sell_value = prtfl.real_estate*(1-real_estate_tax/100)
        else:
            real_estate_sell_value = real_estate*(1-real_estate_tax/100) # sell this one first
            real_estate = 0
        balance = balance + real_estate_sell_value # add to balance
        prtfl.sell(capital = 0, real_estate = real_estate_sell_value) # sell real_estate
        

    # Buy or sell investements at end of year:
    if buy:   
        prtfl.buy(capital = balance, capital_allocation = capital_allocation, capital_tax = capital_tax, 
                  real_estate = real_estate, real_estate_tax = real_estate_tax)
    else:
        prtfl.sell(capital = balance, capital_allocation = capital_allocation, capital_tax = capital_tax,
                   real_estate = real_estate, real_estate_tax = real_estate_tax)

    return prtfl

    

def fire_calculator(age, 
                      end_age = 85,
                      f_age = 45,
                      r_age = 67,
                      u_age = None, # age at which usufruct drops 
                      
                      swr = 3.5,
                      rebalance = True,
                      inflation_adjusted_vals = True, 
                      
                      pf_inflation = 3.5,
                      pf_salary = 0, pf_salary_increase = 1, 
                      pf_extra_income = 0, pf_extra_income_increase = 0,
                      pf_expenses = 0, pf_expenses_increase = 0,
                      pf_capital = 0, 
                      pf_capital_allocation = [60, 30, 10, 0, 0], 
                      pf_capital_returns = [5,1,0,5,0],
                      pf_capital_tax = [0.35, 0.12, 0, 1.32 + 0.22, 0],
                      pf_real_estate = 0, 
                      pf_real_estate_returns = 0, 
                      pf_real_estate_tax = 0,
                      
                      f_inflation = None,
                      f_salary = 0, f_salary_increase = 0, 
                      f_extra_income = 0, f_extra_income_increase = 0,
                      f_expenses = 0, f_expenses_increase = 0,  
                      f_capital = 0, 
                      f_capital_allocation = None, 
                      f_capital_returns = None,
                      f_capital_tax = None,
                      f_real_estate = 0, 
                      f_real_estate_returns = None, 
                      f_real_estate_tax = None,
                      
                      r_inflation = None,
                      r_salary = 0, r_salary_increase = 0, 
                      r_extra_income = 0, r_extra_income_increase = 0,
                      r_expenses = 0, r_expenses_increase = 0,
                      r_capital = 0, 
                      r_capital_allocation = None, 
                      r_capital_returns = None,
                      r_capital_tax = None,
                      r_real_estate = 0, 
                      r_real_estate_returns = None, 
                      r_real_estate_tax = None,
                      
                      u_capital = 0, 
                      u_capital_allocation = None, 
                      u_capital_returns = None,
                      u_capital_tax = None,
                      u_real_estate = 0, 
                      u_real_estate_returns = None, 
                      u_real_estate_tax = None,
                      
                      verbosity = 2, ax = None,
                      leg_dict = {'loc':'upper left'}):
        
    # Use same values as for prefire period:
    if f_inflation is None: f_inflation = pf_inflation
    if f_capital_allocation is None: f_capital_allocation = pf_capital_allocation
    if f_capital_returns is None: f_capital_returns = pf_capital_returns
    if f_capital_tax is None: f_capital_tax = pf_capital_tax
    if f_real_estate_returns is None: f_real_estate_returns = pf_real_estate_returns
    if f_real_estate_tax is None: f_real_estate_tax = pf_real_estate_tax
    
    if r_inflation is None: r_inflation = pf_inflation
    if r_capital_allocation is None: r_capital_allocation = pf_capital_allocation
    if r_capital_returns is None: r_capital_returns = pf_capital_returns
    if r_capital_tax is None: r_capital_tax = pf_capital_tax
    if r_real_estate_returns is None: r_real_estate_returns = pf_real_estate_returns
    if r_real_estate_tax is None: r_real_estate_tax = pf_real_estate_tax

    if u_capital_allocation is None: u_capital_allocation = pf_capital_allocation
    if u_capital_returns is None: u_capital_returns = pf_capital_returns
    if u_capital_tax is None: u_capital_tax = pf_capital_tax
    if u_real_estate_returns is None: u_real_estate_returns = pf_real_estate_returns
    if u_real_estate_tax is None: u_real_estate_tax = pf_real_estate_tax
    
    if pf_salary_increase is None: pf_salary_increase = pf_inflation
    if f_salary_increase is None: f_salary_increase = f_inflation
    if r_salary_increase is None: r_salary_increase = r_inflation

    # prepare plot:
    if verbosity > 1:
        if ax is None: 
            fig, ax = plt.subplots(1,1, sharey = True, figsize = (9,9))
        #ax2 = ax.twinx()
        ax.set_xlabel('Age (years)')
        if inflation_adjusted_vals: 
            ax.set_ylabel('Value (euro) (inflation adjusted)')
        else:
            ax.set_ylabel('Value (euro)')
        #ax2.set_ylabel('Compounded inflation')
        #leg_dict2 = copy.deepcopy(leg_dict)

    years = end_age - age 
    
    fire_number = f_expenses / (swr/100)
    
    if u_age is None: u_age = 1000
    f_prtfl_u = lambda year, inflation: _get_portfolio_for_year(year = year,
                                              inflation = inflation,
                                              salary = 0, 
                                              salary_increase = 0, 
                                              extra_income = 0,
                                              extra_income_increase = 0,
                                              expenses = 0, 
                                              expenses_increase = 0,
                                  
                                              portfolio = prtfl,
                                              capital = u_capital, 
                                              capital_allocation = u_capital_allocation, 
                                              capital_returns = u_capital_returns,
                                              capital_tax = u_capital_tax,
                                              real_estate = u_real_estate, 
                                              real_estate_returns = u_real_estate_returns, 
                                              real_estate_tax = u_real_estate_tax,
                                              access_usufruct = True,
                                              rebalance = rebalance)
    prtfls = []
    bankrupt_without_re = False
    bankrupt_with_re = False
    fire_without_re = False
    fire_with_re = False
    for year in range(int(years + 1)):
        
        if year == 0:
            prtfl = None
            infl_cmpnd = 1
        
        # prefire period:
        if age + year < f_age:
              inflation = pf_inflation
              prtfl = _get_portfolio_for_year(year = year,
                                              inflation = inflation,
                                              salary = pf_salary, 
                                              salary_increase = pf_salary_increase, 
                                              extra_income = pf_extra_income,
                                              extra_income_increase = pf_extra_income_increase,
                                              expenses = pf_expenses, 
                                              expenses_increase = pf_expenses_increase,
                                  
                                              portfolio = prtfl,
                                              capital = pf_capital, 
                                              capital_allocation = pf_capital_allocation, 
                                              capital_returns = pf_capital_returns,
                                              capital_tax = pf_capital_tax,
                                              real_estate = pf_real_estate, 
                                              real_estate_returns = pf_real_estate_returns, 
                                              real_estate_tax = pf_real_estate_tax,
                                              access_usufruct = False,
                                              rebalance = rebalance)
              if age + year == u_age:
                  prtfl = prtfl + f_prtfl_u(year, inflation)
                  
        # fire period:
        if (age + year >= f_age) & (age + year < r_age):
              inflation = f_inflation
              prtfl = _get_portfolio_for_year(year = year,
                                              inflation = inflation,
                                              salary = f_salary, 
                                              salary_increase = f_salary_increase, 
                                              extra_income = f_extra_income,
                                              extra_income_increase = f_extra_income_increase,
                                              expenses = f_expenses, 
                                              expenses_increase = f_expenses_increase,
                                  
                                              portfolio = prtfl,
                                              capital = f_capital, 
                                              capital_allocation = f_capital_allocation, 
                                              capital_returns = f_capital_returns,
                                              capital_tax = f_capital_tax,
                                              real_estate = f_real_estate, 
                                              real_estate_returns = f_real_estate_returns, 
                                              real_estate_tax = f_real_estate_tax,
                                              access_usufruct = False,
                                              rebalance = rebalance)
              if age + year == u_age:
                  prtfl = prtfl + f_prtfl_u(year, inflation)
                  
        # retirement period:
        if (age + year >= r_age):
              r_inflation = inflation
              prtfl = _get_portfolio_for_year(year = year,
                                              inflation = inflation,
                                              salary = r_salary, 
                                              salary_increase = r_salary_increase, 
                                              extra_income = r_extra_income,
                                              extra_income_increase = r_extra_income_increase,
                                              expenses = r_expenses, 
                                              expenses_increase = r_expenses_increase,
                                  
                                              portfolio = prtfl,
                                              capital = r_capital, 
                                              capital_allocation = r_capital_allocation, 
                                              capital_returns = r_capital_returns,
                                              capital_tax = r_capital_tax,
                                              real_estate = r_real_estate, 
                                              real_estate_returns = r_real_estate_returns, 
                                              real_estate_tax = r_real_estate_tax,
                                              access_usufruct = False,
                                              rebalance = rebalance)
              if age + year == u_age:
                  prtfl = prtfl + f_prtfl_u(year, inflation)
                  
        if inflation_adjusted_vals: 
            prtfl_inflation_adjusted = prtfl.get_scaled(scale = (1/((1+inflation/100)**year))) 
            prtfls.append(prtfl_inflation_adjusted)
        else:
            prtfls.append(prtfl)
            
        infl_cmpnd = ((1 + inflation/100)**year)
        fire_number_infl_adj = fire_number / infl_cmpnd
        
        
        if verbosity > 0:
            # Report status:
            if (prtfls[-1].total_with_real_estate < 0) & (bankrupt_with_re == False): 
                print('Bankrupt at age = {:1.0f} years (with real_estate)'.format(age+year))
                bankrupt_with_re = True
            if (prtfls[-1].total_without_real_estate < 0) & (bankrupt_without_re == False): 
                print('Bankrupt at age = {:1.0f} years (without real_estate)'.format(age+year))
                bankrupt_without_re = True
                
            if (prtfls[-1].total_with_real_estate >= fire_number) & (fire_with_re == False): 
                print('Fire at age = {:1.0f} years (with real_estate)'.format(age+year))
                fire_with_re = True
            if (prtfls[-1].total_without_real_estate >= fire_number) & (fire_without_re == False): 
                print('Fire at age = {:1.0f} years (without real_estate)'.format(age+year))
                fire_without_re = True

        
        # make plot:
        if verbosity > 1:
            label_stocks = 'Stocks' if (year == 0) else '_no_legend_'
            label_bonds = 'Bonds' if (year == 0) else '_no_legend_'
            label_cash = 'Cash' if (year == 0) else '_no_legend_'
            label_etfs = 'ETFs' if (year == 0) else '_no_legend_'
            label_crypto = 'Crypto' if (year == 0) else '_no_legend_'
            label_capital = 'Capital' if (year == 0) else '_no_legend_'
            label_re = 'Real estate' if (year == 0) else '_no_legend_'
            label_cap_re = 'Capital + real estate' if (year == 0) else '_no_legend_'
            label_infl_cmpnd = 'Compounded inflation' if (year == 0) else '_no_legend_'
            label_FIRE_infl_cmpnd = 'FIRE number inflation adjusted' if (year == 0) else '_no_legend_'
            ax.plot(age + year, prtfls[-1].stocks, color = 'tab:purple', marker = '>', linestyle = '--', linewidth = 1, label = label_stocks)
            ax.plot(age + year, prtfls[-1].bonds, color = 'tab:blue', marker = '<', linestyle = '--',  linewidth = 1, label = label_bonds)
            ax.plot(age + year, prtfls[-1].cash, color = 'tab:cyan', marker = 'v', linestyle = '--',  linewidth = 1, label = label_cash)
            ax.plot(age + year, prtfls[-1].etfs, color = 'tab:green', marker = '^', linestyle = '--',  linewidth = 1, label = label_etfs)
            ax.plot(age + year, prtfls[-1].crypto, color = 'tab:olive', marker = 's', linestyle = '--',  linewidth = 1, label = label_crypto)
            ax.plot(age + year, prtfls[-1].capital, color = 'tab:orange', marker = 'o', linestyle = '-',  linewidth = 2, label = label_capital)
            ax.plot(age + year, prtfls[-1].real_estate, color = 'tab:red', marker = 'd', linestyle = '-',  linewidth = 2, label = label_re)
            ax.plot(age + year, prtfls[-1].total_with_real_estate, color = 'tab:gray', marker = '.', linestyle = '-',  linewidth = 2, label = label_cap_re)
            #ax2.plot(age + year, infl_cmpnd, color = 'tab:brown', marker = 'x', linestyle = ':', linewidth = 2, label = label_infl_cmpnd)
            ax.plot(age + year, fire_number_infl_adj, color = 'k', marker = '+', linestyle = ':', linewidth = 2, label = label_FIRE_infl_cmpnd)
        
        
            
    if verbosity > 1:
        ax.hlines(0, xmin = age, xmax = end_age, color = 'k', linestyle = '-')
        ax.hlines(fire_number, xmin = age, xmax = end_age, color = 'r', linestyle = '--', label = 'SWR = {:1.2f}% -> fire at {:1.0f} euro'.format(swr,fire_number))
        if leg_dict is not None: 
            ax.legend(**leg_dict)
            #leg_dict2['bbox_to_anchor'] = (1.0, 0.0)
            #ax2.legend(**leg_dict2)

    return prtfls



if __name__ == '__main__': 
    
    run_fire_calculator_ = False
    run_fire_calculator = True 
    
    if run_fire_calculator_:
        fire_age = 47 # 45: 1457, 46: 1542, 47:1614, 48: 1684, 49: 1742 50: 1813; + 170 euro studiejaren
        retirement_income_net = (1610)*12
        capital_allocation = [60,30,10,0] # stocks, bonds, cash, etf
        capital_returns = [5,1,0,5]
        capital_extra_age = None
        inflation = [2,4,6]
        fig, ax = plt.subplots(1,len(inflation), figsize=(19,5))
        
        for i, inflation_i in enumerate(inflation):
            print('\ninflation = {:1.0f}%:'.format(inflation_i))
            ax[i].set_title('Inflation = {:1.0f}%:'.format(inflation_i))
            portfolio0 = fire_calculator_(fire_age = fire_age, 
                                         retirement_income_net = retirement_income_net,
                                        capital_allocation = capital_allocation,
                                        capital_returns = capital_returns,
                                        capital_extra_age = capital_extra_age,
                                        inflation = inflation_i, ax = ax[i])
        
        
        
    if run_fire_calculator:
        age = 45 # 45: 1457, 46: 1542, 47:1614, 48: 1684, 49: 1742 50: 1813; + 170 euro studiejaren
        f_age = 48
        u_age = None
        swr = 3.5
        rebalance = True
        inflation_adjusted_vals = True
        pf_inflation = [2,4,6,3.5]
        pf_salary = 40000
        pf_salary_increase = 1
        pf_expenses = 16000
        pf_capital = 450000
        pf_real_estate = 0
        f_expenses = 20000
        r_salary = (1610+170)*12 # + 170 for study years
        r_salary_increase = None
        r_expenses = 20000
        
        fig, ax = plt.subplots(1,len(pf_inflation), sharex = True, sharey = True, 
                               figsize = (19,5))
        
        for i,pf_inflation_i in enumerate(pf_inflation):
            print('\nInflation = {:1.2f}'.format(pf_inflation_i))
            if i < len(pf_inflation) -1: 
                leg_dict = None
            else:
                leg_dict = {'loc' : 'upper left', 'bbox_to_anchor' : (1.0,1.0), 'ncol' : 1}
            prtfls = fire_calculator(age, 
                              end_age = 85,
                              f_age = f_age,
                              r_age = 67,
                              u_age = u_age, # age at which usufruct drops 
                              
                              swr = swr,
                              rebalance = rebalance,
                              inflation_adjusted_vals = inflation_adjusted_vals,
                              
                              pf_inflation = pf_inflation_i,
                              pf_salary = pf_salary, pf_salary_increase = pf_salary_increase, 
                              pf_extra_income = 0, pf_extra_income_increase = 0,
                              pf_expenses = pf_expenses, pf_expenses_increase = 0,
                              pf_capital = pf_capital, 
                              pf_capital_allocation = [60, 30, 10, 0, 0], 
                              pf_capital_returns = [5,1,0,5,0],
                              pf_capital_tax = [0.35, 0.12, 0, 1.32 + 0.22, 0],
                              pf_real_estate = pf_real_estate, 
                              pf_real_estate_returns = 0, 
                              pf_real_estate_tax = 0,
                              
                              f_inflation = None,
                              f_salary = 0, f_salary_increase = 0, 
                              f_extra_income = 0, f_extra_income_increase = 0,
                              f_expenses = f_expenses, f_expenses_increase = 0,  
                              f_capital = 0, 
                              f_capital_allocation = None, 
                              f_capital_returns = None,
                              f_capital_tax = None,
                              f_real_estate = 0, 
                              f_real_estate_returns = None, 
                              f_real_estate_tax = None,
                              
                              r_inflation = None,
                              r_salary = r_salary, r_salary_increase = r_salary_increase, 
                              r_extra_income = 0, r_extra_income_increase = 0,
                              r_expenses = r_expenses, r_expenses_increase = 0,
                              r_capital = 0, 
                              r_capital_allocation = None, 
                              r_capital_returns = None,
                              r_capital_tax = None,
                              r_real_estate = 0, 
                              r_real_estate_returns = None, 
                              r_real_estate_tax = None,
                              
                              u_capital = 0, 
                              u_capital_allocation = None, 
                              u_capital_returns = None,
                              u_capital_tax = None,
                              u_real_estate = 0, 
                              u_real_estate_returns = None, 
                              u_real_estate_tax = None,
                              verbosity = 2, ax = ax[i],
                              leg_dict = leg_dict)
            
            ax[i].set_title('FIRE age: {:1.0f}, Inflation = {:1.2f}%'.format(f_age, pf_inflation_i))
            plt.tight_layout()
 
