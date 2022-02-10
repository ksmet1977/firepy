# -*- coding: utf-8 -*-
"""
FIRE calculator streamlit app
-----------------------------

Created on Wed Feb  9 11:32:27 2022

@author: KSmet
"""
from firepy import fire_calculator
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# Streamlit app
#------------------------------------------------------------------------------

import streamlit as st

# App title:
st.set_page_config(page_title="FIRE Calculator")
st.title("FIRE Calculator")

st.header("**Investments**")
col_inv_val, col_re_val = st.columns(2)
with col_inv_val:
    investment_val = st.number_input("Initial capital (euro): ", value = 500000.0, min_value = 0.0, step = 10000.0, format = '%f', key = 'invest_value')
with col_re_val: 
    real_estate_val = st.number_input("Real estate value (euro): ", value = 0.0, min_value = 0.0, step = 10000.0, format = '%f', key = 'real_estate_val') 

col_inv_allo, col_inv_ret, col_inv_tax = st.columns(3)  
with col_inv_allo:
    st.subheader('Allocation')
    #exp_allo = st.expander('Investment Allocation')
    #with exp_allo:
    stocks_allo = st.number_input("Stock (acc.) allocation (%): ", value = 60.0, min_value = 0.0, step = 1.0, format = '%f', key = 'stock_allo')
    bonds_allo = st.number_input("Bond allocation (%): ", value = 30.0, min_value = 0.0, step = 1.0, format = '%f', key = 'bond_allo')
    cash_allo = st.number_input("Cash allocation (%): ", value = 10.0, min_value = 0.0, step = 1.0, format = '%f', key = 'cash_allo')
    etfs_allo = st.number_input("ETF allocation (%): ", value = 0.0, min_value = 0.0, step = 1.0, format = '%f', key = 'etf_allo')
    crypto_allo = st.number_input("Crypto allocation (%): ", value = 0.0, min_value = 0.0, step = 1.0, format = '%f', key = 'crypto_allo')
             
with col_inv_ret:
    st.subheader('Returns')
    #exp_ret = st.expander('Investment Returns')
    #with exp_ret:
    stocks_ret = st.number_input("Stock (acc.) return (%): ", value = 5.0, min_value = 0.0, step = 1.0, format = '%f', key = 'stock_ret')
    bonds_ret = st.number_input("Bond return (%): ", value = 1.0, min_value = 0.0, step = 1.0, format = '%f', key = 'bond_ret')
    cash_ret = st.number_input("Cash return (%): ", value = 0.0, min_value = 0.0, step = 1.0, format = '%f', key = 'cash_ret')
    etfs_ret = st.number_input("ETF return (%): ", value = 0.0, min_value = 0.0, step = 1.0, format = '%f', key = 'etf_ret')
    crypto_ret = st.number_input("Crypto return (%): ", value = 0.0, min_value = 0.0, step = 1.0, format = '%f', key = 'crypto_ret')
    real_estate_ret = st.number_input("Real estate return (%): ", value = 0.0, min_value = 0.0, step = 1.0, format = '%f', key = 'real_estate_ret') 

with col_inv_tax:
    # enkel beurstaxen op aandelen, obligaties en funds (etfs), veronderstel geen uitkering van dividenden en obligatiecoupons, alsook verwaarloosbare onroerende voorheffing op bouwgrond
    st.subheader('Stock market tax')
    #exp_tax = st.expander('Stock market tax')
    #with exp_tax:
    stocks_tax = st.number_input("Stock (acc.) tax (%): ", value = 0.35, min_value = 0.0, step = 0.1, format = '%f', key = 'stock_tax')
    bonds_tax = st.number_input("Bond tax (%): ", value = 0.12, min_value = 0.0, step = 0.1, format = '%f', key = 'bond_tax')
    cash_tax = st.number_input("Cash tax (%): ", value = 0.0, min_value = 0.0, step = 0.1, format = '%f', key = 'cash_tax')
    etfs_tax = st.number_input("ETF tax (%): ", value = 1.32, min_value = 0.0, step = 0.1, format = '%f', key = 'etf_tax')
    crypto_tax = st.number_input("Crypto tax (%): ", value = 0.0, min_value = 0.0, step = 0.1, format = '%f', key = 'crypto_tax')
    real_estate_tax = st.number_input("Real estate tax (%): ", value = 0.0, min_value = 0.0, step = 0.1, format = '%f', key = 'real_estate_tax') 


st.header("**Usufruct investments**")
col_u_inv_val, col_u_re_val,col_u_age = st.columns(3)
with col_u_inv_val: 
    u_investment_val = st.number_input("Initial capital (euro): ", value = 0.0, min_value = 0.0, step = 10000.0, format = '%f', key = 'u_invest_value')
with col_u_re_val: 
    u_real_estate_val = st.number_input("Real estate value (euro): ", value = 0.0, min_value = 0.0, step = 10000.0, format = '%f', key = 'u_real_estate_val') 
with col_u_age: 
    u_age = st.number_input("Usufruct stop age (years) ", value = 1000.0, min_value = 0.0, step = 1.0, format = '%f', key = 'u_age')


st.header('**Income and expenses**')
col_pf, col_f, col_ret = st.columns(3)

with col_pf:
    st.subheader("Pre-FIRE period")
    #exp_pf = st.expander("Pre-FIRE period")
    #with exp_pf:
    age = st.number_input("Current age (years): ", value = 45.0,  min_value = 0.0, step = 1.0, format = '%f', key = 'pf_age')
    pf_inflation = st.number_input("Inflation (%): ", value = 3.5, min_value = 0.0, step = 0.5, format = '%f', key = 'pf_infl')
    pf_salary = st.number_input("Yearly net salary (euro): ", min_value = 0.0, step = 1000.0, format = '%f', key = 'pf_sal')
    pf_salary_increase = st.number_input("Yearly salary increase (%): ", step = 1.0, format = '%f', key = 'pf_sal_inc')
    pf_extra_income = st.number_input("Yearly extra income (euro): ", min_value = 0.0, step = 1000.0, format = '%f', key = 'pf_exti')
    pf_extra_income_increase = st.number_input("Yearly extra income increase (%): ", step = 1.0, format = '%f', key = 'pf_exti_inc')
    pf_expenses = st.number_input("Yearly expenses (euro): ", min_value = 0.0, step = 1000.0, format = '%f', key = 'pf_exp')
    pf_expenses_increase = st.number_input("Yearly expense increase (%): ", step = 1.0, format = '%f', key = 'pf_exp_inc')

with col_f:
    st.subheader("FIRE period")
    #exp_f = st.expander("FIRE period")
    #with exp_f:
    f_age = st.number_input("FIRE age (years): ", value = 50.0, min_value = 0.0, step = 1.0, format = '%f', key = 'f_age')
    f_inflation = st.number_input("Inflation (%): ", value = 3.5, min_value = 0.0, step = 0.5, format = '%f', key = 'f_infl')
    f_salary = st.number_input("Yearly net salary (euro): ", value = 0.0, min_value = 0.0, step = 1000.0, format = '%f', key = 'f_sal')
    f_salary_increase = st.number_input("Yearly salary increase (%): ", value = 0.0, step = 1.0, format = '%f', key = 'f_sal_inc')
    f_extra_income = st.number_input("Yearly extra income (euro): ", value = 0.0, min_value = 0.0, step = 1000.0, format = '%f', key = 'f_exti')
    f_extra_income_increase = st.number_input("Yearly extra income increase (%): ", value = 0.0, step = 1.0, format = '%f', key = 'f_exti_inc')
    f_expenses = st.number_input("Yearly expenses (euro): ", min_value = 0.0, step = 1000.0, format = '%f', key = 'f_exp')
    f_expenses_increase = st.number_input("Yearly expense increase (%): ", step = 1.0, format = '%f', key = 'f_exp_inc')
    
    
with col_ret:
    st.subheader("Retirement period")
    #exp_r = st.expander("Retirement period")
    #with exp_r:
    r_age = st.number_input("Retirement age (years) ", value = 67.0, min_value = 0.0, step = 1.0, format = '%f', key = 'ret_age')
    r_inflation = st.number_input("Inflation (%): ", value = 3.5, min_value = 0.0, step = 0.5, format = '%f', key = 'ret_infl')
    r_salary = st.number_input("Yearly net salary (euro): ", value = 0.0, min_value = 0.0, step = 1000.0, format = '%f', key = 'r_sal')
    r_salary_increase = st.number_input("Yearly salary increase (%): ", value = 0.0, step = 1.0, format = '%f', key = 'r_sal_inc')
    r_extra_income = st.number_input("Yearly extra income (euro): ", min_value = 0.0, step = 1000.0, format = '%f', key = 'ret_income')
    r_extra_income_increase = st.number_input("Yearly extra income increase (%): ", step = 1.0, format = '%f', key = 'ret_income_inc')
    r_expenses = st.number_input("Yearly expenses (euro): ", min_value = 0.0, step = 1000.0, format = '%f', key = 'r_exp')
    r_expenses_increase = st.number_input("Yearly expense increase (%): ", step = 1.0, format = '%f', key = 'r_exp_inc')

st.header("Misc. options")
col_swr, col_end_age, col_infl_adj = st.columns(3)  
with col_swr:
    swr = st.number_input("Safe Withdrawal Rate (%): ", value = 3.5, min_value = 0.0, step = 0.1, format = '%f', key = 'swr')
with col_end_age: 
    end_age = st.number_input("Max. age (years) to calculate: ", value = 85.0, min_value = 0.0, step = 1.0, format = '%f', key = 'last_age')
with col_infl_adj:
    inflation_adjusted_vals = st.checkbox('Plot portfolio values adjusted for inflation  (=current euros)', value = False)
        
if st.button('Calculate', key='calc', help="Run FIRE calculation"):    
    fig, ax = plt.subplots(1,1, figsize = (8,8))
    leg_dict = leg_dict = {'loc' : 'upper left', 'bbox_to_anchor' : (1.0,1.0), 'ncol' : 1}
    prtfls = fire_calculator(age, 
                        end_age = end_age,
                        f_age = f_age,
                        r_age = r_age,
                        u_age = u_age, # age at which usufruct drops 
                        
                        swr = swr,
                        rebalance = True,
                        inflation_adjusted_vals = inflation_adjusted_vals,
                        
                        pf_inflation = pf_inflation,
                        pf_salary = pf_salary, pf_salary_increase = pf_salary_increase, 
                        pf_extra_income = pf_extra_income, pf_extra_income_increase = pf_extra_income_increase,
                        pf_expenses = pf_expenses, pf_expenses_increase = pf_expenses_increase,
                        pf_capital = investment_val, 
                        pf_capital_allocation = [stocks_allo, bonds_allo, cash_allo, etfs_allo, crypto_allo], 
                        pf_capital_returns = [stocks_ret, bonds_ret, cash_ret, etfs_ret, crypto_ret],
                        pf_capital_tax = [stocks_tax, bonds_tax, cash_tax, etfs_tax, crypto_tax],
                        pf_real_estate = real_estate_val, 
                        pf_real_estate_returns = real_estate_ret, 
                        pf_real_estate_tax = real_estate_tax,
                        
                        f_inflation = f_inflation,
                        f_salary = f_salary, f_salary_increase = f_salary_increase, 
                        f_extra_income = f_extra_income_increase, f_extra_income_increase = f_extra_income_increase,
                        f_expenses = f_expenses, f_expenses_increase = f_expenses_increase,  
                        f_capital = 0, 
                        f_capital_allocation = [stocks_allo, bonds_allo, cash_allo, etfs_allo, crypto_allo], 
                        f_capital_returns = [stocks_ret, bonds_ret, cash_ret, etfs_ret, crypto_ret],
                        f_capital_tax = [stocks_tax, bonds_tax, cash_tax, etfs_tax, crypto_tax],
                        f_real_estate = 0, 
                        f_real_estate_returns = real_estate_ret, 
                        f_real_estate_tax = real_estate_tax,
                        
                        r_inflation = r_inflation,
                        r_salary = r_salary, r_salary_increase = r_salary_increase, 
                        r_extra_income = r_extra_income, r_extra_income_increase = r_extra_income_increase,
                        r_expenses = r_expenses, r_expenses_increase = r_expenses_increase,
                        r_capital = 0, 
                        r_capital_allocation = [stocks_allo, bonds_allo, cash_allo, etfs_allo, crypto_allo], 
                        r_capital_returns = [stocks_ret, bonds_ret, cash_ret, etfs_ret, crypto_ret],
                        r_capital_tax = [stocks_tax, bonds_tax, cash_tax, etfs_tax, crypto_tax],
                        r_real_estate = 0, 
                        r_real_estate_returns = real_estate_ret, 
                        r_real_estate_tax = real_estate_tax,
                        
                        u_capital = u_investment_val, 
                        u_capital_allocation = [stocks_allo, bonds_allo, cash_allo, etfs_allo, crypto_allo], 
                        u_capital_returns = [stocks_ret, bonds_ret, cash_ret, etfs_ret, crypto_ret],
                        u_capital_tax = [stocks_tax, bonds_tax, cash_tax, etfs_tax, crypto_tax],
                        u_real_estate = u_real_estate_val, 
                        u_real_estate_returns = real_estate_ret, 
                        u_real_estate_tax = real_estate_tax,
                        verbosity = 2, ax = ax,
                        leg_dict = leg_dict)

    st.pyplot(fig)
 
