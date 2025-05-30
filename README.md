# Open Quant Desk Project

### Quick notes

1. This is a development of a framework to promote analytical studies of financial markets. Specifically, we target options in our first release. 
2. This project is not a ready-made algo/hft platform. You will have to adapt it to your needs.
3. If you use this software, its component libraries, etc and lose money, it's your own fault.
4. By using this software, code base, or any version herein, you waive the creator, and the Open Quant Desk, Inc foundation of any liability including but not limited to financial loss.

### What is this project?

Our goal in creating this is to provide a transparent platform where developers and highly technical traders can employ advanced strategies in the Options market. We will be expanding to 
other areas, but options were our providential motivation. Within this codebase, you will find a well structured and (hopefully) performant trading system to work with a variety of
brokerage types. We specifically target Tradier in this initial release as that's who provides our Company account, there is no reason why another broker cannot be easily adapted into this application.

### Who is this project for?

In short, this project is for those who wish to use more quantative means to compute the probability of an options trade success. This is by no means an exact science, but we hope to work collaboratively and democratize this area of trading.

### Why not use Python?

Our philosophy has always been "prototype in Python, deploy in C++". There are a number of Python suites out there to do trading in any number of ways, but due to the inherent memory management complexity of the C derived languages, 
there are not many C derived publicly available trading frameworks out there. Besides that, C/C++ permit us to get intimate with the underlying hardware. In `math/core` for example, you will find specific threading optimizations geared towards doing math fast. We also 
place general program wide optimizations for CPU extenstions and GPU Extensions (currently via OpenCL) in the `optimizations/` folder.

In short, while Python is good for a lot of things, the amount of work to do in order to make it "low latency" results in rewriting it into C++ anyways, so we skip the annoyance of trying to make Python fast and just do the right thing.

# The Greeks

### First Order Greeks

First order Greeks are the most common. They are the partial derivatives of the option price, denoted as **C** for Call and **P** for put. Collectively, the Greeks describe how the option price changes linearly.

1. **Delta** is defined as the rate of change in price.
2. **Theta** is defined as the measure of rate of change with respect to **T**ime.
3. **Vega** is defined as the measure of rate of change with respect to volatility (more on volatility later)
4. **Rho** is defined as the measure of the rate of change with respect to the risk free rate (more on risk free later also!)

Together, these four values provide a picture of a given option's value and potential value. Let us consider an option chain:

This sample is sourced from the [NASDAQ Visa (V)](https://www.nasdaq.com/market-activity/stocks/v/option-chain-greeks) website. 

```
Calls 	  	  	  	  	  	  	Puts 	  	  	  	  	 
Delta	  Gamma	  Rho	      Theta	  Vega	  IV	      Strike	      Delta	  Gamma	   Rho	      Theta	  Vega	   IV
0.9793	0.00274	0.00095	-4.21095	0.00312	2.91338	    327.5	    -0.00447	0.00093	-0.85206	  0.00082	2.26183
0.98419	0.00248	0.00096	-2.9411	  0.00248	2.55163	    330	      -0.00377	0.00088	-0.66138	  0.0007	2.04494
0.9798	0.00318	0.00097	-3.48842	0.00305	2.46005	    332.5	    -0.00511	0.00121	-0.00001	-0.83021	0.00092	1.95597
0.97183	0.00428	0.00097	-4.51741	0.00404	2.41358	    335	      -0.00551	0.0014	-0.00001	-0.81813	0.00099	1.80324
```
In practice, most retail traders stop here. They pay attention only to these and are fine with it, however we can better udnerstand the market and the various strategies by looking further.

### Second Order Greeks

Sometime after the first order Greeks became a thing, someone said "Well, if one is good, two is better right?" and they took a derivative of the first order Greeks (meaning the second derivative of the price action). 
My maths teacher in school always told us "Derivatives are a curve, period", and that's what these offer. The second order Greeks allow us to analyze the curvature or convexity in the option’s behavior. 


1. **Gamma** measures the rate of change of delta with respect to the underlying price, it's calculated the same for both calls and puts. 
2. **Vanna** shows how Delta changes as volatility changes
3. **Vomma** or sometimes called **Volga** indicates how sensitive Vega is to changes in volatility
4. **Charm** is basically Delta Decay and measures the rate of change of Delta with respect to time.

Well our intrepid analysts don't stop there, we will also introduce the reader to **Third-Order** Greeks and that is where will conclude. After the second order, it becomes critical to have quantative data for models, which can be hard to reliably get. 

### Third Order Greeks

We discuss this third order only because it is helpful in constructing the surface calculations for modeling market behavior. In practice, retail traders do not need to understand or use third order Greeks because there is so much latency
and noise between our workstation and the market that the nuances these will help capture become less reliable. However, they will be included in a later release of our software after we have a complete CUDA and OpenCL framework configured.

1. **Speed** describes the rate of change of Gamma with respect to the underlying asset price.
2. **Zomma* describes the rate of change of Gamma with respect to volatility.
3. **Color** or **Gamma Decay** describes the rate of change of Gamma with respect to time.
4. **Ultima** describes the rate of change of Vomma.

So you can think of it as the rate of change for the curves from the second order.

# Volatility

If you look up the definition of what volatility is, you will read it defined as the annualized standard deviation of an underlying asset’s price returns, quantifies the expected magnitude of price fluctuations over a given period. 
In options pricing, volatility is a critical parameter in models such as Black-Scholes, as it directly influences the likelihood of an option expiring in-the-money, thereby affecting its premium. 

Most traders who use options intuitively understand that higher volatility increases the option’s value by enhancing the probability of significant price movements. 

The Greeks—first-, second-, and third-order derivatives of the option price—provide a structured framework to analyze volatility’s impact on pricing and risk. At the first-order level. **Vega** measures the sensitivity of the option price to a 1% change in implied volatility. An important note is that Vega is highest for at-the-money options, reflecting their greater sensitivity to volatility changes. At the second-order, **Vomma** quantifies the rate of change of Vega with respect to volatility by capturing the convexity of **Vega**, which is crucial for assessing how volatility shocks alter an option’s sensitivity to further volatility changes. Finally, at the third-order, **Ultima** measures the rate of change of Vomma with respect to volatility. **Ultima** provides insight into the higher-order curvature of volatility’s impact, which is particularly relevant for managing complex portfolios under volatile market conditions.


# Risk Free Rate

The risk-free rate is the theoretical return on an investment with zero risk of default, to get an idea, thing of the most stable investments imaginable ... Government Bonds and Notes. When compared with the market, the Risk Free Rate represents the baseline cost of capital in a risk-neutral world. To make this more clear, just ask "If I had bought a 3-month Treasure Note, how much money would I make versus this investment". The answer to that question is the Risk Free Rate (hereafter RFR)!

However, in finance, things are not so clear. Like every other simple concept, the nuances exist. The RFR is broken down into a few components:

1. **Time Horizon** which is simply the time until expiration. You compare like to like, so a 0.5% 3-month Treasure Bill would compare with a 90-days to expiration option.
2. **Nominal Rate** is simply the continuously compounded annualized percentage
3. **Economic Context** is the reflection of the state of the **macro**economy (the Fed's policy declarations, inflation, liquidity, etc)
4. **Risk Neutral Assumption** basically holds that the rate at the time of calculation will track for the duration of the investment

Open Quant Desk does not calibrate an assumption in our math, we calculate according the actual rates daily to re-assess our portfolio perforamnce. 






### Who are you?

Open Quant Desk, Inc is a financial software development company with pending Non-Profit status. We aim to provide a framework upon which developers can create and trade their ideas performantly. 

