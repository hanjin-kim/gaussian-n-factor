# Additive Gaussian N factor Model (GN++)
QuantLib extension for model calibration

Followings are theoretical document of this project.
The original document is rendered by the Markdown and Latex, you'd better copy&paste the following document into your favorite Markdown&Latex viewer, at least for the moment.

# Gaussian Affine Model
```
Author : Hanjin Kim
Generated date : 2018-Jan-3rd
Last modified date : 2018-Jan-16th
Modification log :
 * (2018-May-30th, Hanjin Kim) Title changed into "Gaussian affine model"  (new reference)
 * (2018-May-30th, Hanjin Kim) Formula modified for generalized correlations
 * (2018-Jan-16th, Hanjin Kim) Section 'Swaptions' added (WIP)
 * (2018-Jan-16th, Hanjin Kim) Page format changed into MathJax + Pacdoc (w/ crossref)
 * (2018-Jan-4th, Hanjin Kim) Sections for zero bond options added
 * (2018-Jan-3rd, Hanjin Kim) Document generated
```

This document describes mathmetical details on the short-rate model which is linearly composed with N gaussian factors as,

$$
r(t) = \sum_{i=0}^{N}x_i(t) + \phi(t)
$$ {#eq:gausssolution}

where $\phi(t)$ is a deterministic term defined so as to be fitted to the initial term structure, and $x_i(t)$ has a form of the Ornstein-Uhlenbeck process :
$$
dx_i(t) = a_i(t)x(t)dt + \sigma_i(t)dW_i(t),\\
\quad x_i(0) = 0,\quad
dW_i(t)dW_j(t) = \rho_{ij}(t)dt
$$ {#eq:ousde}
where $s$ satisfies $0\le s\le t$, and $dW_i(t)$ is a Wiener process.
Here $\rho_{ij}(t)$ is a time-dependent instantaneous correlation taking values in [-1,1],
and the instantaneous correlation matrix should be non-netative semi-definite for any t.

## Solution to the Ornstein-Uhlenbeck process
To solve the @eq:ousde, first we define,
$$
E_i(s,t) = e^{-\int_s^ta(u)du}
\\\
E_i(t) = E_i(0,t)
$$

This yields a couple of properties, such as :
$$
E(u,t) = \frac{E(s,t)}{E(s,u)}\\
E(t,s) = \frac{1}{E(s,t)}
$$

The solution of the @eq:ousde can be derived by variation of parameters. Changing varable
$$
f(x(t),t) = x(t)E(s,t)
$$
we get
$$
\begin{align}
df(x,t)
 & = -xa(t)E(s,t)dt + E(s,t)dx\\
 & = -xa(t)E(s,t)dt + E(s,t)[a(t)x(t)dt + \sigma(t)dW(t)]\\
 & = \sigma(t)E(s,t)dW^Q(t)
\end{align}
$$

The solution for $f(x(t),t)$ is immediately obtained by  Itô-integrating both sides from $s$ to $t$ :
$$
f(x(t),t) = f(x(s),s) + \int_s^t \sigma(u)E(s,u)dW(u)
$$

Reversing the variation of parameters, we have :
$$
\begin{align}
x(t)
&= \frac{1}{E(s,t)}x(s) + \int_s^t \frac{\sigma(u)}{E(u,t)}dW(u)\\
\end{align}
$$ {#eq:ousolution}

## N-Gaussian factors short-rate model

### Properties of the short rate

Using the @eq:ousolution we can obtain the short rate :
$$
r(t) = \phi(t) + \sum_{i=0}^{N} \left( \frac{x_i(s)}{E_i(s,t)} + \int_s^t \frac{\sigma_i(u)}{E_i(u,t)}dW(u) \right)
$$

Its mean and variance can be combined immediately :
$$
\bar r(t) = \mathbb{E}\left[r(t)|\mathcal{F}_s\right]=\phi(t) + \sum_i \frac{x_i(s)}{E_i(s,t)}
$$ {#eq:meanofr}
and the variance is, for $s\le t\le T$ :
$$
\begin{align}
V_r(s,t)
& =  \mathbb{Var}\left[r(t)|\mathcal{F}_s\right] \\
& = \mathbb{E}\left[ \left(r(t) - \bar r(t) | \mathcal{F}_s\right)^2\right]\\
& = \mathbb{E}\left[ \left(\sum_{i=0}^{N}  \int_s^t \frac{\sigma_i(u)}{E_i(u,t)}dW_i(u)|\mathcal{F}_s \right)^2\right] \\
& = \mathbb{E}\left[
\int_s^t \left(\sum_{i=0}^{N} \frac{\sigma_i(u)}{E_i(u,t)}\right)^2du | \mathcal{F}_s
\right]\\\
& = \sum_{i,j} \int_s^t \frac{ \rho_{ij}(u)\sigma_i(u)\sigma_j(u)}{E_i(u,t)E_j(u,t)}du\\
& = \sum_{i,j} \frac{1}{E_i(t)E_j(t)} \int_s^t \rho_{ij}(u)\sigma_i(u)\sigma_j(u)E_i(u)E_j(u)du
\end{align}
$$ {#eq:varofr}

### Zero-coupon bond price

Having the @eq:meanofr, all we have to do to derived the shor-rate $r(t)$ is to obtain the deterministic solution, $\phi(t)$. And this can be done by analysing the zero-bond prices with the model. We deote by $P(t,T)$ the price at time $t$ of a zero-coupon bond maturing at $T$ and with unit face value such that :
$$
P(t,T) = \mathbb{E}^Q  \left[ e^{-\int_t^T r(u)du} | \mathcal{F}_t \right]
$$ {#eq:bondprice}
where $\mathbb{E}^Q$ denotes the expectation under the risk-adjusted measure $Q$.

Before proceed, we define a form of integral of $E(s,t)$ :
$$
\begin{align}
B_i(t,T)
& = \int_t^T\frac{du}{E_i(t,u)}=E_i(t)\int_t^T\frac{du}{E_i(u)}\\
\end{align}
$$ {#eq:BtT}

which has a few of useful properties :

$$
\begin{align}
\frac{\partial}{\partial t}B(t,T) & =a(t)B(t,T)-1 \\
\frac{\partial}{\partial T}B(t,T) & =\frac{1}{E(t,T)} \\
B(t,S) - B(t,T) & = \frac{B(T,S)}{E(t,T)}
\end{align}
$$ {#eq:BtTprop}

Define the integral of the short rate :
$$
\begin{align}
I(t,T)
& = \int_t^T r(u)du \\
& = \int_t^T\phi(u)du + \sum_{i=0}^{N} \left( \int_t^T \frac{x(s)}{E_i(s,u)}du + \int_t^T\int_s^u \frac{\sigma_i(z)}{  E_i(z,u)}dW_i(z)du \right)\\
& = \int_t^T\phi(u)du + \sum_{i=0}^{N}\frac{B_i(t,T)}{E_i(s,t)}x_i(s)  + \sum_{i=0}^{N} \int_t^T\int_s^u \frac{\sigma_i(z)}{  E_i(z,u)}dW_i(z)du
\end{align}
$$ {#eq:ItT}


Note that $I(t,T)=\int_t^T r(u)du$ is also normal [^1]. Thus, it is entirely defined by its mean and variance. For now, we denote them as $\bar I(t,T)$ and $V_I(s,t)$, respectively. The mean :

[^1]: https://quant.stackexchange.com/questions/17841/why-does-the-short-rate-in-the-hull-white-model-follow-a-normal-distribution/17863

$$
\begin{align}
\bar I(t,T)
& = \mathbb{E}\left[\int_t^T r(u)du|\mathcal{F}_t\right]\\
& = \int_t^T\phi(u)du + \sum_{i=0}^{N}\frac{B_i(t,T)}{E_i(t,t)}x_i(t)\\
& = \int_t^T\phi(u)du + \sum_{i=0}^{N} B_i(t,T)x_i(t)
\end{align}
$$

and the variance :

$$
\begin{align}
V_I(t,T)
& =  \mathbb{Var}\left[I(t,T)|\mathcal{F}_t\right] \\
& = \mathbb{E}\left[  \left( \sum_{i=0}^{N} \int_t^T \left(\int_t^u \frac{\sigma_i(s)}{E_i(s,u)}dW_i(s)\right)du | \mathcal{F}_t \right)^2\right] \\
& = \mathbb{E}\left[  \left( \sum_{i=0}^{N} \int_t^T \left(\int_s^T \frac{\sigma_i(s)}{E_i(s,u)}du\right)dW_i(s) | \mathcal{F}_t\right)^2 \right] \\
& =
\mathbb{E}\left[
\int_t^T
\left(
\sum_{i} \sigma_i(s)\int_s^T \frac{1}{E_i(s,u)}du
\right)^2ds | \mathcal{F}_t
\right]\\
& = \mathbb{E}\left[
\int_t^T
\left(
\sum_{i} \sigma_i(s) B(s,T)
\right)^2ds | \mathcal{F}_t
\right]\\
& = \sum_{i,j}\int_t^T \rho_{ij}(u)\sigma_i(u)\sigma_j(u)B_i(u,T)B_j(u,T)du
\end{align}
$$ {#eq:varItT}

Recall that if $Z$ is normal with mean $m$ and variance $V$, then $\mathbb{E}[
(Z)]=\text{exp}(m-\frac{1}{2}V)$. Likewise, since $-I(t,T)$ is normal,  we can derive the eq.7  :
$$
\begin{align}
P(t,T)
& = \mathbb{E}^Q  \left[ e^{-\int_t^T r(u)du} | \mathcal{F}_t \right]\\
& = \text{exp}\left(-\bar I(t,T) + \frac{1}{2}V_I(t,T)\right)
\end{align}
$$ {#eq:bondpriceexpanded}

Given the actual market price, $P^M(0,T)$ we can relate it with the @eq:bondpriceexpanded:
$$
\begin{align}
P^M(0,T)
& = \text{exp}\left(-\bar I(0,T) + \frac{1}{2}V_I(0,T)\right)\\
& = \text{exp}\left(-\int_0^T\phi(u)du+ \frac{1}{2}V_I(0,T)\right)
\end{align}
$$
Which leads to :
$$
\text{exp}\left(-\int_0^T\phi(u)du\right) = P^M(0,T)\text{exp}\left(-\frac{1}{2}V_I(0,T)\right)
$$

Writing the similar releation for $t\le T\le T^\star$, :
$$
\begin{align}
\text{exp}(-\int_t^T\phi(u)du)
& = \text{exp}\left(-\int_0^T\phi(u)du + -\int_0^t\phi(u)du\right)\\
& = \text{exp}\left(\text{ln}\left(\frac{P^M(0,T)}{P^M(0,t)}\right) -\frac{1}{2}(V_I(0,T)-V_I(0,t))\right)\\
\end{align}
$$

We can obtain the deterministic parameter, $\phi(t)$ :
$$
\begin{align}
\phi(t)
& = f(0,t) + \frac{1}{2}\frac{\partial V_I(0,t)}{\partial t}\\
& = f(0,t) + \frac{1}{2}\sum_{i,j}\int_0^t\rho_{ij}(u)\sigma_i(u)\sigma_j(u)(\frac{B_j(u,t)}{E_i(u,t)}+\frac{B_i(u,t)}{E_j(u,t)})du
\end{align}
$$ {#eq:phi}

where the instantaneous forward rate seen from $t$ with maturity $T$ is denoted by $f(t,T)$. With releations shown above, we can obtain the bond price in an affine form :

$$
P(t,T) = A(t,T)\text{exp}\left(-\sum_i B_i(t,T)x_i(t)\right)
$$ {#eq:affineform}

where $A(t,T)$ is defined as :

$$
\begin{align}
A(t,T) & = \frac{P^M(0,T)}{P^M(0,t)}e^{-\frac{1}{2}(V_I(0,T)-V_I(t,T)-V_I(0,t))}
\end{align}
$$

## Zero-coupon bond option

A zero-coupon bond option can be priced using the Black formula. It takes the variance of the bond price ratio to obtain the volatility. This section wil be devoted to obtain the variance of the bond price ratio.

We begin from the SDE of the zero-coupon bond price. From the eq.7 :
$$
dP(t,T) = r(t)P(t,T)dt -\sum_i \sigma_i(t)B(t,T)P(t,T)dW^Q_i(t)
$$

where the risk neutral measure is denoted by $Q$.

To calculate closedd forms, we are particularly interested in the bond ratio with fixing and paying times $T_F$ and $T_P$($t\le T_F\le T_P$), which has the dynamics in the $T_P$-forward measure :

$$
d\frac{P(t,T_F)}{P(t,T_P)}=\frac{P(t,T_F)}{P(t,T_p)}\sum_i\sigma_i(t)\left(B_i(t,T_P)-B_i(t,T_F)\right)dW_i^{T_P}(t)
$$

with integrated variance :

$$
\begin{align}
V_{P-ratio}(t,T_F,T_P)
& = \mathbb{E}\left[ \left( \sum_i\int_t^{T_F}\sigma_i(u)\left(B_i(u,T_P)-B_i(u,T_F)\right)dW_i(u) | \mathcal{F_{T_P}}\right)^2 \right]\\
& = \sum_{i,j}\int_t^{T_F}\rho_{ij}(u)\sigma_i(u)\sigma_j(u) \frac{B_i(T_F,T_P)}{E_i(u,T_F)}\frac{B_j(T_F,T_P)}{E_j(u,T_F)} du\\
& = \sum_{i,j}B_i(T_F,T_P)B_j(T_F,T_P)\int_t^{T_F} \frac{\rho_{ij}(u)\sigma_i(u)\sigma_j(u)}{E_i(u,T_F)E_j(u,T_F)}du
\end{align}
$$ {#eq:Vpratio}

For the case $N=1$, which is equivalent with the Hull-White 1 factor model, this can be reduced to :
$$
V_{P-ratio}(t,T_F,T_P) = V_r(t,T_F)B^2(T_F,T_P)
$$

Now consider an zero-bond put option with with strike $X$, fixing time $T_F$, and paying time $T_P$. It can pe priced by the Black Formula with the variance of the bond price ratio in eq.15 :
$$
\begin{align}
\text{ZBP}(T_F,T_P,X) & = XP(0,T_F)\mathcal{N}(d_+) - P(0,T_P)\mathcal{N}(d_-)\\\\
d_{\pm} & = \frac{\text{ln}\left(P(0,T_F)X/P(0,T_P)\right)}{\sqrt{V_{P-ratio}(0,T_F,T_P)}}
\pm \frac{1}{2}\sqrt{V_{P-ratio}(0,T_F,T_P)}
\end{align}
$$ {#eq:zerobondputprice}

## Swaptions

Estimating the price of swaptions is essential especially when the N-factor Gaussian model is calibrated to eg. the market data. For the 1-factor case, the well-known __Jamshidian's Trick__ (@Jamshidian). In this section we try to obtain a closed from of the swaption price for the general N-factor model.

### T-Forward Dynamics

Conversion of the drift of a process $X$ under the numerarire $U$ into another numerarire $S$ can be done by a relation :
$$
\mu_T^U(X_t) = \mu_T^S(X_t) - \sigma_t(X_t)\rho\left(\frac{\sigma_t^S}{S_t}-\frac{\sigma_t^U}{U_t}\right)^\intercal
$$ {#eq:MeasureConversion}

with the numeraries evolving under $Q^U$ according to :
$$
\begin{align}
dS_t & = (...)dt + \sigma_t^SCdW_t^U,\quad Q^U\\
dU_t & = (...)dt + \sigma_t^UCdW_t^U,\quad Q^U
\end{align}
$$

where $\sigma_t^{U/S}$ are $1\times N$ vectors, $W^U$ is the usual *n-dimensional drifless (under $Q^U$) standard Brownian motion* and $CC^\intercal=\rho$ refering to the correlation matrix within the process $X$. Its proof can be found in the sectuin 2.3 of "@Brigo:2006".

Since @eq:ousde is of the bank-acount measure, its T-forward dynamics can be obtained with having $U(t) = E(0,t)$ and $S(t) = P(t,T)$. Using the @eq:affineform, we have :
$$
\frac{\partial x_i(t)}{\partial t} = -a_i(t)x_i(t) - \sum_j\rho_{ij}\sigma_i(t)\sigma_j(t)B_j(t,T)
$$ {#eq:ouTforwarddrift}

Then the T-forward dynamics of a gaussian factor becomes :
$$
dx_i(t) = \left[-a_i(t)x_i(t) + \sum_j\rho_{ij}\sigma_i(t)\sigma_j(t)B_j(t,T)\right]dt
          + \sigma dW_i^T(t)
$$ {#eq:ouTforwardsde}

which leads to the T-forward solutions :
$$
x_i(t) = \frac{1}{E(s,t)}x(s) - M_i^T(s,t) + \int_s^t \frac{\sigma_i(u)}{E_i(u,t)}dW_i^T(u)
$$ {#eq:ouTforwardSolution}

with
(WIP)
$$
M_i^T(s,t) = \sum_j\int_s^t\rho_{ij}\sigma_i(u)\sigma_j(u)\frac{B_j(u,T)}{E_i(u,t)}du
$$

# Efficient Simulation Schema for Zero Bond Pricing


# References
