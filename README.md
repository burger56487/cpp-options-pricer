# C++ Options Pricing Engine

A C++ option-pricing engine for European and American vanilla options,
implementing three independent methods with dividend yield, Greeks, input
validation and a built-in test suite.

## Methods
- **Black-Scholes** — closed-form price and Greeks (Delta, Gamma, Vega, Theta,
  Rho), with continuous dividend yield; Greeks reported in market convention
  (Vega/Rho per 1%, Theta per day)
- **CRR binomial tree** — European and American options (early-exercise premium),
  with incremental price recursion and risk-neutral probability checks
- **Monte Carlo** — risk-neutral simulation with Welford online variance and a
  95% confidence interval

## Validation (built-in test suite)
Running the program executes automated tests:
- Put-call parity
- Binomial price converges to Black-Scholes as steps increase
- Monte-Carlo price within the 95% confidence interval of Black-Scholes
- American put >= European put
- Boundary cases (expiry payoff, invalid inputs throw)

All tests pass.

## Sample Output (S=100, K=100, r=5%, q=2%, σ=20%, T=1)
| Method        | Call Price |
|---------------|------------|
| Black-Scholes | 9.2270     |
| Binomial(500) | 9.2231     |
| Monte Carlo   | 9.2907 (95% CI [9.23, 9.35]) |

## Build & Run
```bash
g++ -std=c++17 -O2 main.cpp -o options
./options
