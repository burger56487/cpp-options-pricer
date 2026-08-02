// main.cpp — Option pricing engine: Black-Scholes (+Greeks, dividend yield),
// CRR binomial (European/American), Monte Carlo (Welford + 95% CI),
// with input validation and a built-in test suite.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <string>

const double PI = 3.14159265358979323846;

enum class OptionType { Call, Put };
enum class Exercise   { European, American };

struct OptionParams {
    double S, K, r, q, sigma, T;   // spot, strike, rate, dividend yield, vol, maturity
    OptionType type;
};
struct Greeks { double price, delta, gamma, vega, theta, rho; };
struct MCResult { double price, std_err, ci_low, ci_high; };

double norm_pdf(double x){ return std::exp(-0.5*x*x)/std::sqrt(2.0*PI); }
double norm_cdf(double x){ return 0.5*std::erfc(-x/std::sqrt(2.0)); }

void validate(const OptionParams& p){
    if (p.S<=0 || p.K<=0) throw std::invalid_argument("S and K must be positive");
    if (p.sigma<0)        throw std::invalid_argument("sigma must be non-negative");
    if (p.T<0)            throw std::invalid_argument("T must be non-negative");
}

// ================= Black-Scholes (with dividend yield) =================
Greeks black_scholes(const OptionParams& p){
    validate(p);
    bool call = (p.type==OptionType::Call);
    // deterministic boundary: expiry or zero vol
    if (p.T==0.0 || p.sigma==0.0){
        double fwd  = p.S*std::exp((p.r-p.q)*p.T);
        double disc = std::exp(-p.r*p.T);
        double intr = call ? std::max(fwd-p.K,0.0) : std::max(p.K-fwd,0.0);
        Greeks g{}; g.price = disc*intr; return g;
    }
    double sqrtT = std::sqrt(p.T);
    double d1 = (std::log(p.S/p.K)+(p.r-p.q+0.5*p.sigma*p.sigma)*p.T)/(p.sigma*sqrtT);
    double d2 = d1 - p.sigma*sqrtT;
    double dq = std::exp(-p.q*p.T), dr = std::exp(-p.r*p.T);
    Greeks g{};
    if (call){
        g.price = p.S*dq*norm_cdf(d1) - p.K*dr*norm_cdf(d2);
        g.delta = dq*norm_cdf(d1);
        g.theta = -(p.S*dq*norm_pdf(d1)*p.sigma)/(2*sqrtT)
                  - p.r*p.K*dr*norm_cdf(d2) + p.q*p.S*dq*norm_cdf(d1);
        g.rho   = p.K*p.T*dr*norm_cdf(d2);
    } else {
        g.price = p.K*dr*norm_cdf(-d2) - p.S*dq*norm_cdf(-d1);
        g.delta = -dq*norm_cdf(-d1);
        g.theta = -(p.S*dq*norm_pdf(d1)*p.sigma)/(2*sqrtT)
                  + p.r*p.K*dr*norm_cdf(-d2) - p.q*p.S*dq*norm_cdf(-d1);
        g.rho   = -p.K*p.T*dr*norm_cdf(-d2);
    }
    g.gamma = dq*norm_pdf(d1)/(p.S*p.sigma*sqrtT);
    g.vega  = p.S*dq*norm_pdf(d1)*sqrtT;
    return g;
}

// ================= CRR binomial tree (with q) =================
double binomial(const OptionParams& p, Exercise ex, int steps=500){
    validate(p);
    if (steps<=0) throw std::invalid_argument("steps must be positive");
    bool call = (p.type==OptionType::Call);
    if (p.T==0.0) return call?std::max(p.S-p.K,0.0):std::max(p.K-p.S,0.0);

    double dt   = p.T/steps;
    double u    = std::exp(p.sigma*std::sqrt(dt)), d = 1.0/u;
    double disc = std::exp(-p.r*dt);
    double prob = (std::exp((p.r-p.q)*dt)-d)/(u-d);
    if (prob<0.0 || prob>1.0)
        throw std::runtime_error("risk-neutral probability out of [0,1] — check inputs/steps");

    std::vector<double> val(steps+1);
    double price = p.S*std::pow(d, steps);           // terminal, 0 up-moves
    for (int j=0;j<=steps;++j){
        val[j] = call?std::max(price-p.K,0.0):std::max(p.K-price,0.0);
        price *= u/d;
    }
    for (int step=steps-1; step>=0; --step){
        double sprice = p.S*std::pow(d, step);       // node price at 0 up-moves (one pow per level)
        for (int j=0;j<=step;++j){
            val[j] = disc*(prob*val[j+1] + (1-prob)*val[j]);
            if (ex==Exercise::American){
                double exval = call?std::max(sprice-p.K,0.0):std::max(p.K-sprice,0.0);
                val[j] = std::max(val[j], exval);
                sprice *= u/d;
            }
        }
    }
    return val[0];
}

// ================= Monte Carlo (Welford + 95% CI) =================
MCResult monte_carlo(const OptionParams& p, long n_paths=200000, unsigned seed=0){
    validate(p);
    if (n_paths<=0) throw std::invalid_argument("n_paths must be positive");
    bool call = (p.type==OptionType::Call);
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::normal_distribution<double> Z(0.0,1.0);
    double drift = (p.r-p.q-0.5*p.sigma*p.sigma)*p.T;
    double vol   = p.sigma*std::sqrt(p.T);
    double disc  = std::exp(-p.r*p.T);

    double mean=0.0, M2=0.0; long n=0;            // Welford online variance
    for (long i=0;i<n_paths;++i){
        double ST = p.S*std::exp(drift+vol*Z(rng));
        double pv = disc*(call?std::max(ST-p.K,0.0):std::max(p.K-ST,0.0));
        ++n; double delta=pv-mean; mean+=delta/n; M2+=delta*(pv-mean);
    }
    double var = (n>1)?M2/(n-1):0.0;
    double se  = std::sqrt(var/n);
    return {mean, se, mean-1.96*se, mean+1.96*se};
}

// ================= Test suite =================
int passed=0, failed=0;
void check(const std::string& name, bool cond){
    (cond?++passed:++failed);
    std::cout << (cond?"  [PASS] ":"  [FAIL] ") << name << "\n";
}
void run_tests(){
    std::cout << "Running tests...\n";
    OptionParams call{100,100,0.05,0.0,0.2,1.0,OptionType::Call};
    OptionParams put {100,100,0.05,0.0,0.2,1.0,OptionType::Put};

    double C=black_scholes(call).price, P=black_scholes(put).price;
    double parity=call.S*std::exp(-call.q*call.T)-call.K*std::exp(-call.r*call.T);
    check("Put-call parity", std::abs((C-P)-parity)<1e-8);
    check("Binomial converges to BS", std::abs(binomial(call,Exercise::European,2000)-C)<0.05);
    MCResult mc=monte_carlo(call,400000,12345);
    check("MC within 95% CI of BS", mc.ci_low<=C && C<=mc.ci_high);
    check("American put >= European put",
          binomial(put,Exercise::American,500) >= binomial(put,Exercise::European,500)-1e-9);
    OptionParams itm{120,100,0.05,0.0,0.2,0.0,OptionType::Call};
    check("T=0 call intrinsic == 20", std::abs(black_scholes(itm).price-20.0)<1e-9);
    bool threw=false;
    try { black_scholes({-1,100,0.05,0.0,0.2,1.0,OptionType::Call}); }
    catch(const std::invalid_argument&){ threw=true; }
    check("Negative S throws", threw);

    std::cout << "Tests: " << passed << " passed, " << failed << " failed.\n\n";
}

// ================= Demo =================
int main(){
    run_tests();

    OptionParams call{100,100,0.05,0.02,0.2,1.0,OptionType::Call};  // q = 2%
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "European Call (S=100,K=100,r=5%,q=2%,sigma=20%,T=1)\n";
    std::cout << "===================================================\n";
    Greeks g = black_scholes(call);
    std::cout << "Black-Scholes : " << g.price << "\n";
    std::cout << "Binomial(500) : " << binomial(call,Exercise::European) << "\n";
    MCResult mc = monte_carlo(call,200000,42);
    std::cout << "Monte Carlo   : " << mc.price
              << "  95% CI [" << mc.ci_low << ", " << mc.ci_high << "]\n\n";

    std::cout << "Greeks (market convention):\n";
    std::cout << "  Delta           : " << g.delta << "\n";
    std::cout << "  Gamma           : " << g.gamma << "\n";
    std::cout << "  Vega (per 1%)   : " << g.vega/100.0 << "\n";
    std::cout << "  Theta (per day) : " << g.theta/365.0 << "\n";
    std::cout << "  Rho  (per 1%)   : " << g.rho/100.0 << "\n\n";

    OptionParams put{100,100,0.05,0.02,0.2,1.0,OptionType::Put};
    double eup=binomial(put,Exercise::European), amp=binomial(put,Exercise::American);
    std::cout << "American vs European put: " << amp << " vs " << eup
              << "  (early-exercise premium " << amp-eup << ")\n";
    return 0;
}
