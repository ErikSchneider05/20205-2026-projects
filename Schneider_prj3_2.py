import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t


def main():
    xdata =  np.array([18,19,20,21,22,23,24,25,26,27]).reshape(-1,1)
    ydata = np.array([23.74, 23.30, 23.47, 23.19, 23.35, 23.22, 23.52, 23.32, 23.65, 23.90]).reshape(-1,1)
    
    onesvec = np.ones(xdata.shape)
    X = np.hstack([onesvec, xdata, (xdata**2)])
    Xpinv = np.linalg.pinv(X)
    beta = Xpinv @ ydata

    xplt = np.linspace(18, 27, 1000)
    yplt = beta[0] + beta[1] * xplt + beta[2] * (xplt**2)
    plt.scatter(xdata, ydata, c="k", label = "Data")
    plt.plot(xplt, yplt, c="r", label = "Regression")
    yopt = np.min(yplt)
    xopt = xplt[np.argmin(yplt)]
    plt.scatter(xopt, yopt, c="b", marker="*", label = "Optimal")
    plt.legend()
    plt.xlabel("RF pressure (psi)")
    plt.ylabel("Lap time (s)")
    plt.grid()
    print("Beta =\t", beta.flatten())
    print(f"Optimal RF pressure (xopt) = {xopt: .5f} (psi),\nExpected lap time with RF pressure {xopt: .5f} (psi) is {yopt: .5f} seconds" )
    
    ## Noise model
    ypred = X @ beta
    res = ydata - ypred
    n, p = X.shape
    sigma2 = (res.T @ res) / (n-p)
    sigma2 = sigma2[0,0] 
    print(f"Sigma^2 is  {sigma2: .5f}")

   #Standard error (SE)
    XtX_inv = np.linalg.inv(X.T @ X )
    cov_beta = sigma2 * XtX_inv
    se_beta = np.sqrt(np.diag(cov_beta))
    print("SE_beta = ", se_beta)

    #t*
    alpha = 0.05
    dof = n - p
    t_star = t.ppf(1 - alpha/2, dof)
    print(f"T* = {t_star: .5f}")
    
    #CI
    beta_hat = beta.flatten()
    beta_CI_low = beta_hat - t_star*se_beta
    beta_CI_high = beta_hat + t_star*se_beta
    print("95% CI for beta:")
    for j in range(p):
        print(f" beta_{j}: [{beta_CI_low[j]:.5f}, {beta_CI_high[j]:.5f}]")
    #ci mean
    x0 = xplt.reshape(-1, 1)
    X0 = np.hstack([np.ones(x0.shape), x0, x0**2]) 
    cov_mean = sigma2 * (X0 @ XtX_inv @ X0.T) 
    var_mean = np.diag(cov_mean)
    se_mean = np.sqrt(var_mean)
    yplt_flat = yplt.flatten()
    mean_CI_low = yplt_flat - t_star * se_mean
    mean_CI_high = yplt_flat + t_star * se_mean
    differ_of_mean = np.mean(mean_CI_high - mean_CI_low)
    print(f"Mean CI high = {np.mean(mean_CI_high):.5f}\nMean CI Low = {np.mean(mean_CI_low) :.5f}\nDifference in mean is +- {differ_of_mean: .5f}")

    ##SECONDPLOT
    plt.fill_between(xplt, mean_CI_low, mean_CI_high, alpha=0.2, label="95% CI ")

    
    plt.legend()
    plt.xlabel("RF pressure(PSI)")
    plt.ylabel("Lap times (s)")
    plt.show()

    
if __name__ == "__main__":
    
    main()