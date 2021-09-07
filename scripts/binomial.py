from scipy.stats import binom
import numpy


def main():
    n = 4096
    p = 1/16411
    k = 10000       # ???

    distribution = binom(n=n, p=p)

    x = [0,1,2,3,4,5]

    cdf = distribution.cdf(x)

    print("analytical binomial CDF")
    for i,probability in enumerate(cdf):
        print(i, probability)

    observed = [12795, 3172, 411, 30, 3]

    print("observed CDF")
    cumulative_probability = 0
    for i,count in enumerate(observed):
        cumulative_probability += count/sum(observed)
        print(i,cumulative_probability)



if __name__ == "__main__":
    main()
