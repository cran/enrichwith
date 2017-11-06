N <- 1000

## Poisson
poisson_family <- poisson()
y <- rpois(N, 100)
m <- rpois(N, 100)/N
mu <- y + rnorm(N, 0, 2)
all.equal(2 * m * (y * log(y) -  y - y * log(mu) + mu ),
          poisson_family$dev.resids(y, mu, m),
          tolerance = 1e-08)

## Normal
gaussian_family <- gaussian()
y <- rnorm(N, 10, 4)
m <- rpois(N, 100)
mu <- rnorm(N, 0, 2)

all.equal(m * (y - mu)^2,
          gaussian_family$dev.resids(y, mu, m))


## Binomial
binomial_family <- binomial()
mu <- rbeta(100, 0.9, 0.7)
z <- 1:100
m <- z + rpois(100, 20)
y <- z/m

all.equal(- 2 * m * ((y * log(mu/(1 - mu)) + log(1 - mu)) -
                 (y * log(y/(1 - y)) + log(1 - y))),
          binomial_family$dev.resids(y, mu, m))

## If y = 0 then dev.resids uses -2 * m * log(1 - mu)
all.equal(binomial_family$dev.resids(0, 0.2, 10),
          -2 * 10 * log(1 - 0.2))

## If y = m then dev.resids uses -2 * m * log(mu)
all.equal(-2 * m * log(mu),
          binomial_family$dev.resids(m/m, mu, m))


## Gamma
gamma_family <- Gamma()
mu <- rpois(N, 10)
m <- rpois(N, 100)

mu <- 14
phi <- 0.3333
nu <- 1/phi
rate <- nu/mu

## The actual variance is mu^2/nu
## mu: mean(rgamma(100000, shape = nu, rate = rate))
## var: var(rgamma(100000, shape = nu, rate = rate))

## For gamma theta = -1/mu,
all.equal(-2 * m * (log(y/mu) - (y - mu)/mu),
          Gamma()$dev.resids(y, mu, m))


