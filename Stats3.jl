"""
    Fn(x, xobs)

Empirical distribution function evaluated at value `x` using the observed random sample in a vector `xobs`.

# Example
```
Fn(1/2, rand(50))
```
"""
Fn(x, xobs) = sum(xobs .≤ x) / length(xobs)

"""
    Tn(interval::String, xobs)

Nonparametric point estimation for the probability of an `interval` using the observed random sample in a vector `xobs`.
# Example
```
Tn("]0.25,0.75]", rand(50))
```
"""
function Tn(interval::String, xobs)
    m = length(interval)
    brackets = ['[', ']']
    if !issubset([interval[1], interval[m]], brackets)
        print("Error: interval must start and end with brackets.")
        return nothing
    end
    icomma = 0
    for i ∈ 1:m
        if interval[i] == ','
            icomma = i
        end
    end
    if icomma == 0
        print("Error: interval |a,b| extremes must be separated by a comma.")
        return nothing
    end
    a = parse(Float64, interval[2:(icomma - 1)])
    b = parse(Float64, interval[(icomma + 1):(m-1)])
    if a > b
        print("Error: interval |a,b| extremes must satisfy a ≤ b.")
        return nothing
    end
    n = length(xobs)
    if interval[1] == ']' && interval[m] == ']'
        tn = sum(a .< xobs .≤ b) / n
    end
    if interval[1] == ']' && interval[m] == '['
        tn = sum(a .< xobs .< b) / n
    end
    if interval[1] == '[' && interval[m] == ']'
        tn = sum(a .≤ xobs .≤ b) / n
    end
    if interval[1] == '[' && interval[m] == '['
        tn = sum(a .≤ xobs .< b) / n
    end
    return tn
end

"""
    EDA(fobj, valmin, valmax; iEnteros = zeros(Int, 0), tamgen = 1000, propselec = 0.3, difmax = 0.00001, maxiter = 1000)

`fobj` A real function of several variables to be minimized, where its argument is a vector or 1-D array.

`valmin, valmax` vectors or 1-D arrays of minimum and maximum values for the first generation.

`iEnteros` index of variables that must take integer values.

`tamgen` size of each generation (1000 by default, if not specified).

`propselec` proportion of population to be selected (0.3 by default, if not specified).

`difmax` error tolerance (0.00001 by default, if not specified).

`maxiter` maximum number of iterations (1000 by default, if not specified).

# Example 1
```
f(x) = (x[1] - 5)^4 - 16(x[1] - 5)^2 + 5(x[1] - 5) + 120
EDA(f, [0], [9])
```

# Example 2
This non-negative function is clearly minimized at (5,-2).
```
f(z) = abs(z[1] - 5) + abs(z[2] + 2)
EDA(f, [-10, -10], [10, 10])
```

# Example 3
The same function but only allowing integer values:
```
EDA(f, [-10, -10], [10, 10], iEnteros = [1, 2])
```
"""
function EDA(fobj, valmin, valmax; iEnteros = zeros(Int, 0), tamgen = 1000,
             propselec = 0.3, difmax = 0.00001, maxiter = 1000)
    numiter = 1
    println("Iterando... ")
    numvar = length(valmin)
    nselec = Int(round(tamgen * propselec))
    G = zeros(tamgen, numvar)
    Gselec = zeros(nselec, numvar)
    for j ∈ 1:numvar
        G[:, j] = valmin[j] .+ (valmax[j] - valmin[j]) .* rand(tamgen)
    end
    if length(iEnteros) > 0
        for j ∈ iEnteros
            G[:, j] = round.(G[:, j])
        end
    end
    d(x, y) = sqrt(sum((x .- y) .^ 2))
    rnorm(n, μ, σ) = μ .+ (σ .* randn(n))
    promedio(x) = sum(x) / length(x)
    desvest(x) = sqrt(sum((x .- promedio(x)) .^ 2) / (length(x) - 1))
    fG = zeros(tamgen)
    maxGselec = zeros(tamgen)
    minGselec = zeros(tamgen)
    media = zeros(numvar)
    desv = zeros(numvar)
    while numiter < maxiter
        # evaluando función objetivo en generación actual:
        print(numiter, "\r")
        for i ∈ 1:tamgen
            fG[i] = fobj(G[i, :])
        end
        # seleccionando de generación actual:
        umbral = sort(fG)[nselec]
        iSelec = findall(fG .≤ umbral)
        Gselec = G[iSelec, :]
        for j ∈ 1:numvar
            maxGselec[j] = maximum(Gselec[:, j])
            minGselec[j] = minimum(Gselec[:, j])
            media[j] = promedio(Gselec[:, j])
            desv[j] = desvest(Gselec[:, j])
        end
        # salir del ciclo si se cumple criterio de paro:
        if d(minGselec, maxGselec) < difmax 
            break
        end
        # y si no se cumple criterio de paro, nueva generación:
        numiter += 1
        for j ∈ 1:numvar
            G[:, j] = rnorm(tamgen, media[j], desv[j])
        end
        if length(iEnteros) > 0
            for j ∈ iEnteros
                G[:, j] = round.(G[:, j])
            end
        end
    end
    println("...fin")
    fGselec = zeros(nselec)
    for i ∈ 1:length(fGselec)
        fGselec[i] = fobj(Gselec[i, :])
    end
    xopt = Gselec[findmin(fGselec)[2], :]
    if length(iEnteros) > 0
        for j ∈ iEnteros
            xopt[j] = round(xopt[j])
        end
    end
    fxopt = fobj(xopt)
    r = (x = xopt, fx = fxopt, iter = numiter)
    if numiter == maxiter
        println("Aviso: se alcanzó el máximo número de iteraciones = ", maxiter)
    end
    return r
end

using Distributions

"""
    Bn(interval::String, obs, g = 0.95)

Nonparametric point and 100g% probability interval estimation for the probability of an `interval` using the observed random sample in a vector `obs` and the bayesian paradigm (g = 0.95 default value if not specified). `interval` is specified a string, where the first and last characters must be an open or closed bracket, that is `]` or `[` and the left and right extremes must be numbers separated by a comma.
# Example
```
W = Normal(-2, 3)
cdf(W, 3) - cdf(W, 0) # P(0 < W < 3) 
wobs = rand(W, 100)
Bn("]0,3[", wobs, 0.95)
```
"""
function Bn(interval::String, obs, γ = 0.95)
    # using: Distributions
    # Dependencies: Tn, EDA
    n = length(obs)
    tn = Tn(interval, obs)
    α, β = 1 + n*tn, 1 + n*(1 - tn) # posterior parameters
    Θ = Beta(α, β) # posterior distribution
    θmedia, θmediana = mean(Θ), median(Θ)
    h(z) = (quantile(Θ, γ + cdf(Θ, z[1])) - z[1]) * Inf^(z[1] > quantile(Θ, 1 - γ))
    sol = EDA(h, [0], [quantile(Θ, 1- γ)])
    θ₁ = sol[1][1]
    θ₂ = quantile(Θ, γ + cdf(Θ, θ₁))
    estimación = (insesgado = tn, media = θmedia, mediana = θmediana, intervalo = (θ₁, θ₂))
    return estimación
end

"""
    cuantil_puntual(α, obs)

Empirical `α` quantile point estimation using the observed random sample in a vector `obs`.

# Example
```
cuantil_puntual(0.97725, randn(100000)) # theoretically it's 2.0 approx
```
"""
function cuantil_puntual(α, obs)
    sort!(obs)
    n = length(obs)
    obsmin = minimum(obs)
    obsmax = maximum(obs)
    orden = (n+1)*α
    if orden < 1
        cuantil = obsmin
        println("Aviso: cuantil fuera de rango muestral, es menor.")
        return cuantil
    end
    if orden > n
        cuantil = obsmax
        println("Aviso: cuantil fuera de rango muestral, es mayor.")
        return cuantil
    end
    if orden == round(orden)
        j = Int(orden)
        cuantil = obs[j]
        return cuantil
    end
    # interpolar:
    j = Int(floor(orden))
    cuantil = (j + 1  - orden)*obs[j] + (orden - j)*obs[j+1]
    return cuantil
end

"""
    GoF(xobs, Fo, simFo; prueba = "AD", numsims = 100_000)

Goodness of fit test for an observed random sample `xobs` (as a vector), a proposed distribution function `Fo` along with a function `simFo` which simulates samples from `Fo` with a given size. `prueba` is a string indicating `AD` for Anderson-Darling (default), `CM` for Cramér-von Mises, and `KS` for Kolmogorov-Smirnov. `numsims` is the number of simulations to aproximate the probability distribution of the chosen test statistic (100,000 by default).

# Example
```
xobs = rand(100) # random sample from a continuous Uniform(0,1) distribution
b = 1.1
Fo(x) = (0 < x < b) * x / b + 1*(x ≥ b) # Uniform(0,b) distribution function
simFo(n) = b * rand(n) # simulates a size n random sample from a Uniform(0,b)
GoF(xobs, Fo, simFo; prueba = "AD")
```
"""
function GoF(xobs, Fo, simFo; prueba = "AD", numsims = 100_000)
    n = length(xobs)
    u(x) = Fo.(sort(x))
    i1 = collect(1:n)
    i0 = collect(0:(n-1))
    KS(x) = max(maximum(i1/n .- u(x)), maximum(u(x) .- i0/n))
    CM(x) = 1/(12*n) + sum((u(x) .- (2 .* i1 .- 1) ./ (2*n)).^2)
    AD(x) = -n - (1/n)*sum((2 .* i1 .- 1).*(log.(u(x)) .+ log.(1 .- sort(u(x), rev = true))))
    T = AD
    autores = "Anderson - Darling"
    if prueba == "KS"
        T = KS
        autores = "Kolmogorov - Smirnov"
    end
    if prueba == "CM"
        T = CM
        autores = "Cramér - von Mises"
    end
    tsim = zeros(numsims)
    for j ∈ 1:numsims
        x = simFo(n)
        tsim[j] = T(x)
    end
    tobs = T(xobs)
    pvalue = sum(tsim .> tobs) / numsims
    println("p-valor de la prueba " * autores)
    return pvalue
end

using Statistics
"""
    bootstrap(datos, estadístico; numrep = 10_000, prob = 0.95)

Bootstrap point estimation for the mean of a given statistic `estadístico` (a function to be aplied to a sample given as a vector) and a probability `prob` (0.95 by default) interval estimation, from a given vector `datos` of observed data. `numrep` is the number of boostrap resamples to consider (10,000 by default).

# Example
```
datos = rand(100) # sampling from a Uniform(0,1) distribution
estadístico(x) = (maximum(x) - minimum(x)) / 2
b = bootstrap(datos, estadístico)
b.valor, b.puntual, b.intervalo, b.amplitud
# b.sims contains a vector of bootstrap values of the statistic
```
"""
function bootstrap(datos, estadístico; numrep = 10_000, prob = 0.95)
    n = length(datos)
    sobs = estadístico(datos)
    s = zeros(numrep) # or fill(0.0, numrep)
    for j ∈ 1:numrep
        muestra = rand(datos, n)
        s[j] = estadístico(muestra)
    end
    estimPuntual = mean(s)
    I = (quantile(s, (1 - prob)/2), quantile(s, (1 + prob)/2)) # intervalo de probabilidad prob
    return (valor = sobs, puntual = estimPuntual, intervalo = I, amplitud = extrema(s), sims = s)  
end

"""
    densidad(xobs, h; kernel = "gaussian")

Kernel density estimation for a given random sample `xobs` (as a vector), a smotthing parameter `h` and a chosen kernel function (gaussian by default).

# Example
```
using Distributions, Plots
X = Gamma(2,1)
n = 1000
h = 0.5
xobs = rand(X, n)
f = densidad(xobs, h, kernel = "epanechnikov")
x = range(0, 10, length = 1000)
plot(x, f.(x), lw = 2, label = "kernel density")
plot!(x, pdf(X, x), label = "Gamma(2,1)")
```
"""
function densidad(xobs, h; kernel = "gaussian")
    n = length(xobs)
    w1(x) = exp((-x^2) / 2) / √(2π)
    f1(x) = (1/(n*h)) * sum(w1.((x .- xobs) ./ h))
    w2(x) = (abs(x) < 1) / 2
    f2(x) = (1/(n*h)) * sum(w2.((x .- xobs) ./ h))
    w3(x) = (abs(x) < 1) * (1 - abs(x))
    f3(x) = (1/(n*h)) * sum(w3.((x .- xobs) ./ h))
    w4(x) = (abs(x) < √5) * (3/4) * (1 - x^2/5) / √5
    f4(x) = (1/(n*h)) * sum(w4.((x .- xobs) ./ h))
    w5(x) = (abs(x) < 1) * (15/16) * (1 - x^2)^2
    f5(x) = (1/(n*h)) * sum(w5.((x .- xobs) ./ h))
    if kernel == "gaussian"
        return f1
    end
    if kernel == "rectangular"
        return f2
    end
    if kernel == "triangular"
        return f3
    end
    if kernel == "epanechnikov"
        return f4
    end
    if kernel == "biweight"
        return f5
    end
end

println("Loaded functions: Fn  Tn  EDA  Bn  cuantil_puntual  GoF  bootstrap  densidad")
println("Package dependencies: Distributions (just for Bn)" )
