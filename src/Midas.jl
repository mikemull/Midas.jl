module Midas

using TimeSeries
using NLsolve
using LsqFit
using Optim

@enum FREQ fday=1 fmonth=2 fquarter=3 funknown=4

export mixfrequencies,
       datafreq, FREQ, fquarter, fmonth,
       beta_weights_es, xweighted, ssr_func, example1, example2,
       jacobian_wx

include("mix.jl")

function xweighted(x, theta1, theta2)

  w = beta_weights_es(size(x,2), theta1, theta2)
  return x * w, repmat(w', size(x,2), 1)
end

function beta_weights_es(n, theta1, theta2)
  u = linspace(1e-6, 1.0 - 1e-6, n)

  beta_vals = u.^(theta1 - 1).*(1 - u).^(theta2 - 1)

  return beta_vals / sum(beta_vals)
end

function jacobian_wx(x, theta1, theta2)
  eps = 1e-6

  xt1p, w = xweighted(x, theta1 + eps / 2, theta2)
  xt1m, w = xweighted(x, theta1 - eps / 2, theta2)
  jt1 = (xt1p - xt1m) / eps

  xt2p, w = xweighted(x, theta1, theta2 + eps / 2)
  xt2m, w = xweighted(x, theta1, theta2 - eps / 2)
  jt2 = (xt2p - xt2m) / eps

  return [jt1 jt2]
end

function ssr_func(x, y)
  function ssr(a)
    xw, w = xweighted(x, a[3], a[4])
    error = y - a[1] - a[2] * xw
    return (error' * error)[1]
  end
end

function ssr_grad_func(x, y)
  function ssr_grad(a, storage)
    jwx = jacobian_wx(x, a[3], a[4])
    xw, w = xweighted(x, a[3], a[4])

    error = y - a[1] - a[2] * xw

    jac_e = [ones(size(xw)) xw (a[2] * jwx)]

    jac = zeros(size(jac_e))
    for i=1:length(a)
        jac[:, i] = -2 * jac_e[:,i] .* error
    end
    grad=sum(jac, 1)'
    for i in eachindex(grad)
      storage[i] = grad[i]
    end
  end
end

function example1()
  hf_data = readtimearray("./data/farmpay.csv", format="yyyy-mm-dd")
  lf_data = readtimearray("./data/gdp.csv", format="yyyy-mm-dd")

  hf_g = percentchange(hf_data, :log) * 100
  lf_g = percentchange(lf_data, :log) * 100

  x = mixfrequencies(lf_g[Date(1985,1,1):Date(2009,1,1)], hf_g, 3, 1)

  xw, w = xweighted(x.values, 1, 5)
  y = lf_g[Date(1985,1,1):Date(2009,1,1)].values[:,1]
  a, b = linreg(xw, y)

  f = ssr_func(x.values, y)
  #nlsolve(f, [0.6, 1.9, 1., 5.])
  g! = ssr_grad_func(x.values, y)

  res = optimize(f, g!, [a, b, 1., 5.], LBFGS())

  return res
end

end # module
