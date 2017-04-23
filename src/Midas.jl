module Midas

using TimeSeries
using LsqFit
using Optim

@enum FREQ fday=1 fmonth=2 fquarter=3 funknown=4

export mixfrequencies,
       datafreq, FREQ, fquarter, fmonth,
       beta_weights_es, xweighted, ssr_func, example1, example2,
       jacobian_wx

include("mix.jl")
include("fit.jl")

function xweighted(x, theta1, theta2)
  w = beta_weights_es(size(x,2), theta1, theta2)
  return x * w, repmat(w', size(x,2), 1)
end

function beta_weights_es(n, theta1, theta2)
  u = linspace(1e-6, 1.0 - 1e-6, n)

  beta_vals = u.^(theta1 - 1).*(1 - u).^(theta2 - 1)

  return beta_vals / sum(beta_vals)
end

function example1()
  hf_data = readtimearray("./data/farmpay.csv", format="yyyy-mm-dd")
  lf_data = readtimearray("./data/gdp.csv", format="yyyy-mm-dd")

  hf_g = percentchange(hf_data, :log) * 100
  lf_g = percentchange(lf_data, :log) * 100

  y, x, yfc, xfc = mixfrequencies(lf_g,
                                  hf_g,
                                  3, 0, 1,
                                  start_date=Date(1985,1,1),
                                  end_date=Date(2009,1,1))

  xw, w = xweighted(x.values, 1, 5)

  yvals = y.values[:, 1]

  a, b = linreg(xw, yvals)

  f = ssr_func(x.values, yvals)
  #nlsolve(f, [0.6, 1.9, 1., 5.])
  g! = ssr_grad_func(x.values, yvals)

  res = optimize(f, g!, [a, b, 1., 5.], LBFGS())

  return res
end

end # module
