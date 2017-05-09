module Midas

using TimeSeries
using LsqFit
using Optim
using MultivariateStats

@enum FREQ fday=1 fmonth=2 fquarter=3 funknown=4

export mixfrequencies, example1, example1data

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

function estimate(y, yl, x)
  xw, w = xweighted(x.values, 1, 5)

  yvals = y.values[:, 1]

  c = llsq([xw yl.values], y.values)

  f = ssr_func(x.values, yvals, yl.values)
  g! = ssr_grad_func(x.values, yvals, yl.values)

  optimize(f, g!, [c[1:2]; [1; 5]; c[3:end]], LBFGS())
end

function forecast(xfc, yfcl, res)
  copt = Optim.minimizer(res)

  xw, w = xweighted(xfc.values, copt[3], copt[4])
  yf = copt[1] +  copt[2] * xw + copt[end] * yfcl.values

  return yf
end

function example1data()
  hf_data = readtimearray("./data/farmpay.csv", format="yyyy-mm-dd")
  lf_data = readtimearray("./data/gdp.csv", format="yyyy-mm-dd")

  hf_g = percentchange(hf_data, :log) * 100
  lf_g = percentchange(lf_data, :log) * 100

  y, yl, x, yfc, yfcl, xfc = mixfrequencies(lf_g,
                                  hf_g,
                                  3, 1, 1,
                                  start_date=Date(1985,1,1),
                                  end_date=Date(2009,1,1))
end

function example1()
  y, yl, x, yfc, yfcl, xfc = example1data()

  res = estimate(y, yl, x)

  yf = forecast(xfc, yfcl, res)

  return yf
end

end # module
