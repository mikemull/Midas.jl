module Midas

# package code goes here
using TimeSeries

@enum FREQ fday=1 fmonth=2 fquarter=3 funknown=4

export mixfrequencies,
       datafreq, FREQ, fquarter, fmonth,
       beta_weights_es, xweighted

function mixfrequencies(lf_data, hf_data, xlag, horizon)

  start_hf = find(x -> x >= lf_data.timestamp[1], hf_data.timestamp)[1]
  start_hf = start_hf - xlag - horizon + 1

  hf_gt = find(x -> x >= lf_data.timestamp[end], hf_data.timestamp)
  if length(hf_gt) == 0
    end_hf = length(hf_data)
  else
    end_hf = hf_gt[1]
  end
  end_hf -= horizon
  # end hf
  # num_obs
  num_obs = length(lf_data)

  hf_lags = zeros(num_obs * xlag)
  for i in eachindex(hf_data.values[start_hf:end_hf])
           hf_lags[i] = hf_data.values[i+start_hf-1]
  end

  hfv = transpose(reshape(hf_lags, (xlag, num_obs)))
  hfv = hfv[:,end:-1:1]

  return  TimeArray(lf_data.timestamp, hfv)
end

function datafreq(tsdata)

  deltas = lag(tsdata).timestamp[1:end] - tsdata.timestamp[1:end-1]

  mind = Int(minimum(deltas))

  if 90 <= mind <= 92
    return fquarter
  elseif 28 <= mind <= 31
    return fmonth
  elseif mind == 1
    return fday
  end
  return funknown
end


function xweighted(x)

  w = beta_weights_es(size(x,2), 1, 5)

  return x * w, repmat(w', size(x,2), 1)
end


function beta_weights_es(n, theta1, theta2)
  u = linspace(1e-6, 1.0 - 1e-6, n)

  beta_vals = u.^(theta1 - 1).*(1 - u).^(theta2 - 1)

  return beta_vals / sum(beta_vals)
end

end # module
