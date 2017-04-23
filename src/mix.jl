
function mixfrequencies(lf_data,
                        hf_data,
                        xlag,
                        ylag,
                        horizon;
                        start_date=nothing,
                        end_date=nothing)

  if start_date == nothing
    start_date = lf_data.timestamp[ylag + 1]
  end

  if end_date == nothing
    end_date = lf_data.timestamp[end]
  end

  min_y_date = lf_data.timestamp[1 + ylag]
  min_x_date = hf_data.timestamp[1 + xlag + horizon]

  if ylag > 0
    ylags = lag(lf_data, 1)
    for lag in 2:ylag
      yl = lag(lf_data, lag)
      rename(yl, "ylag$(lag)")
      merge(ylags, yl)
    end
  end

  y = lf_data[start_date:end_date]
  y_forecast = lf_data[end_date:lf_data.timestamp[end]]
  if ylag > 0
    ylags_forecast = ylags[end_date:lf_data.timestamp[end]]
    ylags = ylags[start_date:end_date].values
  else
    ylags_forecast = Array{Float64}(length(y_forecast), 0)
    ylags = Array{Float64}(length(y), 0)
  end

  hfv = [ylags xlags(y, hf_data, xlag, horizon)]
  hfv_forecast = [ylags_forecast xlags(y_forecast, hf_data, xlag, horizon)]

  return  y, TimeArray(y.timestamp, hfv), y_forecast, TimeArray(y_forecast.timestamp, hfv_forecast)
end


function xlags(y, hf_data, xlag, horizon)
  hf_lags = zeros(length(y), xlag)
  irow = 1
  for (t,v) in y
    start_hf = find(x -> x >= t, hf_data.timestamp)[1]
    hf_lags[irow, :] = hf_data.values[start_hf - horizon:-1:start_hf - horizon - xlag + 1]
    irow += 1
  end
  return hf_lags
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
