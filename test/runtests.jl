using Midas
using TimeSeries
using Base.Test

# write your own tests here
months = [Date(2000, m, 1) for m in 1:12]
quarters = [Date(2000, 3, 31), Date(2000, 6, 30),
Date(2000, 9, 30), Date(2000, 12, 31)]
lf_data = TimeArray(quarters, [1.0, 2.0, 3.0, 4.0])
hf_data = TimeArray(months, collect(1.0:12.0))

@test mixfrequencies(lf_data, hf_data, 3, 1).values[1, :] == [4.0, 3.0, 2.0]
@test datafreq(lf_data) == fquarter
@test datafreq(hf_data) == fmonth


hf_data = readtimearray("./data/farmpay.csv", format="yyyy-mm-dd")
lf_data = readtimearray("./data/gdp.csv", format="yyyy-mm-dd")

@test mixfrequencies(lf_data, hf_data, 3, 0).values[1, :] == [43545.0,43379.0,43396.0]
@test mixfrequencies(lf_data, hf_data, 3, 0).values[end, :] == [130974.0,130757.0,130563.0]

@test isapprox(collect(beta_weights_es(3, 1, 5)), [0.941176,0.0588238,9.4118e-25], rtol=1e-6)

x = mixfrequencies(lf_data, hf_data, 3, 1)
w = collect(beta_weights_es(3, 1, 5))
z = x.values * w
