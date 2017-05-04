
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

function ssr_func(x, y, yl)
  function ssr(a)
    xw, w = xweighted(x, a[3], a[4])
    if yl == nothing
      error = y - a[1] - a[2] * xw
    else
      error= y - a[1] - a[2] * xw - yl * a[5]
    end
    return (error' * error)[1]
  end
end

function ssr_grad_func(x, y, yl)
  function ssr_grad(a, storage)
    jwx = jacobian_wx(x, a[3], a[4])
    xw, w = xweighted(x, a[3], a[4])

    if yl == nothing
      error = y - a[1] - a[2] * xw
      jac_e = [ones(size(xw)) xw (a[2] * jwx)]
    else
      error = y - a[1] - a[2] * xw - yl * a[5]
      jac_e = [ones(size(xw)) xw (a[2] * jwx) yl]
    end

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
