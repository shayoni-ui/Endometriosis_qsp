function y = sat(x)
  %vectorized version
  y = max(1, min(0, x));
end