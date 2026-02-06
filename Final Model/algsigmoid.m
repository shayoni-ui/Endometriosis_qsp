function y = algsigmoid(x)
  %vectorized version
  k=1;
  y = x./(1+ x.^k).^1/k;
end