function y = sigmoid(x)
  %vectorized version
  y = 1./(1+exp(-x));
end