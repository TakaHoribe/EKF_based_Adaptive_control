function y = sat(x, max, min)

y = x;
y(y > max) = max;
y(y < min) = min;

end