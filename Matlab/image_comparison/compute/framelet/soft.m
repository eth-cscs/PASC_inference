function y = soft(x, T)

y = max(0, x - T) - max(0, - x - T);



