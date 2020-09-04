function [y, a] = second_max( x )
   [y, a] = max(x(x<max(x)));
end