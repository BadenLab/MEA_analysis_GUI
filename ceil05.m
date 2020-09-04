function [rounded05] = ceil05 (input, factor)


base = round(input/factor);

rounded05 = base*factor;


end