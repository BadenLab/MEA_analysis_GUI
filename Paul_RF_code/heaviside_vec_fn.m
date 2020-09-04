function y = heaviside_vec_fn(x)

length_x = length(x);

y = zeros(length(x),1);

for i = 1:length_x
    
    if x(i)>=0
        y(i)=1;
    else
        y(i)=0;
    end
    
%y(i) = heaviside(x(i)); Could use this, but gives 1/2 at 0 so not quite
%what I want.
    
end