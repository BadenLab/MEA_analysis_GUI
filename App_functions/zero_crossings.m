function out = zero_crossings(in)


size_in = numel(in);
c = 1;
d = 1;
for ii = 1:size_in-1
    a = in(ii);
    b = in(ii+c);
    
    if a == 0
        continue
    end
    
    while b == 0
        try
        c = c+1;
        b = in(ii+c);
        catch
           break
        end
        
        
            
        
    end
    c = 1;
    
    if abs(diff([a,b])) == 2
        zero_crossing(d) = ii;
        d = d+1;
    end
end

try
out = zero_crossing;
catch
    out = 0;
end
    





end