function true_unique = unique_test (input)
input_size = size(input,2);
true_unique = zeros(1,input_size);
%test_size = 0.05*input_size;

if input_size<1
    true_unique = ones(1,input_size);
else


for ii = 1:input_size
    input_test = ceil05(input / input(ii),0.1);
    input_test = nnz(input_test == 1);
    input_test
   
    
    
    if input_test == 1
        true_unique(ii) = true;
    elseif input_test > 1
        true_unique(ii) = false;
    end
    
    
end
end
out = num2str(nnz(true_unique)+1);
disp([out,' stimuli have been found']);

        



end