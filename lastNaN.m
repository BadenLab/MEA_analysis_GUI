function [last_idx, out] = lastNaN (in,varargin)

p = inputParser;
default_dim = 1;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%validMatrixShape = @(x) isnumeric(x) && (size(x,1) ==1) && (size(x,2) == 2);
addRequired(p,'in',@ismatrix);
addOptional(p,'dim',default_dim,validScalarPosNum);


parse(p,in,varargin{:}); 

in = p.Results.in;
dim = p.Results.dim;


size_in = size(in);
test_size = size_in(1)*size_in(2);
if test_size == 0
    error('array is empty')
end


if numel(size_in) == 3
    error('Cant compute last NaN for 3D array')
end


loop_out = zeros(1,size_in(2));
for ii = 1:size_in(2)
    if size_in(1) == 1
        test_array = in(1,:);
        nan_test = ~isnan(test_array);
        loop_out(1,ii) = nnz(nan_test);
    elseif size_in(2) == 1
        test_array = in(:,1);
        nan_test = ~isnan(test_array);
        loop_out(1,ii) = nnz(nan_test);
    elseif numel(size_in) == 2 && dim == 2
        test_array = in(ii,:);
        nan_test = ~isnan(test_array);
        loop_out(1,ii) = nnz(nan_test);
    elseif numel(size_in) == 2 && dim == 1
        test_array = in(:,ii);
        nan_test = ~isnan(test_array);
        loop_out(1,ii) = nnz(nan_test);
    end
   
        
        
end

last_idx = max(loop_out);



    if size_in(1) == 1
       out = in(1,1:last_idx);
    elseif size_in(2) == 1
       out = in(1:last_idx,1);
    elseif numel(size_in) == 2 && dim == 2
       out = in(:,1:last_idx);
    elseif numel(size_in) == 2 && dim == 1
       out = in(1:last_idx,:);
    end
   



end