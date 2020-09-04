length_str = length(Kernel_info);
ii = Kernel_info.true_kernel_log;
ii_size = size(ii,2);
out = zeros(length_str,ii_size);

for i = 1:length_str
    out(i,:) = Kernel_info(i).true_kernel_log;
end
sum_array = sum(out,2);
idx_array = find(sum_array);

centers_selected = centres(idx_array,:);
figure
for ii = 1:size(idx_array,1)
    b = num2str(idx_array(ii));
    c = cellstr(b);
    dx = 0.2; dy = 0.2;
    scatter(centers_selected(ii,1),centers_selected(ii,2));
    text(centers_selected(ii,1)+dx, centers_selected(ii,2)+dy, c);
    hold on
end
hold off
figure
for ii = 1:size(idx_array,1)
    pixel = nonzeros(Kernel_info(idx_array(ii)).true_kernel_idx);
    pixel = unique(pixel);
    pixel_s = length(pixel);
    x_value(1:pixel_s) = idx_array(ii); 
    test(ii,1:pixel_s) = pixel;
    scatter(x_value,pixel)
    x_value = [];
    hold on
    
end



