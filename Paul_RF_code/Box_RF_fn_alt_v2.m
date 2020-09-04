function [RF_coords,Num_RF_rows,Num_RF_cols,Num_RF_pixels] = Box_RF_fn_alt_v2(Centre_RF_coords,p)
%%% Find the box RF coordinates & stats
%
% Input:
% Centre_RF_coords: coordinates of most significant pixel for each colour
%
% Output:
% RF_coords
% Num_RF_rows
% Num_RF_cols
% Num_RF_pixels

% Choose box around these pixels
min_RF_row_temp = min(Centre_RF_coords(:,1));
max_RF_row_temp = max(Centre_RF_coords(:,1));
min_RF_col_temp = min(Centre_RF_coords(:,2));
max_RF_col_temp = max(Centre_RF_coords(:,2));

if min_RF_row_temp - p.RF_layers >= 2
    min_RF_row = min_RF_row_temp - p.RF_layers;
else
    min_RF_row = 1;
end

if min_RF_col_temp - p.RF_layers >= 2
    min_RF_col = min_RF_col_temp - p.RF_layers;
else
    min_RF_col = 1;
end

if max_RF_row_temp - p.RF_layers <= p.stim_rows - 1
    max_RF_row = max_RF_row_temp + p.RF_layers;
else
    max_RF_row = p.stim_rows;
end

if max_RF_col_temp - p.RF_layers <= p.stim_columns - 1
    max_RF_col = max_RF_col_temp + p.RF_layers;
else
    max_RF_col = p.stim_columns;
end


Num_RF_rows   = max_RF_row - min_RF_row + 1;
Num_RF_cols   = max_RF_col - min_RF_col + 1;
Num_RF_pixels = Num_RF_rows*Num_RF_cols;

% Find row and column coords of RF
RF_coords = NaN((max_RF_row-min_RF_row+1)*(max_RF_col-min_RF_col+1),2);

for i = min_RF_row:max_RF_row
    for j = min_RF_col:max_RF_col
       RF_coords((i-min_RF_row)*Num_RF_cols + j-min_RF_col + 1,:) = [i,j]; 
    end
end