function [Neigh_coords,Num_Neigh_pixels] = Neighbour_fn(centre_RF_row,centre_RF_col,p)
%%% Find the box RF coordinates & stats
%
% Input:
% centre_RF_row
% centre_RF_col
% p
%
% Output:
% Neigh_coords
% Num_Neigh_pixels


Num_Layers = 1; % Use 1 layer as looking for direct neighbour.


% Choose box around this pixel
if centre_RF_col - Num_Layers >= 2
    min_RF_col = centre_RF_col - Num_Layers;
else
    min_RF_col = 1;
end

if centre_RF_col + Num_Layers <= p.stim_columns - 1
    max_RF_col = centre_RF_col + Num_Layers;
else
    max_RF_col = p.stim_columns;
end

if centre_RF_row - Num_Layers >= 2
    min_RF_row = centre_RF_row - Num_Layers;
else
    min_RF_row = 1;
end

if centre_RF_row + Num_Layers <= p.stim_rows - 1
    max_RF_row = centre_RF_row + Num_Layers;
else
    max_RF_row = p.stim_rows;
end

Num_RF_rows      = max_RF_row - min_RF_row + 1;
Num_RF_cols      = max_RF_col - min_RF_col + 1;
Num_Neigh_pixels = Num_RF_rows*Num_RF_cols - 1;

RF_coords_centre = [centre_RF_row,centre_RF_col];

% Find row and column coords of RF
RF_coords = NaN((max_RF_row-min_RF_row+1)*(max_RF_col-min_RF_col+1),2);

for i = min_RF_row:max_RF_row
    for j = min_RF_col:max_RF_col
       RF_coords((i-min_RF_row)*(max_RF_col-min_RF_col+1) + j-min_RF_col + 1,:) = [i,j]; 
    end
end

% Find non-centre coordinates
[~,RF_centre_index,~]           = intersect(RF_coords,RF_coords_centre,'rows');
Neigh_coords                    = RF_coords;
Neigh_coords(RF_centre_index,:) = [];

