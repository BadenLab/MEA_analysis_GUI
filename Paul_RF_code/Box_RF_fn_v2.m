function [RF_coords,RF_coords_centre,RF_coords_noncentre,Num_RF_rows,Num_RF_cols,Num_RF_pixels] = Box_RF_fn_v2(Sig_mat,p,RF_Ident_Method)
%%% Find the box RF coordinates & stats
%
% Input:
% Sig_mat (significance matrix) which could be: i) STA_SD, ii) LCA, or iii) MIA.
%
% Output:
% RF_coords
% RF_coords_centre
% RF_coords_noncentre
% Num_RF_rows
% Num_RF_cols
% Num_RF_pixels 

% Find (centre) pixel with max std dev
if RF_Ident_Method == 1 || RF_Ident_Method == 2 || RF_Ident_Method == 3
    max_Sig_mat                   = max(Sig_mat,[],'all');
    [centre_RF_row,centre_RF_col] = find(Sig_mat==max_Sig_mat);
    length_centre_RF_row = length(centre_RF_row);
    if length_centre_RF_row>1
        %         mean_vec = NaN(length_centre_RF_row,1);
        %         for i = 1:length_centre_RF_row
        %             mean_vec(i) = mean(Sig_mat(max(centre_RF_row-1,1):min(centre_RF_row+1,p.stim_rows),max(centre_RF_col-1,1):min(centre_RF_col+1,p.stim_columns)),'all');
        %         end
        %         max_index = find(mean_vec == max(mean_vec));
        %         centre_RF_row = centre_RF_row(max_index);
        %         centre_RF_col = centre_RF_col(max_index);
        centre_RF_row = centre_RF_row(1);
        centre_RF_col = centre_RF_col(1);
    end
elseif RF_Ident_Method == 4
    min_Sig_mat                   = min(Sig_mat,[],'all');
    [centre_RF_row,centre_RF_col] = find(Sig_mat==min_Sig_mat);
    length_centre_RF_row = length(centre_RF_row);
    if length_centre_RF_row>1
        centre_RF_row = centre_RF_row(1);
        centre_RF_col = centre_RF_col(1);
    end
end


% Choose box around this pixel
if centre_RF_col - p.RF_layers >= 2
    min_RF_col = centre_RF_col - p.RF_layers;
else
    min_RF_col = 1;
end

if centre_RF_col + p.RF_layers <= p.stim_columns - 1
    max_RF_col = centre_RF_col + p.RF_layers;
else
    max_RF_col = p.stim_columns;
end

if centre_RF_row - p.RF_layers >= 2
    min_RF_row = centre_RF_row - p.RF_layers;
else
    min_RF_row = 1;
end

if centre_RF_row + p.RF_layers <= p.stim_rows - 1
    max_RF_row = centre_RF_row + p.RF_layers;
else
    max_RF_row = p.stim_rows;
end

Num_RF_rows   = max_RF_row - min_RF_row + 1;
Num_RF_cols   = max_RF_col - min_RF_col + 1;
Num_RF_pixels = Num_RF_rows*Num_RF_cols;

RF_coords_centre = [centre_RF_row,centre_RF_col];

% Find row and column coords of RF
RF_coords = NaN((max_RF_row-min_RF_row+1)*(max_RF_col-min_RF_col+1),2);

for i = min_RF_row:max_RF_row
    for j = min_RF_col:max_RF_col
       RF_coords((i-min_RF_row)*(max_RF_col-min_RF_col+1) + j-min_RF_col + 1,:) = [i,j]; 
    end
end

% Find non-centre coordinates
[~,RF_centre_index,~]                     = intersect(RF_coords,RF_coords_centre,'rows');
RF_coords_noncentre                    = RF_coords;
RF_coords_noncentre(RF_centre_index,:) = [];

