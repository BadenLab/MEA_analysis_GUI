function [RF_row_neighbour,RF_col_neighbour,Num_neigh] ...
    = ASP_neighbour_fn(RF_row,RF_col,RF_row_centre,RF_col_centre)

% Input:
% RF_row: row indices of all pixels passing significance threshold
% RF_col: col indices of all pixels passing significance threshold
% RF_row_centre: significant pixel whose neighbours we're investigating
% RF_col_centre: significant pixel whose neighbours we're investigating
%
% Output:
% RF_row_neighbour: row indices of significant neighbouring pixels
% RF_col_neighbour: col indices of significant neighbouring pixels
% Num_neigh: number of significant neighbouring pixels

% Logical vectors to find the rows and columns of neighbouring sig. pix.
row_logical_inner = abs(RF_row-RF_row_centre)>0;
col_logical_inner = abs(RF_col-RF_col_centre)>0;
row_logical_outer = abs(RF_row-RF_row_centre)<=1;
col_logical_outer = abs(RF_col-RF_col_centre)<=1;

% Indices of neighbouring sig. pix.
Neighbour_logical = logical((row_logical_inner|col_logical_inner).*(row_logical_outer&col_logical_outer));

% Number of neighbouring sig. pix.
Num_neigh = sum(Neighbour_logical);

% Add new sig. pix. to RF_row_ASP and RF_col_ASP
if Num_neigh > 0
    RF_row_neighbour = RF_row(Neighbour_logical);
    RF_col_neighbour = RF_col(Neighbour_logical);
else 
    RF_row_neighbour = [];
    RF_col_neighbour = [];
end

