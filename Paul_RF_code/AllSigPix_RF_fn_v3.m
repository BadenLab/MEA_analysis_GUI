function [RF_coords,Num_Pixels] = AllSigPix_RF_fn_v3(Sig_mat,Sig_Thresh,RF_Ident_Method,Initial_Num_Pix)
%%% Find the all significant pixels RF coordinates & stats
%
% Input:
% Sig_mat (significance matrix) which could be: i) STA_SD, ii) LCA, iii) MIA, or iv) SC.
% Sig_Thresh (significance threshold) which could be: i) STA_SD_Sig_Thresh, ii) LC_Sig_Thresh, iii) MI_Sig_Thresh or iv) SC_Sig_Thresh.
% RF_Ident_Method: i) STA_SD, ii) LCA, iii) MIA, or iv) SC.
% Initial_Num_Pix: number of pixels satisfying the significance threshold
%
% Output:
% RF_coords
% Num_Pixels

% Step 1: Find all significant pixels (ASP)
% Step 2: Find most significant pixel (MSP)
if RF_Ident_Method == 1 || RF_Ident_Method == 2 || RF_Ident_Method == 3
    [RF_row,RF_col] = find(Sig_mat>Sig_Thresh);
    RF_val = NaN(Initial_Num_Pix,1);
    for i = 1:Initial_Num_Pix
        RF_val(i) = Sig_mat(RF_row(i),RF_col(i));
    end
    MSP_index = find(RF_val==max(RF_val));
    RF_row_MSP = RF_row(MSP_index);
    RF_col_MSP = RF_col(MSP_index);
elseif RF_Ident_Method == 4
    [RF_row,RF_col] = find(Sig_mat<Sig_Thresh);
    RF_val = NaN(Initial_Num_Pix,1);
    for i = 1:Initial_Num_Pix
        RF_val(i) = Sig_mat(RF_row(i),RF_col(i));
    end
    MSP_index = find(RF_val==min(RF_val));
    RF_row_MSP = RF_row(MSP_index);
    RF_col_MSP = RF_col(MSP_index);
end

% Take first pixel in list if multiple are most significant
if length(RF_row_MSP)>1
    RF_row_MSP = RF_row_MSP(1);
    RF_col_MSP = RF_col_MSP(1);
end

% Step 3: Loop over known significant pixels to find all remaining significant pixels
%Identify which pixels are connected to the MSP

% Preallocate RF ASP vectors
RF_row_ASP = NaN(Initial_Num_Pix,1);
RF_col_ASP = NaN(Initial_Num_Pix,1);

% Insert MSP into ASP vectors
RF_row_ASP(1) = RF_row_MSP;
RF_col_ASP(1) = RF_col_MSP;

% Define RF_row_remaining and RF_col_remaining and remove confirmed significant pixel
RF_row_remaining = RF_row;
RF_col_remaining = RF_col;

Pix_index_temp   = (RF_row_remaining==RF_row_MSP)&(RF_col_remaining==RF_col_MSP);

RF_row_remaining(Pix_index_temp) = [];
RF_col_remaining(Pix_index_temp) = [];

% Loop variable to keep track of the current index in which to enter new sig pix into RF_row_ASP and RF_col_ASP
sigpix_loop = 2;

% Preallocate vectors of pixels whose neighbours should be checked (or do as matrix)
Pixel_row_to_check = NaN(Initial_Num_Pix,1);
Pixel_col_to_check = NaN(Initial_Num_Pix,1);

% Insert MSP into pixel row to check vectors
Pixel_row_to_check(1) = RF_row_MSP;
Pixel_col_to_check(1) = RF_col_MSP;
Num_Pixels_to_check   = 1; % sum(~isnan(Pixel_row_to_check));

while Num_Pixels_to_check > 0 && sigpix_loop <= Initial_Num_Pix
    
    % Find neighbouring pixels
    [RF_row_neighbour,RF_col_neighbour,num_neigh_loop] ...
        = ASP_neighbour_fn(RF_row_remaining,RF_col_remaining,Pixel_row_to_check(1),Pixel_col_to_check(1));
    
    % Add new sig. pix. to RF_row_ASP and RF_col_ASP
    if num_neigh_loop > 0
        
        RF_row_ASP(sigpix_loop:sigpix_loop+num_neigh_loop-1) = RF_row_neighbour;
        RF_col_ASP(sigpix_loop:sigpix_loop+num_neigh_loop-1) = RF_col_neighbour;
        
        for i=1:num_neigh_loop
            
            Pix_index_temp = (RF_row_remaining==RF_row_neighbour(i))&(RF_col_remaining==RF_col_neighbour(i));
            
            RF_row_remaining(Pix_index_temp) = [];
            RF_col_remaining(Pix_index_temp) = [];
            
        end
        
        Pixel_row_to_check(Num_Pixels_to_check+1:Num_Pixels_to_check+1+num_neigh_loop-1) = RF_row_neighbour;
        Pixel_col_to_check(Num_Pixels_to_check+1:Num_Pixels_to_check+1+num_neigh_loop-1) = RF_col_neighbour;
        
        Num_Pixels_to_check = Num_Pixels_to_check + num_neigh_loop;
        
        sigpix_loop = sigpix_loop + num_neigh_loop;
        
    end
    
    Pixel_row_to_check(1) = [];
    Pixel_col_to_check(1) = [];
    Num_Pixels_to_check   = Num_Pixels_to_check - 1; %sum(~isnan(Pixel_row_to_check));
    
end

% Remove NaN elements
RF_row_ASP(isnan(RF_row_ASP)) = [];
RF_col_ASP(isnan(RF_col_ASP)) = [];

% Output Variables
RF_coords  = [RF_row_ASP,RF_col_ASP];
Num_Pixels = size(RF_coords,1);

