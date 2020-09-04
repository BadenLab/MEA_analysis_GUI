function [RF_coords_output,Num_RF_pixels,Gaussian_mean,Gaussian_covar,eval_descend_vec,evec_descend_mat,...
    Angle_Major_Axis,Angle_Minor_Axis,Axis_eval_Ratio,Rel_Ellipse_Area,Abs_Ellipse_Area,gmPDF,Thresh_height]...
    = Gaussian_RF_fn_v3(RF_coords_input,Num_Sig_Pix,p)
%%% Find the Gaussian RF coordinates & stats
%
% Input:
% RF_coords_input: ASP coords
% Num_Sig_Pix: Number of significant pixels
% p.RF_SDs: number of SDs at which to place limits of Gaussian
%
% Output:
% RF_coords_output
% Num_RF_pixels
% Gaussian_mean
% Gaussian_covar
% Angle_Major_Axis
% Angle_Minor_Axis
% Axis_eval_Ratio
% Rel_Ellipse_Area
% Abs_Ellipse_Area
% gmPDF
% Thresh_height

if Num_Sig_Pix>2
    
    Sig_RF_row = RF_coords_input(:,1);
    Sig_RF_col = RF_coords_input(:,2);
    
else % Num_Sig_Pix<=2
    
    Sig_RF_row_temp = RF_coords_input(:,1);
    Sig_RF_col_temp = RF_coords_input(:,2);
    
    Sig_RF_row = NaN(4*Num_Sig_Pix,1);
    Sig_RF_col = NaN(4*Num_Sig_Pix,1);
    
    for i = 1:Num_Sig_Pix
        
        Sig_RF_row(4*(i-1)+1:4*i) =  [Sig_RF_row_temp(i) - 0.25; Sig_RF_row_temp(i) - 0.25; Sig_RF_row_temp(i) + 0.25; Sig_RF_row_temp(i) + 0.25];
        Sig_RF_col(4*(i-1)+1:4*i) =  [Sig_RF_col_temp(i) - 0.25; Sig_RF_col_temp(i) + 0.25; Sig_RF_col_temp(i) - 0.25; Sig_RF_col_temp(i) + 0.25];
        
    end
    
end

% Fit a Gaussian to the significant pixels
GMModel           = fitgmdist([Sig_RF_col,Sig_RF_row],1,'CovarianceType','full',...
                    'RegularizationValue',1e-5,'Replicates',20,'Options',statset('MaxIter',1e4));
Gaussian_mean  = GMModel.mu;
Gaussian_covar = GMModel.Sigma;

[evec_mat,eval_mat] = eig(Gaussian_covar); % The eigenvectors are normalized so that the 2-norm of each is 1.
eval_vec            = diag(eval_mat);
[eval_descend_vec,eval_sort_indices_vec] = sort(eval_vec,'descend');
evec_descend_mat                         = evec_mat(:,eval_sort_indices_vec);
Angle_Major_Axis = atan(evec_descend_mat(2,1)/evec_descend_mat(1,1)); % Inverse tangent in radians
Angle_Minor_Axis = atan(evec_descend_mat(2,2)/evec_descend_mat(1,2)); % For real values of X, atan(X) returns values in the interval [-π/2, π/2].
Axis_eval_Ratio  = eval_descend_vec(1)/eval_descend_vec(2);    % Divide bigger by smaller, so values in range [1,inf)
Rel_Ellipse_Area = pi*eval_descend_vec(1)*eval_descend_vec(2); % Not absolute area for chosen CI or SD, but measure which allows comparison between Gauss RFs
Abs_Ellipse_Area = Rel_Ellipse_Area*p.RF_SDs^2; % Area up to chosen SD contour.

% Find the height of the Gaussian at chosen num of SDs out from mean
sigma_X       = sqrt(GMModel.Sigma(1,1));
sigma_Y       = sqrt(GMModel.Sigma(2,2));
rho           = GMModel.Sigma(1,2)/(sigma_X*sigma_Y);
Thresh_height = (1/(2*pi*sigma_X*sigma_Y*sqrt(1 - rho^2)))*exp(-(p.RF_SDs^2)/(2*(1 - rho^2)));

% Find which pixels are within the threshold SD countour
[x_mat,y_mat]    = meshgrid(1:1:p.stim_columns,1:1:p.stim_rows);
mu_X             = GMModel.mu(1);
mu_Y             = GMModel.mu(2);

Pixel_height_mat = ((x_mat - mu_X).^2)/sigma_X^2 + ((y_mat - mu_Y).^2)/sigma_Y^2 ...
                   - (2*rho*(x_mat - mu_X).*(y_mat - mu_Y))/(sigma_X*sigma_Y);
               
[RF_row,RF_col]  = find(Pixel_height_mat<=p.RF_SDs^2);           
RF_coords_output = [RF_row,RF_col];
Num_RF_pixels = length(RF_row);

gmPDF = @(x,y)pdf(GMModel,[x y]);

