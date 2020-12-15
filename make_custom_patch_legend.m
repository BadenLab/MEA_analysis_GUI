function [] = make_custom_patch_legend(colors, labels, varargin)
% [] = MAKE_CUSTOM_PATCH_LEGEND(colors, labels)
% MAKE_CUSTOM_PATCH_LEGEND will plot a custom legend on the current figure.
% It uses a hack involving invisible patch objects that I don't understand.
% COLORS: either a cell array of n valid MATLAB color indicators (color string ('r'), colormap num, or 1x3 RGB vector); 
% or a column vector of n scalars, which will index into the current colormap; 
% or a n x 3 matrix of RGB values for each patch.
% LABELS: a cell array of n string labels for your legend.
% VARARGIN: name/value arguments to pass to legend().

%{
Here are the default MATLAB colors, for reference.
matlab_colors =  [0    0.4470    0.7410; % the matlab colors!
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0.3          0.3        0.3;
    0          0        0    ];
%}

% EXAMPLE: 
%{
xvals = 0:20;
slopes = 1:3;
num_reps = 5;
slopes_to_plot = repelem(slopes, num_reps);
yvals = xvals .* slopes_to_plot' + rand(numel(slopes_to_plot), numel(xvals))*3;
reds =  [0.9961    0.8980    0.8510;
        0.9882    0.6824    0.5686;
        0.9843    0.4157    0.2902;
        0.8706    0.1765    0.1490;
        0.6471    0.0588    0.0824];
blues =  [0.9373    0.9529    1.0000;
        0.7412    0.8431    0.9059;
        0.4196    0.6824    0.8392;
        0.1922    0.5098    0.7412;
        0.0314    0.3176    0.6118];
greens = [0.9294    0.9725    0.9137;
        0.7294    0.8941    0.7020;
        0.4549    0.7686    0.4627;
        0.1922    0.6392    0.3294;
             0    0.4275    0.1725];
color_vals = [reds; blues; greens];
figure
hold on
for i = 1:length(slopes_to_plot)
  plot(yvals(i,:), 'Color', color_vals(i,:))
end
make_custom_patch_legend({'r', 'b', 'g'}, string(slopes))
%}

% Warn the user if it looks like they might be making a mistake.
if isnumeric(colors) && all(size(colors) == [1 3])
   warning('Size of inputted colors argument is 1x3. If you intend for three scalars, use a column vector (3x1).')
end


p = [];

if iscell(colors)
    for ii = 1:length(colors)
        p(ii) = patch(NaN, NaN, colors{ii});
    end
elseif size(colors,2) == 1
    for ii = 1:length(colors)
        p(ii) = patch(NaN, NaN, colors(ii));
    end
elseif size(colors,2) == 3
    for ii = 1:size(colors,1)
        p(ii) = patch(NaN, NaN, colors(ii,:));
    end
else
    error('Bad type for argument colors')
end

legend(p, labels, varargin{:});


