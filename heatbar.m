function out = heatbar(X,Z,varargin)
%Plots data that is too big for classical bar plots in a surface plot with
%or without gaps between datapoints

%% Input control

default_colourmap = gray;
default_gap = true;
p = inputParser;
validMatrixShape = @(x) isnumeric(x) && (size(x,2) == 3);
addRequired(p,'X',@isvector);
addRequired(p,'Z',@ismatrix);
addParameter(p,'Y',@ismatrix);
addParameter(p,'gap',default_gap,@islogical);
addParameter(p,'cmap',default_colourmap,validMatrixShape);
addParameter(p,'figure_handle',@isfigure);


parse(p,X,Z,varargin{:});

%Just for debugging
% out = p;
%% Set variables

%This will remove all variables from the input structure which are not
%beeing used
fields = fieldnames(p.Results);
variables = p.Results;



for ii = 1:length(fields)
   if isa(variables.(fields{ii}),'function_handle')
       variables = rmfield(variables,fields{ii});
   end
       
end

%% Main
X = variables.X;
Z = variables.Z;
if variables.gap
    gap = 2;
else
    gap = 1;
end


%create mesh
dim1 = size(Z,1);
dim2 = size(Z,2);

if dim2 < dim1
    warning("Check Z array dimensions")
end

%Prepare X value for the mesh
if variables.gap
     X = repmat(X,[gap*dim1 ,1]);
     %Extend Z array for higher resolution of the mesh
     Z = repelem(Z,gap,1);
else
     X = repmat(X,[dim1 ,1]);
end


%Prepare Y value for mesh
if isfield(variables, 'Y')
    Y = variables.Y;
    %Check if array has the right orientation
    sz_Y = size(Y);
    if sz_Y(2) == dim1
        Y = Y';
    end
    Y = repmat(Y,[1,dim2]);
    
    if variables.gap
        Y = repelem(Y,gap,1);
    end
        
    
    
else
    [~,Y] = meshgrid(1:dim2,1/gap:1/gap:dim1);
end
%Prepare C data (which defines the colour maps)
C = Z;
if variables.gap
for ii = 1:dim1*2
    if mod(ii,2) == 0
        C(ii,:) = NaN;
    end
    
end
end

if isfield(variables, 'figure_handle')
    heatbar_fig = figure(variables.figure_handle);
    hold on
    heatbar_ax = surf(X,Y,Z,C,'EdgeColor','none','FaceColor','flat');
    colorbar
    colormap(flipud(variables.cmap));
    
    
else
    heatbar_fig = figure;
    hold on
    heatbar_ax = surf(X,Y,Z,C,'EdgeColor','none','FaceColor','flat');
    colorbar
    colormap(flipud(variables.cmap));
end


out.fig = heatbar_fig;
out.ax = heatbar_ax;
end
