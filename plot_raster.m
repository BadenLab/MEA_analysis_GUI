function plot_raster(in,varargin)

%% Manage inputs

default_dim = 1;
p = inputParser;
valid_diminput = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x<3);
matrix_check = @(x) ismatrix(x);
figure_check = @(x) isfigure(x);
addRequired(p,'in', matrix_check);
addOptional(p,'dim', default_dim, valid_diminput)
addOptional(p,'ax', figure_check);

parse(p,in,varargin{:});
fields = fieldnames(p.Results);
variables = p.Results;


for ii = 1:length(fields)
   if isa(variables.(fields{ii}),'function_handle')
       variables = rmfield(variables,fields{ii});
   end
       
end


%% Main
in = variables.in;
dim = variables.dim;

if isfield(variables,'ax')
   
    ax = variables.ax;
    if dim == 2
        in = in';
    end
    
    yvalue = ones(length(in),1)*0.25;


    
    for ii = 1:size(in,2)
        scatter(ax,in(:,ii),yvalue+0.25*(ii-1),'.','k');
        hold(ax,'on')

    end
    ylim(ax,[0,1+ii*0.25])
    xlim(ax,[-0.5,max(in,[],'all')+1])
    
    hold(ax,'off')
    
    
    
    
    
else

   
    if dim == 2
        in = in';
    end
    
    yvalue = ones(length(in),1)*0.25;


    
    for ii = 1:size(in,2)
        scatter(in(:,ii),yvalue+0.25*(ii-1),'.','k');
        hold on

    end
    ylim([0,1+ii*0.25])
    xlim([-0.5,max(in,[],'all')+1])
    
    hold off
    
end

end