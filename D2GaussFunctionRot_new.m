function F = D2GaussFunctionRot_new(x,xdata)
%% x = [Amp, x0, wx, y0, wy, fi]
%% x = [Amp, fi, wx, wy, x0, y0]
%[X,Y] = meshgrid(x,y) 
%  xdata(:,:,1) = X
%  xdata(:,:,2) = Y           
% Mrot = [cos(fi) -sin(fi); sin(fi) cos(fi)]
%%
xdatarot(:,:,1)= xdata(:,:,1)*cos(x(2)) - xdata(:,:,2)*sin(x(2));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(2)) + xdata(:,:,2)*cos(x(2));
x0rot = x(5)*cos(x(2)) - x(6)*sin(x(2));
y0rot = x(5)*sin(x(2)) + x(6)*cos(x(2));

F = x(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(4)^2) )    );

% figure(3)
% alpha(0)
% imagesc(F)
% colormap('gray')
% figure(gcf)%bring current figure to front
% drawnow
% beep
% pause %Wait for keystroke
