function out = calculate_rec_field (savepath, add_info)





stim_idx = add_info.stim_idx;
%For gui later
one_picture = 0;
try
    
     M = matfile(findfile_app(stim_idx,savepath,'cell_kernel_overview'),...
         'Writable',false);
     test = size(M.cell_kernel_overview);
     
     if sum(test) == 0
         out = 0; % Return 0 because field doesnt exist 
         return 
     end
     
catch
    %If the structure doesnt exist, function gets returned here
    disp('Kernel overview info not found, check if kernel_analysis.m has run')
    out = 0;
    return
end



cell_kernel_overview = M.cell_kernel_overview;
clear M

%Delet the cells that didnt pass the peak test
a = 1;
for ii = 1:numel(cell_kernel_overview)
    try
       if isnan(cell_kernel_overview(a).luminance_kernel)
           cell_kernel_overview(a) = [];
       else 
           a = a+1;
       end
    end
end
    
nr_cells = length(cell_kernel_overview);


%Load Kernel location
S = load(findfile_app(stim_idx,savepath,'Kernel_location'));
Kernel_location = S.Kernel_location;
Kernel_idx = str2double(Kernel_location(:,2));


%% Test
gaussian1 = fspecial('Gaussian', 20, 30);
gaussian2 = fspecial('Gaussian', 20, 4);
dog = gaussian1 - gaussian2;

%Preallocate structure
clear Receptive_field_overview
for ss = 1:nr_cells
    Receptive_field_overview(ss).picture = [];
    Receptive_field_overview(ss).fitresult = [];
    Receptive_field_overview(ss).zfit = [];
    Receptive_field_overview(ss).fiterr = [];
    Receptive_field_overview(ss).zerr = [];
    Receptive_field_overview(ss).resnorm = [];
    Receptive_field_overview(ss).rr = [];
end



for ii = 1:nr_cells
    ii
    detailed_info = cell_kernel_overview(ii).detailed_info;
    
    %Check which colours are active
    
    active_colour = find(detailed_info.active_colour);
%     %Check which colour is dominant
%     dominant_colour = find(detailed_info.colour_dominant);
%     Receptive_field_overview(ii).dom_colour = dominant_colour;
%     
%     %Check if the on or off peak is higher
%     if detailed_info.ON_dominance(dominant_colour) == 1
%         Receptive_field_overview(ii).peak_loc = detailed_info.position_on(dominant_colour);
%     else
%         Receptive_field_overview(ii).peak_loc = detailed_info.position_off(dominant_colour);
%     end
    
    %Load Kernel location
     cell_idx = cell_kernel_overview(ii).cell_idx;
    
     
     
    idx = Kernel_idx == cell_idx;
    location = Kernel_location(idx,1);
    S = load(location);
    Kernel = S.Kernels;
   
    Kernel_size = size(Kernel);
    %Normalize all traces 
    channel_id = cell_kernel_overview(ii).channel_idx;
    channel_id_pixel = channel_id_to_pixel(channel_id,Kernel_size(2));

   %% Calculate receptive field based on one picture at maximum kernel
   if one_picture == 1
    
    %Get the picture from the moment set by peaks loc
%     idx = Receptive_field_overview(ii).peak_loc;
%     picture = squeeze(Kernel(idx,:,:));
   
    
    %znormalize picture
    %for this we take a subset of  20% of the closest to mean active
    %channels calculate mean and std for them and znormalize based on that
    picture_norm = NaN(Kernel_size(2),Kernel_size(3));
    
    for kk = active_colour'
        %Check which colour is dominant
        if isnan(detailed_info.ON_dominance(kk))
            idx = detailed_info.position_off(kk);
        elseif isnan(detailed_info.OFF_dominance(kk))
           idx = detailed_info.position_on(kk);
        elseif detailed_info.ON_dominance(kk) > detailed_info.OFF_dominance(kk)
            idx = detailed_info.position_on(kk);
        else
            idx = detailed_info.position_off(kk);
        end
        
       picture = squeeze(Kernel(idx,:,:));
       picture_mean = mean(picture(:,kk));
       picture_std =  std(picture(:,kk));
       
       picture_norm(:,kk) = znormalise(picture(:,kk),picture_mean,picture_std);
           
    end
    
    %reshape picture to image as it was
    reshape_factor = sqrt(size(Kernel,2));
    picture_reshape = reshape(picture_norm,[reshape_factor,reshape_factor,Kernel_size(3)]);
    picture_size = max(size(picture_reshape));
    Receptive_field_overview(ii).picture = picture_reshape;
    
   end
   
   
   %% Calculate picture based on std through time
   if one_picture == 0
       
      %Calculate std through time
      
      picture_norm = NaN(Kernel_size(2),Kernel_size(3));
      
      for kk = active_colour'
          
      picture_norm(:,kk) = squeeze(std(Kernel(:,:,kk),1,'omitnan'));
      
      end
      %reshape picture to image as it was
      reshape_factor = sqrt(size(Kernel,2));
      picture_reshape = reshape(picture_norm,[reshape_factor,reshape_factor,Kernel_size(3)]);
      picture_size = max(size(picture_reshape));
      Receptive_field_overview(ii).picture_std = picture_reshape;
      
       
       
   end
        
    
    %% Interpolate and gaussian model
    MdataSize = picture_size-1;
    [X,Y] = meshgrid(0:MdataSize);
    [X_new,Y_new] = meshgrid(0:0.25:MdataSize);
    
    picture_interp = zeros(length(X_new),length(X_new),Kernel_size(3));
    
    
    fitresult = NaN(Kernel_size(3),7);
    fiterr = NaN(Kernel_size(3),7);
    zfit = NaN(length(X_new),length(X_new),Kernel_size(3));
    zerr = NaN(length(X_new),length(X_new),Kernel_size(3));
    resnorm = NaN(Kernel_size(3),1);
    rr = NaN(Kernel_size(3),1);
    for kk = 1:Kernel_size(3)
        
        %skip NaN picture
        if nnz(isnan(picture_norm(:,kk))) > 10
            continue
        else
        
        picture_interp(:,:,kk) = interp2(X,Y,picture_reshape(:,:,kk),X_new,Y_new);
        [fitresult(kk,:), zfit(:,:,kk), fiterr(kk,:), zerr(:,:,kk), resnorm(kk,1),...
            rr(kk,1)] = fmgaussfit(X_new,Y_new,abs(picture_interp(:,:,kk)),channel_id_pixel(1),...
            channel_id_pixel(2));
        end
            
    end
    
    
    
    
    %% Evaulate the results
    %Check if the fiterror is over 1
    
    
    
    
    %Create a table with fitresults (for better readability)
    
    
    fitresult_table = array2table(fitresult,...
        'VariableNames',{'Amplitude', 'Angle', 'Size_x', 'Size_y', 'X0', 'Y0','Z0'});
    %Save to table
    Receptive_field_overview(ii).fitresult = fitresult_table;
    Receptive_field_overview(ii).zfit = zfit;
    Receptive_field_overview(ii).fiterr = fiterr;
    Receptive_field_overview(ii).zerr = zerr;
    Receptive_field_overview(ii).resnorm = resnorm;
    Receptive_field_overview(ii).rr = rr;
    Receptive_field_overview(ii).channel_id_pixel = channel_id_pixel;
    Receptive_field_overview(ii).Cell_idx = cell_kernel_overview(ii).cell_idx;
    cell_kernel_overview(ii).Receptive_field_overview = Receptive_field_overview(ii);
    
   
    
    
    
   
    
end


%% Save data
[~] = sf_organizer(stim_idx,savepath,'variable_name','cell_kernel_overview',...
    'variable',cell_kernel_overview,'overwrite',true);

out = 1;


% 
% 
% 
% 
% 
% %%
% colour_string = {'m','b','g','r','k'};
% figure
% axis([-0.5 MdataSize+0.5 -0.5 MdataSize+0.5])
% hold on
% for ii = 1:nr_cells
%     %This could be a potential function by itself
%     detailed_info = cell_kernel_overview(ii).detailed_info;
%     
%         
%     
% 
% 
% 
% 
% for kk = 1:4
%     zfit_test = squeeze(Receptive_field_overview(ii).zfit(:,:,kk));
%     fitr = mean(Receptive_field_overview(ii).fitresult,1,'omitnan');
%     
%     if prod(isnan(fitr)) == 1
%         continue
%     else 
%         el_ax = ellipse(fitr(5),fitr(6),fitr(3),fitr(4), fitr(2)*pi/180);
%         el_ax.Color = colour_string{kk};
%         
%     end
% end
% end
% % figure
% % C = del2(zfit_test);
% % mesh(X_new,Y_new,zfit_test,C) %plot data
% % xdatahr(:,:,1) = X_new;
% % xdatahr(:,:,2) = Y_new;
% % hold on
% % 
% % surface(X_new,Y_new,D2GaussFunctionRot_new(fitr,xdatahr),'EdgeColor','none') %plot fit
% % axis([-0.5 MdataSize+0.5 -0.5 MdataSize+0.5 -noise noise+fitr(1)])
% % alpha(0.2)  
% % hold off
%      
% 
% 
% 
% 
% %% -----Plot profiles----------------
% hf2 = figure;
% set(hf2, 'Position', [20 20 950 900])
% alpha(0)
% profil_ax = subplot(4,4, [5,6,7,9,10,11,13,14,15]);
% imagesc(X_new(1,:),Y_new(:,1)',zfit_test)
% set(gca,'YDir','reverse')
% colormap('jet')
% 
% string1 = ['       Amplitude', '    Angle', '    X-Width','    Y-Width', '    X-Coordinate', '    Y-Coordinate',];
% string3 = ['Fit      ',num2str(fitr(1), '% 100.3f'),'       ',num2str(fitr(2), '% 100.3f'),'      ',num2str(fitr(3), '% 100.3f'),'      ',num2str(fitr(4), '% 100.3f'),'        ',num2str(fitr(5), '% 100.3f'),'          ',num2str(fitr(6), '% 100.3f')];
% 
% text(0,+MdataSize*1.10,string1,'Color','red')
% text(0,+MdataSize*1.15,string3,'Color','red')
% 
% 
% 
% %% -----Calculate cross sections-------------
% % generate points along horizontal axis
% InterpolationMethod = 'nearest';
% m = tan(fitr(2)*pi/180);% Point slope formula
% b = (-m*fitr(5) + fitr(6));
% xvh = 0:MdataSize;
% yvh = xvh*m + b;
% %hPoints = interp2(X_new,Y_new,zfit_test,xvh,yvh,InterpolationMethod);
% % generate points along vertical axis
% mrot = -m;
% brot = (mrot*fitr(6) - fitr(5));
% yvv = 0:MdataSize;
% xvv = yvv*mrot - brot;
% %vPoints = interp2(X_new,Y_new,zfit_test,xvv,yvv,InterpolationMethod);
% 
% hold on % Indicate major and minor axis on plot
% 
% % % plot pints 
% % plot(xvh,yvh,'r.') 
% % plot(xvv,yvv,'g.')
% 
% % plot lins 
% % plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r') 
% % plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'g') 
% ellipse(fitr(5),fitr(6),fitr(3),fitr(4), fitr(2)*pi/180)
% hold off
% axis([-0.5 MdataSize+0.5 -0.5 MdataSize+0.5])
% 
% 
% %%
% 
% ymin = - noise * fitr(1);
% ymax = fitr(1)*(1+noise);
% xdatafit = linspace(-0.5,MdataSize+0.5,length(X_new));
% hdatafit = fitr(1)*exp(-(xdatafit-fitr(5)).^2/(2*fitr(3)^2));
% vdatafit = fitr(1)*exp(-(xdatafit-fitr(6)).^2/(2*fitr(4)^2));
% cross_ax1 = subplot(4,4, [1:3]);
% %xposh = (xvh-fitr(5))/cos(fitr(2)*pi/180)+fitr(5);% correct for the longer diagonal if fi~=0
% % plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black')
% plot(xdatafit,hdatafit,'black')
% axis([-0.5 MdataSize+0.5 ymin*1.1 ymax*1.1])
% cross_ax2 = subplot(4,4,[8,12,16]);
% % xposv = (yvv-fitr(6))/cos(fitr(2)*pi/180)+fitr(6);% correct for the longer diagonal if fi~=0
% % plot(vPoints,xposv,'g.',vdatafit,xdatafit,'black')
% plot(vdatafit,xdatafit,'black')
% axis([ymin*1.1 ymax*1.1 -0.5 MdataSize+0.5])
% set(gca,'YDir','reverse')
% figure(gcf) % bring current figure to front
% 
% linkaxes([profil_ax,cross_ax1],'x')
% linkaxes([profil_ax,cross_ax2],'y')
%     end
% 
% end
% 
end
%     
% 
% 
% 
% 
% 
