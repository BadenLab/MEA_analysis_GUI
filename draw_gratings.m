
function out = draw_gratings (savepath, add_info)
%This script draws the results of the grating analysis on the GUI


%Input Parameters
stim_idx = add_info.stim_idx;
cells = add_info.Grating_info.Cell;
cells
max_draw = 100;

%If too many cells are selected a random subset will be plotted
if length(cells)> max_draw
    too_much = length(cells)-30;
    delete_idx = randperm(length(cells),too_much);
    cells(delete_idx) = [];
end
    


%% Load data
S = load(findfile_app(stim_idx,savepath,'Data_circular.mat'));
Data_circular = S.Data_circular;

S = load(findfile_app(stim_idx,savepath,'dir.mat'));
dir = S.dir;

S = load(findfile_app(stim_idx,savepath,'Grating_QC'));
Grating_QC = S.Grating_QC;
clear S

panel_id = add_info.Grating_info.panel_id;
panel_id2 = add_info.Grating_info.panel_id2;
delete(panel_id.Children);
delete(panel_id2.Children);
dir_radians = Data_circular(1).dir_radians;
diff_dir = diff(dir_radians(1:2));
%% Preallocate and stuff
%Depending on the input all or a subset of cells will be selected
nr_cells = length(cells);
true_idx = zeros(1,nr_cells);
true_idx_QC = zeros(1,nr_cells);
for ii = 1:nr_cells
    try
        true_idx(ii) = find([Data_circular.cell_idx] == cells(ii));
        true_idx_QC(ii) = find([Grating_QC.cell_idx] == cells(ii));
    catch
        true_idx(ii) = NaN;
        true_idx_QC = NaN;
    end
end
true_idx(isnan(true_idx)) = [];

Data_new = Data_circular(true_idx);
Grating_QC_new = Grating_QC(true_idx_QC);

%Plotting

sub_nr = numSubplots(length(Data_new));



for jj = 1:length(Data_new)
  %try
  ax(jj) = subplot(sub_nr(1),sub_nr(2),jj,'Parent',panel_id);
  %ax(jj) = subplot(sub_nr(1),sub_nr(2),jj);
  % compute and plot mean resultant vector length and direction
  mw = max(Data_new(jj).spikes_deg);
  mw_w = max(Data_new(jj).spikes_deg_w);
  
  r = circ_r(dir_radians',Data_new(jj).spikes_deg,diff_dir) * mw;
  r_w = circ_r(dir_radians',Data_new(jj).spikes_deg_w,diff_dir) * mw_w;
  
  phi = circ_mean(dir_radians',Data_new(jj).spikes_deg);
  phi_w = circ_mean(dir_radians',Data_new(jj).spikes_deg_w);
  
  
  hold(ax(jj));
  zm = r*exp(1i*phi');
  zm_w = r_w*exp(1i*phi_w');
  
  plot(ax(jj),[0 real(zm)], [0, imag(zm)],'r','linewidth',1.5)
  plot(ax(jj),[0 real(zm_w)], [0, imag(zm_w)],'b','linewidth',1.5)
  
 
  
  % plot the tuning function of the three neurons 
  polar(ax(jj),[dir_radians dir_radians(1)], ...
      [Data_new(jj).spikes_deg' Data_new(jj).spikes_deg(1)],'k')
  polar(ax(jj),[dir_radians dir_radians(1)], ...
      [Data_new(jj).spikes_deg_w' Data_new(jj).spikes_deg_w(1)],'--k')
 
  
  % draw a unit circle
  zz = exp(1i*linspace(0, 2*pi, 101)) * mw;
  zz_w = exp(1i*linspace(0, 2*pi, 101)) * mw_w;
  
  plot(ax(jj),real(zz),imag(zz),'k:')
  plot(ax(jj),[-mw mw], [0 0], 'k:', [0 0], [-mw mw], 'k:')
  
  QC_temp = Grating_QC_new(jj).QC; 
  title_text = ['Cell ',num2str(Data_new(jj).cell_idx),'; QI: ',num2str(QC_temp)];
  title(ax(jj),title_text)
  
  
  if ~Grating_QC_new(jj).QC_pass
      set(ax(jj),'color',[0.5,0.5,0.5]);
  end
    

  formatSubplot(ax(jj),'ax','square','box','off','lim',[-mw mw -mw mw])
  set(ax(jj),'xtick',[])
  set(ax(jj),'ytick',[])
  hold(ax(jj),'off');
%     catch
%         continue
%     end

end




if jj == 1
nr_dir = length(dir);
nr_repeats = length(dir(1).spike_nr_repeats(1,1,:));

for ii = 1:length(Data_new)

sub_nr = numSubplots(nr_dir);
for kk = 1:nr_dir
    
   %Collect all spikes
   
   ix(kk) = subplot(sub_nr(1),sub_nr(2),kk,'Parent',panel_id2);
   for cc = 1
   
   spikes = squeeze(dir(kk).spike_nr_repeats(:,true_idx,:));
   spikes(spikes == 0) = NaN;
   yvalue = ones(1,length(spikes));
   hold(ix(kk));
   for rr = 1:nr_repeats
       last_idx = lastNaN(spikes(:,rr),1);
       scatter(ix(kk),spikes(1:last_idx,rr),(yvalue(1:last_idx)*0.1*rr)+(cc-1)*yvalue(1:last_idx),'k','.');
      
   end
   end
    ylim(ix(kk),[0,cc+0.05])
end
   hold(ax(jj),'off');
   %linkprop(ix,{'XLim','YLim'})
   %linkaxes(ix,'xy');
end
    
end
out = 1;



end
