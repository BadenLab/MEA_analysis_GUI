function spiketimestamps = load_spiketimestamps_app (savepath, add_info, varargin)
%This function returns the spiketimestamps for a given stimulus and a
%either all cells in the recording or a subset. 
%Inputs:
%savepath the location of the mat file with the data
%structure wich .stim_begin and .stim_end
%cell subset (returns only a subset of cells)

%Manage Input arguments
p = inputParser;
validMatrixShape = @(x) isnumeric(x) && (size(x,1) ==1) && (size(x,2) == 2);
addRequired(p,'savepath',@ischar);
addRequired(p,'add_info',@isstruct);
addOptional(p,'cell_subset',validMatrixShape);



parse(p,savepath,add_info,varargin{:});

savepath = p.Results.savepath;
add_info = p.Results.add_info;
stx_temp = p.Results.cell_subset;
%% Main




stim_end = add_info.stim_end;
stim_begin = add_info.stim_begin; %Extract information about stimulus begin
%and stimulus end

%Load spiketimestamps from the .mat file
spiketimestamps_temp = load(savepath,'-mat','spiketimestamps');
spiketimestamps_temp = spiketimestamps_temp.spiketimestamps;
sty = size(spiketimestamps_temp,1);


if isnumeric(stx_temp)
    spiketimestamps_temp = spiketimestamps_temp(:,stx_temp(1):stx_temp(2));
end    

stx = size(spiketimestamps_temp,2);



% spiketimestamps_temp = spiketimestamps_temp(:,1:5000);

spiketimestamps = NaN(sty,stx);
%spiketimestamps_temp(1,1)

%Get the right spiketimestamps
y = zeros(numel(stim_end));
for ii = 1:numel(stim_begin)
    spiketimestamps(:,:,ii) = spiketime_extract_app(spiketimestamps_temp,stim_begin(ii),stim_end(ii));
    y(ii) = lastNaN(spiketimestamps(:,:,ii));
    
end
spiketimestamps = spiketimestamps(1:max(y),:,:);

