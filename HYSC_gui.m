function varargout = HYSC_gui(varargin)
% HYSC_GUI MATLAB code for HYSC_gui.fig
%      HYSC_GUI, by itself, creates a new HYSC_GUI or raises the existing
%      singleton*.
%
%      H = HYSC_GUI returns the handle to a new HYSC_GUI or the handle to
%      the existing singleton*.
%
%      HYSC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYSC_GUI.M with the given input arguments.
%
%      HYSC_GUI('Property','Value',...) creates a new HYSC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HYSC_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HYSC_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HYSC_gui

% Last Modified by GUIDE v2.5 11-May-2018 11:20:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HYSC_gui_OpeningFcn, ...
    'gui_OutputFcn',  @HYSC_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before HYSC_gui is made visible.
function HYSC_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HYSC_gui (see VARARGIN)

% Choose default command line output for HYSC_gui
handles.output = hObject;
% handles.runAdvSet = 0;
% set(handles.geneFiltering,'Value',1);
% set(handles.log,'Value',1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HYSC_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HYSC_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
else
    disp('');
    disp('*************************');
    disp('*  Program  Terminated  *');
    disp('*************************');
    for i=1:65
        varargout{i}=0;
    end;
    close(gcf);
end



function dataPath_Callback(hObject, eventdata, handles)
% hObject    handle to dataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataPath as text
%        str2double(get(hObject,'String')) returns contents of dataPath as a double
handles.dataPath = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dataPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex]  = uigetfile('*.xlsx', 'Select Data file');
set(handles.dataPath,'String',fullfile(PathName, FileName));
% handles.dataPathStr = fullfile(PathName, FileName);
guidata(hObject, handles);

% --- Executes on button press in log.
function log_Callback(hObject, eventdata, handles)
% hObject    handle to log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of log
handles.log.Value = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in geneFiltering.
function geneFiltering_Callback(hObject, eventdata, handles)
% hObject    handle to geneFiltering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of geneFiltering
handles.geneFiltering.Value = get(hObject,'Value');
guidata(hObject,handles);


% function HYSCLayers_Callback(hObject, eventdata, handles)
% % hObject    handle to HYSCLayers (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hints: get(hObject,'String') returns contents of HYSCLayers as text
% %        str2double(get(hObject,'String')) returns contents of HYSCLayers as a double
%
%
% % --- Executes during object creation, after setting all properties.
% function HYSCLayers_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to HYSCLayers (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
%
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

% --- Executes on button press in advSet.
function advSet_Callback(hObject, eventdata, handles)
% hObject    handle to advSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.runAdvSet = 1;
hAdvSet = AdvancedSetting();
handles.AdvSet = guidata(hAdvSet);
guidata(hObject,handles);

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
disp('HYSC Running...');
[~, ~, raw] = xlsread(handles.dataPath.String);
data = cell2mat(raw(2:end, 2:end));
cell_IDs = raw(1, 2:end); % cell IDs
geneList = raw(2:end, 1);
%% change cell_IDs to numbers if not
if ~iscellstr(geneList)
    A = cellfun(@ischar,geneList);
    idx = find(A==0);
    for i = 1:length(idx)
        geneList{idx(i)} = num2str(geneList{idx(i)});
    end
end
if ~isnumeric(cell_IDs{1})
    cell_type_unique = unique(cell_IDs);
    cell_types = zeros(size(cell_IDs));
    for j = 1:length(cell_type_unique)
        idx1 = strcmp(cell_IDs, cell_type_unique{j});
        idx2 = find(idx1==1); % idx2==1, matched cell type
        cell_types(idx2) = j;
    end
else
    cell_types = cell2mat(cell_IDs);
end


%% log transform
if (handles.log.Value == 1)
    log2data = log2(data+1);
else
    log2data = data;
end

%% gene filtering
if (handles.geneFiltering.Value == 1)
    Idx = [];
    for i = 1:size(log2data,1)
        num = length(find(log2data(i,:)~=0));
        if num>str2double(handles.AdvSet.filteringNo.String)
            if var(log2data(i,:))>str2double(handles.AdvSet.filteringVar.String)
                Idx = [Idx, i];
            end
        end
    end
    dataBuffer = log2data(Idx,:);
    geneListBuffer = geneList(Idx);
else
    dataBuffer = log2data;
    geneListBuffer = geneList;
end


if handles.runAdvSet==1 && handles.cellAnno.Value==1
    output = HYSC(dataBuffer, geneListBuffer, ...
        'minClusterSize_number',handles.AdvSet.minClusterSize_number.Value, ...
        'dimClustering',handles.AdvSet.dimClustering.Value, 'r2cutoff',handles.AdvSet.r2cutoff.Value, ...
        'pcutoff',handles.AdvSet.pcutoff.Value, 'perp', handles.AdvSet.Perplexity.Value,...
        'maxHYSCLayer', handles.maxLayers.Value, 'cores',handles.AdvSet.cores.Value, ...
        'tSNEScores', 1:5, 'annoPlot', cell_types);
elseif handles.runAdvSet==1 && handles.cellAnno.Value==0
    output = HYSC(dataBuffer, geneListBuffer, ...
        'minClusterSize_number',handles.AdvSet.minClusterSize_number.Value, ...
        'dimClustering',handles.AdvSet.dimClustering.Value, 'r2cutoff',handles.AdvSet.r2cutoff.Value, ...
        'pcutoff',handles.AdvSet.pcutoff.Value, 'perp', handles.AdvSet.Perplexity.Value,...
        'maxHYSCLayer', handles.maxLayers.Value, 'cores',handles.AdvSet.cores.Value, ...
        'tSNEScores', 1:5);
else
    output = HYSC(dataBuffer, geneListBuffer);
end
% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
set(gcf,'visible','off');
quit;

% --- Executes on button press in cellAnno.
function cellAnno_Callback(hObject, eventdata, handles)
% hObject    handle to cellAnno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cellAnno
handles.cellAnno.Value = get(hObject,'Value');
guidata(hObject,handles);



function maxLayers_Callback(hObject, eventdata, handles)
% hObject    handle to maxLayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxLayers as text
%        str2double(get(hObject,'String')) returns contents of maxLayers as a double
handles.maxLayers.Value = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function maxLayers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxLayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function output = HYSC(data, gene_symbol, varargin)

% run Hybrid-clustering for Single-cell Categrization (HYSC)
%
% varargin:
% Input parameters:
% ClusterNum               - the given number of clusters. If ClusterNum = 0, the algorithm automatically estimates the number of clusters
% Nmax                     - the maximum number of cell clusters. By default, Nmax = 50
% minClusterSize_number    - the minimum number of samples (cells) in each cluster. By default, minClusterSize_number = 10
% dimClustering            - the number of principle components adopted by k-means clustering. By default, dimClustering = 30
% r2cutoff                 - OPLSDA parameters, r2 cutoff value to identify discriminatory variables. By default 0.5
% pcutoff                  - OPLSDA parameters, p cutoff value to identify discriminatory variables. By default 0.05 after bonferroni correction
% perp                     - perplexity parameter of tSNE. By default, perp = 30
% maxHYSCLayer             - the maximum number of HYSC layers. By default, maxHYSCLayer = 5
% cores                    - the number of cores used in parallel computation. By default, cores = 8
% tSNEScores               - the index of tSNE scores adopted for scores plot. By default, the first two components are illustrated.
% annoPlot                 - priori information, e.g., time information of cells, to be incorporated in the output plot. By default, []
%
% Output:
% 'HYSC Output.xlsx':
% Sheet 1 - the cell categorization indexs
% Sheet 2 - tSNE scores for scatter plot
% Sheet 3 - gene expression data
% 'HYSC markerGenes.xlsx': each sheet contains the gene markers of a cell cluster
%
% Author Jie Hao & Xin Zou, SJTU, China
% Copyright Jie Hao & Xin Zou

warning off

varList = {'ClusterNum', 'Nmax', 'minClusterSize_number',...
    'dimClustering', 'nc', 'ncox', 'ncoy', 'r2cutoff', 'pcutoff', 'preprocessing',...
    'perp', 'maxHYSCLayer', 'cores', 'tSNEScores', 'annoPlot'};
convgThrh = {0, 50, 10, 5,  1, 1, 0, 0.5, 0.05, 'uv', 30, 5, 8, 1:2, []}; % default values
if ~isempty(varargin)
    L = length(varargin);
    for i = 1:2:L-1
        idx = find(strcmp(varargin{i}, varList));
        convgThrh{idx} = varargin{i+1};
    end
end

ClusterNum = convgThrh{1}; % 'maxClusterNum'
Nmax = convgThrh{2}; % 'Nmax'
minClusterSize_number = convgThrh{3}; % 'minClusterSize_number'
dim = convgThrh{4}; % 'dimClustering'
oplsda_para.nc = convgThrh{5};
oplsda_para.ncox = convgThrh{6};
oplsda_para.ncoy = convgThrh{7};
oplsda_para.r2cutoff = convgThrh{8};
oplsda_para.pcutoff = convgThrh{9};
oplsda_para.preprocessing = convgThrh{10};
perp = convgThrh{11}; % tSNE perplexity
maxHYSCLayer = convgThrh{12}; % maximum number of HYSC layers
cores = convgThrh{13}; % for parallel computation
tSNEScores = convgThrh{14}; % index of scores for plot demonstration
annoPlot = convgThrh{15}; % priori information
oplsda_para.NP = 0; % by default no permutation test is performed

%% HYSC parameters
paraList_HYSC = {'ClusterNum',ClusterNum, 'Nmax',Nmax, 'minClusterSize_number', minClusterSize_number,...
    'dimClustering', dim, 'nc', oplsda_para.nc, 'ncox', oplsda_para.ncox, ...
    'ncoy', oplsda_para.ncoy, 'r2cutoff',oplsda_para.r2cutoff...
    'pcutoff', oplsda_para.pcutoff, 'preprocessing', oplsda_para.preprocessing, 'perp', perp};

% HYSC cluster symbles
types = cell(1,26);
s=(1:26)+96;
str=char(s);
for i = 1:length(str)
    types{i} = str(i);
end

clusterMarkParal = cell(1,cores);
unsepClusterMarkParal = cell(1,cores);
%%
startmatlabpool(cores);
biomarkers = cell(1, cores);
parfor paral = 1:cores
    % for paral = 1:cores
    clusterMark = cell(size(data, 2), 1);
    clusterMark(:) = cellstr(types{1});
    unsepClusterMark = [];
    while 1
        uniqueClusterMark = unique(clusterMark);
        NumOfClusters = length(uniqueClusterMark);
        L = zeros(1,NumOfClusters);
        for i = 1:NumOfClusters
            L(i) = length(char(uniqueClusterMark{i}));
        end
        HYSCLayers = max(L); %  current HYSC layer
        if NumOfClusters < Nmax && HYSCLayers <= maxHYSCLayer
            for n = 1:NumOfClusters
                if ismember(uniqueClusterMark{n}, unsepClusterMark) == 0
                    idx = find(ismember(clusterMark, uniqueClusterMark{n}));
                    if length(idx) >= minClusterSize_number % the size of each cell cluster should be no smaller than the pre-defined minimum cluster size
                        dataBuffer = data(:, idx);
                        temp = sum(dataBuffer~=0,2);
                        Idx = find(temp>2);
                        dataBuffer = dataBuffer(Idx,:);
                        outputHYSC = HYSC_s(dataBuffer, paraList_HYSC); % HYSC for single layer
                        if max(outputHYSC.groupIdx) == 1 % the current subset is unsplittable
                            unsepClusterMark = [unsepClusterMark, uniqueClusterMark(n)];
                        else
                            for j = 1:max(outputHYSC.groupIdx)
                                idx2 = find(outputHYSC.groupIdx==j);
                                for k = 1:length(idx2)
                                    clusterMark{idx(idx2(k))} = [clusterMark{idx(idx2(k))}, types{j}];
                                end
                            end
                            biomarkers{paral} = [biomarkers{paral}; Idx(outputHYSC.biomarkers)];
                        end
                    else
                        unsepClusterMark = [unsepClusterMark, uniqueClusterMark(n)];
                    end
                end
            end
        end
        if length(unsepClusterMark) == NumOfClusters && NumOfClusters == 1
            break;
        end
        if ~isempty(unsepClusterMark)  && length(unsepClusterMark) == NumOfClusters || NumOfClusters >= Nmax || HYSCLayers == maxHYSCLayer
            break;
        end
    end
    unsepClusterMarkParal{paral} = unsepClusterMark;
    clusterMarkParal{paral} = clusterMark;
    biomarkers{paral} = unique(biomarkers{paral});
end


NumOfClusters_paral = zeros(1,cores);
clusterIdx = cell(1,cores);
uniqueClusterMark = cell(1,cores);
for paral = 1:cores
    uniqueClusterMark{paral} = unique(clusterMarkParal{paral});
    NumOfClusters_paral(paral) = length(uniqueClusterMark{paral});
    clusterIdx{paral} = cell(1,NumOfClusters_paral(paral));
    for j = 1:NumOfClusters_paral(paral)
        idx = find(strcmp(uniqueClusterMark{paral}{j},clusterMarkParal{paral}));
        clusterIdx{paral}{j} = idx;
    end
end

[freqs, centers] = hist(NumOfClusters_paral, unique(NumOfClusters_paral));
[~, IX] = max(freqs);
idx = find(NumOfClusters_paral==centers(IX));
idx = idx(1);


clusterMark = clusterMarkParal{idx};
unsepClusterMark = unsepClusterMarkParal{idx};
NumOfClusters = NumOfClusters_paral(idx);
biomarkers_final = biomarkers{idx};

closematlabpool;

if length(unsepClusterMark) == NumOfClusters && NumOfClusters == 1
    output.biomarkers = [];
    output.clusterIdx = ones(size(data, 2), 1);
else
    uniqueClusterMark = unique(clusterMark);
    NumOfClusters = length(uniqueClusterMark);
    biomarkers_idx = biomarkers_final;
    IDX = zeros(size(clusterMark));
    
    for n = 1:NumOfClusters
        IDX(ismember(clusterMark, uniqueClusterMark{n})) = n;
    end
    
    
    output.clusterIdx = IDX;
    output.biomarkers = biomarkers_final;

    
    
    % identify gene markers for each cluster and sort by p-values
    markers_idx = cell(1,max(IDX));
    markers_genes = cell(1,max(IDX));
    markers_genes_list = [];
    Temp = [];
    sorted_gene_symbol =[];
    for i = 1:max(IDX)
        tempData1 = data(:,IDX==i)';
        tempData2 = data(:,IDX~=i)';
        model = OPLSDA(tempData1, tempData2, 1:size(data,1), oplsda_para);
        Sign = sign(mean(data(model.sig_idx,IDX==i),2)-mean(data(model.sig_idx,IDX~=i),2));
        markers_idx{i} = model.sig_idx(Sign>0);
        [B, IX] = sort(model.P(markers_idx{i}),'ascend');
        markers_idx{i} = markers_idx{i}(IX);
        markers_genes{i} = gene_symbol(markers_idx{i});
        markers_genes_list = [markers_genes_list; markers_genes{i}];
    end
    markers_genes_list = unique(markers_genes_list);
    markers_genes_xls = markers_genes;
    
    % remove duplicate gene markers between cell clusters
    for i = 1:length(markers_genes_list)
        Locb = [];
        LocGroup = [];
        for j = 1:max(IDX)
            [~,locTemp] = ismember(markers_genes_list{i}, markers_genes{j});
            if locTemp ~= 0
                Locb = [Locb, locTemp];
                LocGroup = [LocGroup, j];
            end
        end
        if length(Locb)>1 % duplicate genes
            mu = zeros(1,length(Locb));
            for k = 1:length(Locb)
                mu(k) = mean(data(markers_idx{LocGroup(k)}(Locb(k)),IDX==LocGroup(k)));
            end
            IX = find(mu<max(mu));
            for k = 1:length(IX) % remove duplicate genes
                markers_idx{LocGroup(IX(k))}(Locb(IX(k))) = [];
                markers_genes{LocGroup(IX(k))}(Locb(IX(k))) = [];
            end
        end
    end
    
    % construct expression profiling data based on the gene markers
    for i = 1:max(IDX)
        IX = 1:length(markers_idx{i});
        Temp = [Temp; data(markers_idx{i}(IX(1:min(floor(70/max(IDX)), length(IX)))),:)];
        sorted_gene_symbol = [sorted_gene_symbol; markers_genes{i}(IX(1:min(floor(70/max(IDX)), length(IX))))];
        if i ~= max(IDX)
            Temp = [Temp; NaN(1,size(Temp,2))];
            sorted_gene_symbol = [sorted_gene_symbol; {''}];
        end
    end
    Temp2 = [];
    for i = 1:max(IDX)
        Temp2 = [Temp2, Temp(:, IDX==i)];
        if i ~= max(IDX)
            Temp2 = [Temp2, NaN(size(Temp2,1),1)];
        end
    end
    figure;
    map = colormap(redbluecmap);
    map = map(3:10,:);
    colormap(map);
    figHandle1 = imagesc(Temp2);
    yticks(1:length(sorted_gene_symbol));
    yticklabels(sorted_gene_symbol);
    set(gca,'xTickLabel','');
    
    
    if isempty(annoPlot) % demonstrate HYSC categorization result if annoPlot is not provided
        for i = 1:max(IDX)
            annoPlot = [annoPlot, i*ones(1, length(find(IDX==i)))];
            if i ~= max(IDX)
                annoPlot = [annoPlot, NaN(size(annoPlot,1),1)];
            end
        end
    else % if annoPlot is provided
        Temp3 = []; % for annoPlot
        for i = 1:max(IDX)
            Temp3 = [Temp3, annoPlot(find(IDX==i))];
            if i ~= max(IDX)
                Temp3 = [Temp3, NaN(size(Temp3,1),1)];
            end
        end
        annoPlot = Temp3;
    end
    figHandle2 = subplot('Position',[0.13 0.934 0.775 0.01]);
    hold on;
    
    Color = [0,0,0.5625;0,0,0.625;0,0,0.6875;0,0,0.7500;0,0,0.8125;0,0,0.8750;0,0,0.9375;...
        0,0,1;0,0.06250,1;0,0.1250,1;0,0.1875,1;0,0.2500,1;0,0.3125,1;0,0.3750,1;0,0.4375,1;...
        0,0.5000,1;0,0.5625,1;0,0.6250,1;0,0.6875,1;0,0.7500,1;0,0.8125,1;0,0.8750,1;0,0.9375,1;...
        0,1,1;0.06250,1,0.9375;0.1250,1,0.8750;0.1875,1,0.8125;0.2500,1,0.7500;0.3125,1,0.6875;...
        0.3750,1,0.6250;0.4375,1,0.5625;0.5000,1,0.5000;0.5625,1,0.4375;0.6250,1,0.3750;0.6875,1,0.3125;...
        0.7500,1,0.2500;0.8125,1,0.1875;0.8750,1,0.1250;0.9375,1,0.06250;1,1,0;1,0.9375,0;1,0.8750,0;...
        1,0.8125,0;1,0.7500,0;1,0.6875,0;1,0.6250,0;1,0.5625,0;1,0.5000,0;1,0.4375,0;1,0.3750,0;1,0.3125,0;...
        1,0.2500,0;1,0.1875,0;1,0.1250,0;1,0.06250,0;1,0,0;0.9375,0,0;0.8750,0,0;0.8125,0,0;0.7500,0,0;0.6875,0,0;...
        0.6250,0,0;0.5625,0,0;0.5000,0,0];
    Step_size = max(annoPlot);
    Step = floor(size(Color,1)/Step_size);
    color_code = 1:Step:Step*Step_size;
    xscale = 1:length(annoPlot);
    for i = 1:length(annoPlot)
        if ~isnan(annoPlot(i))
            plot(xscale(i),0,'.', 'MarkerEdgeColor', Color(color_code(annoPlot(i)),:), 'MarkerFaceColor', Color(color_code(annoPlot(i)),:),  'MarkerSize',15, 'Parent', figHandle2);
        else
            plot(xscale(i),0,'.', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w',  'MarkerSize',5, 'Parent', figHandle2);
        end
    end
    xlim([0.5 length(annoPlot)+0.5]);
    axis off;
    
    
    subplot('Position',[0.13 0.08 0.775 0.02]);
    xlim([0.5 size(Temp2,2)+0.5]);
    cate_IDX = [];
    for i = 1:max(IDX)
        cate_IDX = [cate_IDX, i*ones(1, length(find(IDX==i)))];
        if i ~= max(IDX)
            cate_IDX = [cate_IDX, NaN(size(cate_IDX,1),1)];
        end
    end
    idx = find(isnan(cate_IDX));
    idx = [1, idx, length(cate_IDX)];
    centres = zeros(1, length(idx)-1);
    for i = 1:length(idx)-1
        centres(i) = (idx(i+1)-idx(i))/2+idx(i);
        t = text(centres(i), 0.5, num2str(i));
        t.FontSize = 15;
    end
    axis off;
    
    subplot('Position',[0.92 0.111 0.025 0.813]);
    Step_size = max(annoPlot);
    Step = 10/Step_size;
    Ys = 0.5:Step:10.5;
    for i = 1:max(annoPlot)
        plot(0, Ys(i), 'o', 'markersize', 8, 'MarkerFaceColor', Color(color_code(i),:), 'MarkerEdgeColor', Color(color_code(i),:));
        t = text(1.5, Ys(i), num2str(i));
        t.FontSize = 15;
        if i == 1
            hold on;
        end
    end
    ylim([0,11]);
    axis off
    
    
    delete('HYSC gene markers.xlsx');
    for i = 1:length(markers_genes_xls)
        if ~isempty(markers_genes_xls{i})
            xlswrite('HYSC gene markers.xlsx', markers_genes_xls{i}, i);
        else
            xlswrite('HYSC gene markers.xlsx', {NaN}, i);
        end
    end
    e = actxserver('Excel.Application'); % # open Activex server
    ewb = e.Workbooks.Open([pwd, '\HYSC gene markers.xlsx']); % # open file (enter full path!)
    for i = 1:length(markers_genes_xls)
        ewb.Worksheets.Item(i).Name = ['Cluster', num2str(i)]; % # rename 1st sheet
    end
    ewb.Save % # save to the same file
    ewb.Close(false)
    e.Quit
    
    
    %% tSNE scatter plot
    X = data(biomarkers_idx,:)';
    [COEFF,SCORE,latent] = pca(X);
    if size(X,1)<5000 % to save execution time
        ydata = tsne(X,[],SCORE(:,1:min(5,size(X,2))), perp);
    else
        ydata = SCORE(:,1:min(5,size(X,2)));
    end
    clusteringCentroids = zeros(max(IDX), min(5,size(X,2)));
    for i = 1:max(IDX)
        ID = find(IDX==i);
        clusteringCentroids(i,:) = mean(ydata(ID,:));
    end
    output.tSNE = ydata;
    output.clusteringCentroids = clusteringCentroids;
    
    figure;
    subplot('Position',[0.14 0.16 0.73 0.76]);
    colormap(jet);
    Color=get(gcf,'Colormap');
    
    % color code according to HYSC hierarchy
    layerInfor = cellfun('length',clusterMark)-1;
    cate_layerd = zeros(length(layerInfor), max(layerInfor));
    for j = 1:max(layerInfor) % for each layer
        buffer = cell(1,length(layerInfor));
        for jj = 1:length(layerInfor)
            if layerInfor(jj)>=j % only apply for the cluters obtained by current layer or subsequent layers
                buffer{jj} = clusterMark{jj}(j+1); % exclude the firt pseudo layer symbol
            end
        end
        idx1 = find(~cellfun('isempty', buffer)); % identify the index of none empty 'buffer' element
        buffer_unique = unique(buffer(idx1));
        temp_idx = cell(length(buffer_unique));
        for jj = 1:length(buffer_unique)
            idx2 = find(strcmp(buffer_unique{jj}, buffer(idx1)));
            cate_layerd(idx1(idx2),j) = jj;
        end
    end
    % assign the color to each sample
    color_cate_layerd = cell(size(cate_layerd,1), min(max(layerInfor),3));
    for j = 1:min(max(layerInfor),3) % for each layer
        if j == 1 % the first layer
            buffer = unique(cate_layerd(:,j));
            linIdx = floor(linspace(1,64,length(buffer)+1)); % linearly spaced index vector
            for jj = 1:length(buffer)
                color_cate_layerd(find(cate_layerd(:,j)==buffer(jj)),j) = {[linIdx(jj),linIdx(jj+1)]};
            end
        elseif j == 2 % the subsequent layers
            buffer = unique(cate_layerd(:,j-1)); % for each upper-layer
            for jj = 1:length(buffer)
                idx = find(cate_layerd(:,j-1)==buffer(jj)); % check if the upper layer has subsequent layers
                if sum(cate_layerd(idx,j))>0 % there is a subsequent layer
                    buffer2 = unique(cate_layerd(idx,j));
                    linIdx = floor(linspace(color_cate_layerd{idx(1), j-1}(1), color_cate_layerd{idx(1), j-1}(2), length(buffer2)+1)); % linearly spaced index vector according the color range of the upper layer
                    for jjj = 1:length(buffer2)
                        idx2 = find(cate_layerd(idx,j)==buffer2(jjj));
                        color_cate_layerd(idx(idx2),j) = {[linIdx(jjj),linIdx(jjj+1)]};
                    end
                end
            end
        elseif j == 3
            buffer1 = unique(cate_layerd(:,1)); % for each upper-layer
            buffer2 = unique(cate_layerd(:,2)); % for each upper-layer
            for jj = 1:length(buffer1)
                for jjj = 1:length(buffer2)
                    idx = find(cate_layerd(:,1)==buffer1(jj) & cate_layerd(:,2)==buffer2(jjj)); % check if the upper layer has subsequent layers
                    
                    if sum(cate_layerd(idx,j))>0 % there is a subsequent layer
                        buffer3 = unique(cate_layerd(idx,j));
                        linIdx = floor(linspace(color_cate_layerd{idx(1), j-1}(1), color_cate_layerd{idx(1), j-1}(2), length(buffer3)+1)); % linearly spaced index vector according the color range of the upper layer
                        for j4 = 1:length(buffer3)
                            idx2 = find(cate_layerd(idx,j)==buffer3(j4));
                            color_cate_layerd(idx(idx2),j) = {[linIdx(j4),linIdx(j4+1)]};
                        end
                    end
                    
                end
            end
        end
        
    end
    
    for i = 1:length(IDX)
        idx = max(find(~cellfun('isempty', color_cate_layerd(i,:))));
        color_code(i) = round(mean(color_cate_layerd{i,idx}));
        plot(ydata(i,tSNEScores(1)), ydata(i,tSNEScores(2)), '.', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 0.5);
        t = text(ydata(i,tSNEScores(1)), ydata(i,tSNEScores(2)),num2str(IDX(i)));
        t.Color = Color(color_code(i),:);
        t.FontSize = 15;
        if i == 1
            hold on;
        end
    end
    set(gca, 'fontsize', 14);
    xlabel('tSNE_1');
    ylabel('tSNE_2');
    
%     subplot('Position',[0.9 0.16 0.05 0.76]);
%     imagesc(1,1:-1/64:0,(1:64)');
%     set(gca,'xTickLabel','');
%     set(gca,'yTickLabel','');
    
    delete('HYSC Output.xlsx');
    xlswrite('HYSC Output.xlsx', output.clusterIdx, 1);
    xlswrite('HYSC Output.xlsx', output.tSNE, 2);
    temp = [gene_symbol, num2cell(data)];
    xlswrite('HYSC Output.xlsx', temp, 3);
    e = actxserver('Excel.Application'); % # open Activex server
    ewb = e.Workbooks.Open([pwd, '\HYSC Output.xlsx']); % # open file (enter full path!)
    ewb.Worksheets.Item(1).Name = 'Cluster indexs'; % # rename 1st sheet
    ewb.Worksheets.Item(2).Name = 'tSNE scores'; % # rename 2nd sheet
    ewb.Worksheets.Item(3).Name = 'Expr data'; % # rename 3rd sheet
    ewb.Save % # save to the same file
    ewb.Close(false)
    e.Quit
end




function model = OPLSDA(data1,data2, ppm, oplsda_para)

% run OPLS-DA and display results
% Input:
% "data1" and "data2" are groups of spectra from two biological classes
% "ppm" is the ppm variables

% Output:
% "model" contains all OPLS-DA results, where the discriminatory variables are in the model.sig_idx.
% To extract the ppm of biomarker variables, type in
% 'biomarker_ppm=ppm(model.sig_idx)' in the command window.

% The parameters are set in the OPLSDA para.txt file
% nrcv - number of folds cross-validation
% nc - number of correlated variables in X
% ncox - number of orthogonal variables in X
% ncoy - number of orthogonal variables in Y
% r2cutoff - r2 cutoff value to identify discriminatory variables
% p-cutoff - p cutoff value to identify discriminatory variables
% np - number of permutation tests
% preprocessing - pre processing methods, 'uv', 'pa', ''mc'

nc = oplsda_para.nc;
ncox = oplsda_para.ncox;
ncoy = oplsda_para.ncoy;
r2cutoff = oplsda_para.r2cutoff;
pCutoff = oplsda_para.pcutoff;
prep = oplsda_para.preprocessing;
pCutoff = pCutoff/length(ppm);
NP = oplsda_para.NP;

%% format data
X = [data1;data2];
Y = zeros(size(X,1),2);
Y(1:size(data1),1) = 1;
Y(size(data1)+1:end,2) = 1;

%%
if NP == 0
    model=o2pls_m(X,Y,prep,nc,ncox,ncoy,0, r2cutoff, pCutoff);
    if isempty(model)
        return;
    end
else
    
end



function model=o2pls_m(X,Y,prep,nc,ncox,ncoy,nrcv, r2cutoff, pCutoff)
%
% O2-PLS function used to build the model
% For further predictions use the o2pls_pred function
%
% Input:
% X - X matrix input
% Y - Y matrix input
% prep - preprocessing available
%           - 'no' : no preprocessing
%           - 'mc' : meancentering
%           - 'uv' : Univariance scaling
%
% nc - number of correlated components, it can be determinated by PCA of Y'X
% ncox - number of orthogonal components in X
% ncoy - number of orthogonal component in Y
% nrcv - number of fold in the cross-validation (full cross validation)
% v1.0, Xin Zou, Kent

opt = 1;
BT = 0;

if nargin<8
    r2cutoff=0.3;
    pCutoff=0.05;
end


muClass1 = mean(X((Y(:,1)==1),:));
muClass2 = mean(X((Y(:,2)==1),:));
Signs = sign(muClass2-muClass1);

if nargin < 7
    nrcv = 0; % no cross validation
end

[nsx,nvx]=size(X);
[nsy,nvy]=size(Y);

model.ns = nsx;

model.nc=nc;
model.ncox=ncox;
model.ncoy=ncoy;

if nsx~=nsy
    disp('Number of samples are different in X and Y')
    model=[];
    return
end

model.MeanX=mean(X);
model.MeanY=mean(Y);

model.StdX=std(X);
model.StdY=std(Y);

model.SSX=sum(sum(X.^2));
model.SSY=sum(sum(Y.^2));
model.CSX=sum(X.^2);
model.CSY=sum(Y.^2);
% Preprocessing

model.modelType = 'pls';



if nrcv == 0 % no CV
    switch prep
        
        case 'no'
            model.preprocessing='no';
            
        case 'mc'
            %         disp('Meancentering')
            model.preprocessing='mc';
            X=X-repmat(model.MeanX,nsx,1);
            Y=Y-repmat(model.MeanY,nsy,1);
            
        case 'uv'
            % disp('Univariance scaling')
            model.preprocessing='uv';
            X=X-repmat(model.MeanX,nsx,1);
            Y=Y-repmat(model.MeanY,nsy,1);
            X=X./repmat(model.StdX,nsx,1);
            Y=Y./repmat(model.StdY,nsy,1);
            
        case 'pa'
            model.preprocessing='pa';
            X=X-repmat(model.MeanX,nsx,1);
            Y=Y-repmat(model.MeanY,nsy,1);
            X=X./repmat(sqrt(model.StdX),nsx,1);
            Y=Y./repmat(sqrt(model.StdY),nsy,1);
            
        otherwise
            disp('Unknown Preprocessing')
            model=[];
            return
            
    end
    X(find(isnan(X))) = 0;
    M=O2pls(X,Y,nc,ncox,ncoy,'standard');
    if isempty(M)
        return;
    end
    M.preprocessing = model.preprocessing;
    %     M.W = M.W;
    M.ns = model.ns;
    
    M.nc = model.nc;
    M.ncox = model.ncox;
    M.ncoy = model.ncoy;
    
    M.MeanX = model.MeanX;
    M.MeanY = model.MeanY;
    
    M.StdX = model.StdX;
    M.StdY = model.StdY;
    M.X_out = X;
    
    M.SSX = model.SSX;
    M.SSY = model.SSY;
    
    M.CSX = model.CSX;
    M.CSY = model.CSY;
    M.Q2Yhatcum = [];
    model = M;
    
    [~,nv]=size(model.W);
    switch prep
        case 'uv'
            WC = model.W';
        case 'pa'
            for k = 1:nv
                WC(k,:) = model.W(:,k)'./sqrt(model.StdX);
            end
        case 'mc'
            for k = 1:nv
                WC(k,:) = model.W(:,k)'./(model.StdX);
            end
    end
    Bool = exist('BT', 'var');
    if Bool == 1 && ~isempty(BT) && BT == 1 % bootstrapping
        
        f = @(x,y)o2pls_bootstrap(x,y,model.preprocessing,1,ncox,0,0);
        % % bootstrap
        [ci,bootstat] = bootci(500,{f,X,Y}, 'type', 'norm'); % bootstrap
        R = bootstat(:,1:size(X,2));
        median_R = median(R);
        [~,IX_R] = sort(abs(median_R),'descend');
        P = bootstat(:,size(X,2)+1:2*size(X,2)); % p-value
        median_P = median(P);
        W2_bt = bootstat(:,2*size(X,2)+1:3*size(X,2)); % W2
        median_W2_bt = median(W2_bt);
        if ~isempty(pCutoff)
            IX_P = find(median_P<pCutoff); %
        end
        
        
    else  % without bootstrapping
        [nsx,nvx]=size(X);
        buffer = repmat(model.T(:,1),[1,nvx]);
        CC = mean((buffer-repmat(mean(buffer),nsx,1)).*(X-repmat(mean(X),nsx,1)))/(nsx-1)*nsx;
        s1 = std(model.T(:,1));
        S1 = repmat(s1,[1,nvx]);
        S2 = std(X);
        median_R = CC./(S1.*S2);
        median_R(isnan(median_R))=0;
        median_P = tTest(median_R,nsx);
        if ~isempty(pCutoff)
            IX_P = find(median_P<pCutoff); %
        end
        
        switch prep
            case 'uv'
                %                 WC = model.W';
                W2=WC.*WC;
                for k = 1:nv
                    W2(k,:) = W2(k,:)/norm(W2(k,:), 'inf');
                end
            case 'pa'
                for k = 1:nv
                    %                     WC(k,:) = model.W(:,k)'./sqrt(model.StdX);
                    W2(k,:)=WC(k,:).*WC(k,:);
                    W2(k,:) = W2(k,:)/norm(W2(k,:), 'inf');
                end
            case 'mc'
                for k = 1:nv
                    %                     WC(k,:) = model.W(:,k)'./(model.StdX);
                    W2(k,:)=WC(k,:).*WC(k,:);
                    W2(k,:) = W2(k,:)/norm(W2(k,:), 'inf');
                end
        end
        median_W2_bt = W2(1,:);
    end
    
    Bool = exist('opt', 'var');
    if ~isempty(pCutoff) && ~isempty(r2cutoff)
        if Bool == 1 && ~isempty(opt) && opt == 0
            idx = find(median_W2_bt>r2cutoff);
        else
            II = find(median_W2_bt>r2cutoff);
            idx = intersect(IX_P, II);
        end
        
    else
        idx = [];
    end
    
    % % save output
    model.sig_idx = idx; % the significant ppm points
    model.R = median_R;
    model.P = median_P;
    %     model.WC_unscaled = WC.*model.StdX;
    model.W2 = median_W2_bt;
    model.WC = WC;
    model.signs = Signs;
    
    
else % CV
    
    % nrcv-fold cross validation
    
    block_num = floor(nsx/nrcv);
    %     Q2Yhatcum = zeros(1,nrcv);
    
    for cv = 1:nrcv % nrcv iterations of CV
        idx_test = nrcv*([1:block_num]-1)+cv;
        idx_tr = 1:nsx;
        idx_tr(idx_test) = [];
        X_test = X(idx_test,:);
        Y_test = Y(idx_test,:);
        X_tr = X(idx_tr,:);
        Y_tr = Y(idx_tr,:);
        
        [nsx_tr,nvx_tr]=size(X_tr);
        [nsy_tr,nvy_tr]=size(Y_tr);
        [nsx_test,nvx_test]=size(X_test);
        [nsy_test,nvy_test]=size(Y_test);
        switch prep
            
            case 'no'
                model.preprocessing='no';
                
            case 'mc'
                %                 disp('Meancentering')
                model.preprocessing='mc';
                X_tr=X_tr-repmat(mean(X_tr),nsx_tr,1);
                Y_tr=Y_tr-repmat(mean(Y_tr),nsy_tr,1);
                X_test=X_test-repmat(mean(X_test),nsx_test,1);
                Y_test=Y_test-repmat(mean(Y_test),nsy_test,1);
            case 'uv'
                % disp('Univariance scaling')
                model.preprocessing='uv';
                X_tr=X_tr-repmat(mean(X_tr),nsx_tr,1);
                Y_tr=Y_tr-repmat(mean(Y_tr),nsy_tr,1);
                X_tr=X_tr./repmat(std(X_tr),nsx_tr,1);
                Y_tr=Y_tr./repmat(std(Y_tr),nsy_tr,1);
                
                X_test=X_test-repmat(mean(X_test),nsx_test,1);
                Y_test=Y_test-repmat(mean(Y_test),nsy_test,1);
                X_test=X_test./repmat(std(X_test),nsx_test,1);
                Y_test=Y_test./repmat(std(Y_test),nsy_test,1);
            case 'pa'
                model.preprocessing='pa';
                X_tr=X_tr-repmat(mean(X_tr),nsx_tr,1);
                Y_tr=Y_tr-repmat(mean(Y_tr),nsy_tr,1);
                X_tr=X_tr./repmat(sqrt(std(X_tr)),nsx_tr,1);
                Y_tr=Y_tr./repmat(sqrt(std(Y_tr)),nsy_tr,1);
                
                X_test=X_test-repmat(mean(X_test),nsx_test,1);
                Y_test=Y_test-repmat(mean(Y_test),nsy_test,1);
                X_test=X_test./repmat(sqrt(std(X_test)),nsx_test,1);
                Y_test=Y_test./repmat(sqrt(std(Y_test)),nsy_test,1);
                
            otherwise
                disp('Unknown Preprocessing')
                model=[];
                return
                
        end
        
        % training
        model=O2pls(X_tr,Y_tr,nc,ncox,ncoy,'standard');
        if isempty(model)
            return;
        end
        %yhat
        modelPredy=O2plsPred(X_test, Y_test, model,ncox,ncoy,'x');
        %         %xhat
        %         modelPredx=O2plsPred(cvSet.xTest, cvSet.yTest, model,ioax-1,ioay-1,'y');
        
        % % the overall Q2
        SSY=sum(sum(Y_test.^2));
        Q2Yhatcum(cv) = 1-sum(sum((modelPredy.Yhat-Y_test).*(modelPredy.Yhat-Y_test)))/SSY;
    end
    model.Q2Yhatcum = mean(Q2Yhatcum(~isnan(Q2Yhatcum)));
    
    
    
end

function p = tTest(r,n)

p = betainc(1-r.^2,n/2-1,1/2);

function modelPred=O2plsPred(X,Y,model,oax,oay,dir)

if(strcmpi(dir,'x'))
    To=[];
    %if(oax>0)
    
    
    for(i = 1:oax)
        to=X* model.Wo(:,i) *inv(model.Wo(:,i)'*model.Wo(:,i));
        To=[To,to];
        X=X-to*model.Pyo(:,i)';
    end
    
    T=X*model.W*inv(model.W'*model.W);
    Yhat=T*model.Bts{oax+1,oay+1}*model.C';
    modelPred.T=T;
    modelPred.Yhat=Yhat;
    modelPred.To=To;
    return;
    %end
end

if(strcmpi(dir,'y'))
    Uo=[];
    %	if(oay>0)
    for(i = 1:oay)
        uo=Y*model.Co(:,i)*inv(model.Co(:,i)'*model.Co(:,i));
        Uo=[Uo,uo];
        Y=Y-uo*model.Pxo(:,i)';
    end
    
    U=Y*model.C*inv(model.C'*model.C);
    Xhat=U*model.Bus{oax+1,oay+1}*model.W';
    modelPred.U=U;
    modelPred.Xhat=Xhat;
    modelPred.Uo=Uo;
    return;
    
    %end
    
end


function model=O2pls(X,Y,pax,oax,oay,mtype,varargin)
%------------------------------------------------
%------------------------------------------------

if(~isempty(varargin)) %double checks here for backward compatibilty
    if(~isempty(varargin{1}))
        splitSize=varargin{1};
    else
        splitSize=300; %this is the Y block size (if mtype='large');
    end
end

ssx=sum(sum(X.^2));
ssy=sum(sum(Y.^2));

if(strcmp(mtype,'large')) %large matrices
    [m,n]=size(Y);
    nsplits=floor(n/splitSize);
    Wtot=[];
    Ctot=[];
    
    h=waitbar(0,['estimating Wp for Y-blocks...']);%,num2str(1),' of ',num2str(nsplits)]);
    for i = 1:nsplits
        %disp(['ysplit',num2str(i)]);
        waitbar((i-1)/nsplits,h);%,['estimating Wp for Y-block: ',num2str(i),' of ',num2str(nsplits)]);
        [Wtmp,Stmp,Ctmp] = svd([Y(:, [((i-1)*splitSize+1):(i*splitSize)] )'*X]',0);
        Wtmp=Wtmp(:,1:pax);
        Stmp=Stmp(1:pax,1:pax);
        Wtot=[Wtot,Wtmp*Stmp];
        
        if(i==nsplits && mod(n,splitSize)~=0)
            %disp(['last Y-block - got remaining part of ',num2str(mod(n,splitSize)),' vars...']);
            [Wtmp,Stmp,Ctmp] = svd([Y(:, [((i*splitSize)+1):end] )'*X]',0);
            if(mod(n,splitSize)<pax)
                paxTmp=mod(n,splitSize);
            else
                paxTmp=pax;
            end
            Wtmp=Wtmp(:,1:paxTmp);
            Stmp=Stmp(1:paxTmp,1:paxTmp);
            Wtot=[Wtot,Wtmp*Stmp];
            
        end
    end
    
    close(h);
    
    
    
    [W,S,c]=svd(Wtot,0);
    clear Wtot;
    W=W(:,1:pax);
    S=S(1:pax,1:pax);
    Ts{1}=X*W;
    T=Ts{1};
    C=Y'*T;
    C=[C*diag(1./sqrt(sum(C(:,1:pax).^2)))];
    
else %conventional
    if ~isempty(find(isinf(Y'*X))) || ~isempty(find(isnan(Y'*X)))
        model = [];
        return;
    end
    [W,S,C] = svd([Y'*X]',0);
    W=W(:,1:pax);
    C=C(:,1:pax);
    S=S(1:pax,1:pax);
    
    %     [~,~,~,~,W,~,C,~] = plsr(X,Y,pax);
    Ts{1}=X*W;
    T=Ts{1};
end


Exy=X-T*W';
Xr=X;

R2Xo=0;
R2Xcorr=sum(sum((Ts{1}*W').^2))/ssx;
R2X=1-sum(sum(Exy.^2))/ssx;

Wo=[];
Pyo=[];
To=[];

for(i = 1:oax)
    [wo,syo,wor]=svd([Exy'*Ts{i}],0);
    wo=wo(:,1); %/sqrt(wo'*wo);
    to=Xr*wo; %/(wo'*wo);
    pyo=Xr'*to/(to'*to);
    Xr=Xr-to*pyo';
    
    Wo=[Wo,wo];
    Pyo=[Pyo,pyo];
    To=[To,to];
    
    Ts{i+1}=Xr*W;
    Exy=X-Ts{i+1}*W'-To*Pyo';
    
    R2Xo=[R2Xo,sum(sum((To*Pyo').^2))/ssx];
    R2Xcorr=[R2Xcorr,sum(sum((Ts{i+1}*W').^2))/ssx];
    R2X=[R2X,1-sum(sum(Exy.^2))/ssx];
    T=Ts{i+1};
end



Yr=Y;
Us{1}=Yr*C;
U=Us{1};
Fxy=Y-Us{1}*C';

R2Yo=0;
R2Ycorr=sum(sum((Us{1}*C').^2))/ssy;
R2Y=1-sum(sum(Fxy.^2))/ssy;

Uo=[];
Pxo=[];
Co=[];

for(i = 1:oay)
    [co,sxo,cor]=svd(Fxy'*Us{i},0);
    co=co(:,1); %/sqrt(co'*co);
    uo=Yr*co; %/(co'*co);
    pxo=Yr'*uo/(uo'*uo);
    Yr=Yr-uo*pxo';
    
    Co=[Co,co];
    Pxo=[Pxo,pxo];
    Uo=[Uo,uo];
    
    Us{i+1}=Yr*C;
    Fxy=Y-Us{i+1}*C'-Uo*Pxo';
    
    R2Yo=[R2Yo,sum(sum((Uo*Pxo').^2))/ssy];
    R2Ycorr=[R2Ycorr,sum(sum((Us{i+1}*C').^2))/ssy];
    R2Y=[R2Y,1-sum(sum(Fxy.^2))/ssy];
    U=Us{i+1};
    
end

Bus=cell(1,1);
bts=cell(1,1);
for(i = 1:oax+1)
    for(j = 1:oay+1)
        Bus{i,j}= inv(Us{j}'*Us{j})*Us{j}'*Ts{i};
        Bts{i,j}= inv(Ts{i}'*Ts{i})*Ts{i}'*Us{j};
    end
end

Bu= inv(U'*U)*U'*T;
Bt= inv(T'*T)*T'*U;






R2Yhat=[];
R2Xhat=[];


%This is OC style - keeping for now... :
for(i = 1:oax+1)
    BtTmp=inv(Ts{i}'*Ts{i})*Ts{i}'*U;
    YhatTmp=Ts{i}* BtTmp * C';
    R2Yhat=[R2Yhat, 1-sum(sum((YhatTmp-Y).^2))/ssy];
end
for(i = 1:oay+1)
    BuTmp=inv(Us{i}'*Us{i})*Us{i}'*T;
    XhatTmp=Us{i}* BuTmp * W';
    R2Xhat=[R2Xhat, 1-sum(sum((XhatTmp-X).^2))/ssx];
    
end


model.T=T;
model.Ts=Ts;
model.W=W;
model.Wo=Wo;
model.Pyo=Pyo;
model.To=To;

model.U=U;
model.Us=Us;
model.C=C;
model.Co=Co;
model.Pxo=Pxo;
model.Uo=Uo;

model.Bt=Bt;
model.Bu=Bu;


model.Bts=Bts;
model.Bus=Bus;

model.R2X=R2X;
model.R2Xcorr=R2Xcorr;
model.R2Xo=R2Xo;
model.R2Xhat=R2Xhat;

model.R2Y=R2Y;
model.R2Ycorr=R2Ycorr;
model.R2Yo=R2Yo;
model.R2Yhat=R2Yhat;

model.ssx=ssx;
model.ssy=ssy;

function [SIGMA, grouped_genes_bin, grouped_genes_idx_bin] = Binarization(data, groupIdx, biomarkersIdx)

% Binarization by biclustering
% 'data'                  - the input data matrix
% 'groupIdx'              - the cell categorization results obtained by HYSC
% 'biomarkersIdx'         - the variable genes obtained by HYSC
% 'SIGMA'                 - the binary output
% 'grouped_genes_bin'     - grouped gene symbols
% 'grouped_genes_idx_bin' - grouped gene index

%% association genes with cell clusters
grouped_genes_bin = [];
X_binary = zeros(length(groupIdx), length(biomarkersIdx));
for i = 1:length(biomarkersIdx)
    temp = data(:,biomarkersIdx(i));
    
    % % binarization by k-means
    IDX = kmeans(temp, 2, 'replicates', 20);
    mu(1) = mean(temp(IDX == 1));
    mu(2) = mean(temp(IDX == 2));
    [~, IX] = max(mu);
    thresh = min(temp(IDX == IX));
    X_binary(temp>=thresh,i) = 1;
    X_binary(temp<thresh,i) = 0;
end
clear P;
temp = [];
for i = floor(log2(max(groupIdx))):2^max(groupIdx) % categorize markers into numbers of groups which are <= max(groupIdx)
    if i == 1
        counts = zeros(2,max(groupIdx));
        for j = 1:max(groupIdx) % for each cell category
            counts(1,j) = length(find(X_binary(groupIdx==j,:)==1));
            counts(2,j) = length(find(X_binary(groupIdx==j,:)==0));
            grouped_genes{1,j} = X_binary(groupIdx==j,:)';
        end
        grouped_genes_idx{1} = biomarkersIdx;
        K = sum(counts);
        n = sum(counts,2);
        M = sum(K);
        for k = 1:size(counts,1)
            P{i}(k,:) = 1-cdf('hyge',counts(k,:), M, K, n(k));
        end
        temp = P{i}(1,:);
        temp(find(P{i}(1,:)<0.05)) = 1;
        temp(find(P{i}(1,:)>=0.05)) = 0;
    else
        IDX = kmeans(X_binary', i, 'replicates', 20, 'Distance', 'hamming');
        temp = [];
        geneNgroup = zeros(1,i);
        for jj = 1:i % for each marker subgroup
            counts = zeros(2,max(groupIdx));
            for j = 1:max(groupIdx) % for each cell category
                counts(1,j) = length(find(X_binary(groupIdx==j,IDX==jj)==1));
                counts(2,j) = length(find(X_binary(groupIdx==j,IDX==jj)==0));
                grouped_genes{jj,j} = X_binary(groupIdx==j,IDX==jj)';
            end
            grouped_genes_idx{jj} = biomarkersIdx(IDX==jj);
            K = sum(counts);
            n = sum(counts,2);
            M = sum(K);
            for k = 1:size(counts,1)
                P{i,jj}(k,:) = 1-cdf('hyge',counts(k,:), M, K, n(k));
            end
            temp_P = P{i,jj}(1,:);
            temp_P(find(P{i,jj}(1,:)<0.05)) = 1;
            temp_P(P{i,jj}(1,:)>=0.05) = 0;
            temp = [temp; temp_P];
            geneNgroup(jj) = length(find(IDX==jj));
        end
        D = pdist(temp, 'hamming');
        if min(D) == 0 || min(geneNgroup)<10
            if exist('temp_ori')
                SIGMA = temp_ori; % rows are gene clusters and columns are cells clusters
                grouped_genes_bin = grouped_genes_ori;
                grouped_genes_idx_bin = grouped_genes_idx_ori;
            else
                SIGMA = [];
                grouped_genes_bin = [];
                grouped_genes_idx_bin = [];
            end
            return;
            
        end
    end
    
    temp_ori = temp;
    grouped_genes_ori = grouped_genes;
    grouped_genes_idx_ori = grouped_genes_idx;
end


function output = HYSC_s(data, varargin)

% run single layer HYSC
% Input:
% data      - M by N matrix, M is the number of variables (genes) and N is the number of cells.
% varargin  -  'ClusterNum'            - the given number of clusters. If not given (ClusterNum = 0), the algorithm automatically estimates the number of clusters
%              'Nmax'                  - the maximum number of cell clusters. By default Nmax = 10
%              'minClusterSize_number' - the minimum number of samples (cells) in each cluster, by default, minClusterSize_number = 3
%              'dimClustering'         - the number of PCs adopted by k-means clustering
%              'nc'                    - OPLSDA parameters, number of correlated variables in X, by default 1
%              'ncox'                  - OPLSDA parameters, number of orthogonal variables in X, by default 1
%              'ncoy'                  - OPLSDA parameters, number of orthogonal variables in Y, by default 0
%              'r2cutoff'              - OPLSDA parameters, r2 cutoff value to identify discriminatory variables, by default 0.5
%              'p-cutoff'              - OPLSDA parameters, p cutoff value to identify discriminatory variables, by default 0.05
%              'preprocessing'         - OPLSDA parameters, pre processing methods, 'uv', 'pa', ''mc', by default 'uv'.
%              'perp'                  - perplexity of tSNE
%
% Output:
% the structure 'output'
% 'output.groupIdx'   - the cells categorization result based on the estimated number of cells clusters
% 'output.biomarkers' - the variable genes identified by hybrid clustering
% 'output.P'          - the p-values of variable genes
%
% Author Jie Hao & Xin Zou, SJTU, China
% Copyright Jie Hao & Xin Zou



varList = {'ClusterNum', 'Nmax', 'minClusterSize_number', 'dimClustering', 'nc', 'ncox', 'ncoy', 'r2cutoff', 'pcutoff', 'preprocessing', 'perp'};
convgThrh = {0, 10, 10, 5, 1, 1, 0, 0.5, 0.05, 'uv', 30}; % default values
if ~isempty(varargin)
    L = length(varargin);
    if L == 1
        varargin = varargin{1};
        L = length(varargin);
    end
    for i = 1:2:L-1
        idx = find(strcmp(varargin{i}, varList));
        convgThrh{idx} = varargin{i+1};
    end
end

ClusterNum = convgThrh{1}; % if being 0, the algorithm automatically determine the convergence, otherwise, ClusterNum is the number of cell clusters
Nmax = convgThrh{2}; % number of cell clusters by default
minClusterSize_number = convgThrh{3}; % 'minClusterSize_number'
dim = convgThrh{4}; % the number of PCs for k-means clustering
oplsda_para.nc = convgThrh{5};
oplsda_para.ncox = convgThrh{6};
oplsda_para.ncoy = convgThrh{7};
oplsda_para.r2cutoff = convgThrh{8};
oplsda_para.pcutoff = convgThrh{9};
oplsda_para.preprocessing = convgThrh{10};
perp = convgThrh{11}; % tSNE perplexity
oplsda_para.NP = 0; % by default no permutation test is performed

%% preprocessing
var = 1:size(data,1);
data = data';
data(data==0) = 10^(-6)*rand(1, length(find(data==0)));
%% HYSC
N = 2; % initial number of clusters
biomarkers_idx_pool_ori2 = []; % for termitation of the following 'while' loop
biomarkers_idx_pool_unsigned = 1:size(data,2);
randVar = [];
while 1 % loop for each given N find biomarkers
    idx_biomarkers = unique(biomarkers_idx_pool_unsigned);
    % % PCA & k-means
    X = data(:,[idx_biomarkers, randVar]);
    X_temp = (X-ones(size(X, 1),1)*mean(X,1));
    [~,SCORE] = pca(X_temp);
    temp = SCORE(:,1:min(dim,size(SCORE,2)));
    IDX = kmeans(temp,N, 'replicates', 20);
    
    % pairwise comparison between clusters
    biomarkers_idx_pool_unsigned = [];
    C = combnk(1:N,2);
    P_temp = ones(size(C,1)+N, size(data,2));
    randVar = [];
    for i = 1:size(C,1)
        tempData1 = data(IDX==C(i,1), :);
        tempData2 = data(IDX==C(i,2), :);
        model = OPLSDA(tempData1, tempData2, var, oplsda_para);
        [~, randTemp] = sort(model.W2);
        randTemp = randTemp(1:floor(length(model.sig_idx)*0.05)+1);
        if isempty(model) || isempty(model.sig_idx) || length(model.sig_idx)<2 || length(find(IDX==C(i,1)))<minClusterSize_number || length(find(IDX==C(i,2)))<minClusterSize_number
            output.groupIdx = ones(size(IDX));
            return;
        end
        
        biomarkers_idx_pool_unsigned = [biomarkers_idx_pool_unsigned, model.sig_idx];
        randVar = [randVar, randTemp];
        P_temp(i,:) = model.P;
    end
    
    
    % Hybrid clustering termination detection
    biomarkers_idx_pool_unsigned = unique(biomarkers_idx_pool_unsigned);
    randVar = unique(randVar);
    if ~isempty(biomarkers_idx_pool_ori2)
        c = intersect(biomarkers_idx_pool_ori2, biomarkers_idx_pool_unsigned);
        if (length(biomarkers_idx_pool_unsigned)-length(c))/length(biomarkers_idx_pool_ori2)<0.05 % if the variable transcripts identified in two consecutive iterations have more than 95% overlap, converge
            break;
        end
    end
    biomarkers_idx_pool_ori2 = biomarkers_idx_pool_unsigned;
end

if size(X,1)<5000
    X = data(:,biomarkers_idx_pool_unsigned);
    [COEFF,SCORE,latent] = pca(X);
    ydata = tsne(X,[],SCORE(:,1:min(dim, length(model.sig_idx))), perp);
    [~, p] = ttest2(ydata(IDX==1, 1), ydata(IDX==2, 1));
    if p>=0.05
        output.groupIdx = ones(size(IDX));
        return;
    end
end

P_min = min(P_temp);

biomarkers_idx_pool_ori = [];
biomarkers_idx_pool_global = unique(biomarkers_idx_pool_unsigned);
[SIGMA, grouped_genes_bin, grouped_genes_idx_bin] = Binarization(data, IDX, biomarkers_idx_pool_global);
% if the maximum number of clusters reached, the algorithm terminates
if N == min(Nmax, ClusterNum)
    groupIdx = IDX;
    output.groupIdx = groupIdx;
    output.biomarkers = biomarkers_idx_pool_global;
    output.P = P_min(biomarkers_idx_pool_global);
    output.SIGMA = SIGMA;
    return;
end

N = N+1;

%% Hybrid clustering iteration
while 1 % loop for choosing N
    if N>2
        IDX_lastIte = IDX;
        SIGMA_lastIte =SIGMA;
        biomarkers_idx_pool_ori = biomarkers_idx_pool_global;
        P_ori = P_min;
    end
    biomarkers_idx_pool_ori2 = []; % for termitation of the following 'while' loop
    
    while 1 % loop for each given N find biomarkers
        idx_biomarkers = unique(biomarkers_idx_pool_unsigned);
        % % PCA & k-means
        X = data(:,[idx_biomarkers, randVar]);
        X_temp = (X-ones(size(X, 1),1)*mean(X,1));
        [~,SCORE] = pca(X_temp);
        temp = SCORE(:,1:min(dim,size(SCORE,2)));
        IDX = kmeans(temp,N, 'replicates', 20);
        
        a = unique(IDX);
        out = [a,histc(IDX(:),a)]; % the number of samples in each cell cluster
        if min(out(:,2))<minClusterSize_number && ClusterNum ==0 % if the size of cluster is too small, algorithm converges.
            groupIdx = IDX_lastIte;
            biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
            output.groupIdx = groupIdx;
            output.biomarkers = biomarkers_idx_pool_global;
            output.P = P_ori(biomarkers_idx_pool_global);
            output.SIGMA = SIGMA_lastIte;
            return;
        end
        
        biomarkers_idx_pool_unsigned = [];
        % pairwise comparison between clusters
        C = combnk(1:N,2);
        P_temp = ones(size(C,1)+N, size(data,2));
        biomarkers_cell = cell(1, size(C,1));
        data_cluster_idx_cell = cell(size(C));
        randVar = [];
        for i = 1:size(C,1)
            tempData1 = data(IDX==C(i,1), :);
            tempData2 = data(IDX==C(i,2), :);
            model = OPLSDA(tempData1, tempData2, var, oplsda_para);
            [~, randTemp] = sort(model.W2);
            randTemp = randTemp(1:floor(length(model.sig_idx)*0.05)+1);
            if isempty(model)
                groupIdx = IDX_lastIte;
                biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
                output.groupIdx = groupIdx;
                output.biomarkers = biomarkers_idx_pool_global;
                output.P = P_ori(biomarkers_idx_pool_global);
                output.SIGMA = SIGMA_lastIte;
                return;
            end
            biomarkers_idx_pool_unsigned = [biomarkers_idx_pool_unsigned, model.sig_idx];
            randVar = [randVar, randTemp];
            P_temp(i,:) = model.P;
            biomarkers_cell{i} = model.sig_idx;
            data_cluster_idx_cell{i, 1} = find(IDX==C(i,1));
            data_cluster_idx_cell{i, 2} = find(IDX==C(i,2));
        end
        
        for i = 1:N
            tempData1 = data(IDX==i, :);
            tempData2 = data(IDX~=i, :);
            model = OPLSDA(tempData1, tempData2, var, oplsda_para);
            [~, randTemp] = sort(model.W2);
            randTemp = randTemp(1:floor(length(model.sig_idx)*0.05)+1);
            if isempty(model)
                groupIdx = IDX_lastIte;
                biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
                output.groupIdx = groupIdx;
                output.biomarkers = biomarkers_idx_pool_global;
                output.P = P_ori(biomarkers_idx_pool_global);
                output.SIGMA = SIGMA_lastIte;
                return;
            end
            P_temp(size(C,1)+1,:) = model.P;
            biomarkers_idx_pool_unsigned = [biomarkers_idx_pool_unsigned, model.sig_idx];
            randVar = [randVar, randTemp];
        end
        P_min = min(P_temp);
        
        
        % Hybrid clustering termination detection
        biomarkers_idx_pool_unsigned = unique(biomarkers_idx_pool_unsigned);
        if ~isempty(biomarkers_idx_pool_ori2)
            c = intersect(biomarkers_idx_pool_ori2, biomarkers_idx_pool_unsigned);
            if (length(biomarkers_idx_pool_unsigned)-length(c))/length(biomarkers_idx_pool_ori2)<0.1 % if the variable transcripts identified in two consecutive iterations have more than 90% overlap, converge
                break;
            end
        end
        biomarkers_idx_pool_ori2 = biomarkers_idx_pool_unsigned;
    end
    
    
    if ~isempty(find(cellfun('isempty',biomarkers_cell))) && ClusterNum ==0 % if no variable transcript observed in any comparison, converge
        groupIdx = IDX_lastIte;
        biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
        output.groupIdx = groupIdx;
        output.biomarkers = biomarkers_idx_pool_global;
        output.P = P_ori(biomarkers_idx_pool_global);
        output.SIGMA = SIGMA_lastIte;
        return;
    end
    
    % % proposed algorithm
    biomarkers_idx_pool_global = unique([biomarkers_idx_pool_unsigned, biomarkers_idx_pool_global]);
    [SIGMA, grouped_genes_bin, grouped_genes_idx_bin] = Binarization(data, IDX, biomarkers_idx_pool_global);
    if N == Nmax  && ClusterNum ==0 || N == ClusterNum % if the maximum number of clusters reached or N has reached the given number of clusters, converge
        groupIdx = IDX;
        biomarkers_idx_pool_global = unique(biomarkers_idx_pool_global);
        output.groupIdx = groupIdx;
        output.biomarkers = biomarkers_idx_pool_global;
        output.P = P_min(biomarkers_idx_pool_global);
        output.SIGMA = SIGMA;
        return;
    end
    
    
    if ~isempty(SIGMA)
        D = pdist(SIGMA', 'hamming');
        if min(D) == 0
            groupIdx = IDX_lastIte;
            biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
            output.groupIdx = groupIdx;
            output.biomarkers = biomarkers_idx_pool_global;
            output.P = P_min(biomarkers_idx_pool_global);
            output.SIGMA = SIGMA_lastIte;
            return;
        end
    else
        groupIdx = IDX_lastIte;
        biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
        output.groupIdx = groupIdx;
        output.biomarkers = biomarkers_idx_pool_global;
        output.P = P_min(biomarkers_idx_pool_global);
        output.SIGMA = SIGMA_lastIte;
        return;
    end
    
    
    N=N+1;
    
end

function ydata = tsne(X, labels, no_dims, initial_dims, perplexity)
%TSNE Performs symmetric t-SNE on dataset X
%
%   mappedX = tsne(X, labels, no_dims, initial_dims, perplexity)
%   mappedX = tsne(X, labels, initial_solution, perplexity)
%
% The function performs symmetric t-SNE on the NxD dataset X to reduce its
% dimensionality to no_dims dimensions (default = 2). The data is
% preprocessed using PCA, reducing the dimensionality to initial_dims
% dimensions (default = 30). Alternatively, an initial solution obtained
% from an other dimensionality reduction technique may be specified in
% initial_solution. The perplexity of the Gaussian kernel that is employed
% can be specified through perplexity (default = 30). The labels of the
% data are not used by t-SNE itself, however, they are used to color
% intermediate plots. Please provide an empty labels matrix [] if you
% don't want to plot results during the optimization.
% The low-dimensional data representation is returned in mappedX.
%
%
% (C) Laurens van der Maaten, 2010
% University of California, San Diego


if ~exist('labels', 'var')
    labels = [];
end
if ~exist('no_dims', 'var') || isempty(no_dims)
    no_dims = 2;
end
if ~exist('initial_dims', 'var') || isempty(initial_dims)
    initial_dims = min(30, size(X, 2));
end
if ~exist('perplexity', 'var') || isempty(perplexity)
    perplexity = 30;
end

% First check whether we already have an initial solution
if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
    perplexity = initial_dims;
else
    initial_solution = false;
end

%     sum_X = sum(X .^ 2, 2);
%     D = real((bsxfun(@plus, sum_X, bsxfun(@plus, sum_X', -2 * (X * X')))).^(1/2));
%     D(logical(eye(size(D))))=0;
%     Var = zeros(size(X,1),1);
%     for i = 1:length(D)
%         temp = D(i,:);
%         temp(temp==0) = [];
%         Var(i) = var(temp);
%     end
%     perplexity = log(det(diag(Var)))/2;


% Normalize input data
X = X - min(X(:));
X = X / max(X(:));
X = bsxfun(@minus, X, mean(X, 1));

% Perform preprocessing using PCA
if ~initial_solution
    disp('Preprocessing data using PCA...');
    if size(X, 2) < size(X, 1)
        C = X' * X;
    else
        C = (1 / size(X, 1)) * (X * X');
    end
    [M, lambda] = eig(C);
    [lambda, ind] = sort(diag(lambda), 'descend');
    M = M(:,ind(1:initial_dims));
    lambda = lambda(1:initial_dims);
    if ~(size(X, 2) < size(X, 1))
        M = bsxfun(@times, X' * M, (1 ./ sqrt(size(X, 1) .* lambda))');
    end
    X = bsxfun(@minus, X, mean(X, 1)) * M;
    clear M lambda ind
end

% Compute pairwise distance matrix
sum_X = sum(X .^ 2, 2);
D = real((bsxfun(@plus, sum_X, bsxfun(@plus, sum_X', -2 * (X * X')))));

% Compute joint probabilities
P = d2p(D, perplexity, 1e-5);                                           % compute affinities using fixed perplexity
clear D

% Run t-SNE
if initial_solution
    ydata = tsne_p(P, labels, ydata);
else
    ydata = tsne_p(P, labels, no_dims);
end

function ydata = tsne_d(D, labels, no_dims, perplexity)
%TSNE_D Performs symmetric t-SNE on the pairwise Euclidean distance matrix D
%
%   mappedX = tsne_d(D, labels, no_dims, perplexity)
%   mappedX = tsne_d(D, labels, initial_solution, perplexity)
%
% The function performs symmetric t-SNE on the NxN pairwise Euclidean
% distance matrix D to construct an embedding with no_dims dimensions
% (default = 2). An initial solution obtained from an other dimensionality
% reduction technique may be specified in initial_solution.
% The perplexity of the Gaussian kernel that is employed can be specified
% through perplexity (default = 30). The labels of the data are not used
% by t-SNE itself, however, they are used to color intermediate plots.
% Please provide an empty labels matrix [] if you don't want to plot
% results during the optimization.
% The data embedding is returned in mappedX.
%
%
% (C) Laurens van der Maaten, 2010
% University of California, San Diego


if ~exist('labels', 'var')
    labels = [];
end
if ~exist('no_dims', 'var') || isempty(no_dims)
    no_dims = 2;
end
if ~exist('perplexity', 'var') || isempty(perplexity)
    perplexity = 30;
end

% First check whether we already have an initial solution
if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
else
    initial_solution = false;
end

% Compute joint probabilities
D = D / max(D(:));                                                      % normalize distances
P = d2p(D .^ 2, perplexity, 1e-5);                                      % compute affinities using fixed perplexity

% Run t-SNE
if initial_solution
    ydata = tsne_p(P, labels, ydata);
else
    ydata = tsne_p(P, labels, no_dims);
end

function ydata = tsne_p(P, labels, no_dims)
%TSNE_P Performs symmetric t-SNE on affinity matrix P
%
%   mappedX = tsne_p(P, labels, no_dims)
%
% The function performs symmetric t-SNE on pairwise similarity matrix P
% to create a low-dimensional map of no_dims dimensions (default = 2).
% The matrix P is assumed to be symmetric, sum up to 1, and have zeros
% on the diagonal.
% The labels of the data are not used by t-SNE itself, however, they
% are used to color intermediate plots. Please provide an empty labels
% matrix [] if you don't want to plot results during the optimization.
% The low-dimensional data representation is returned in mappedX.
%
%
% (C) Laurens van der Maaten, 2010
% University of California, San Diego


if ~exist('labels', 'var')
    labels = [];
end
if ~exist('no_dims', 'var') || isempty(no_dims)
    no_dims = 2;
end

% First check whether we already have an initial solution
if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
else
    initial_solution = false;
end

% Initialize some variables
n = size(P, 1);                                     % number of instances
momentum = 0.5;                                     % initial momentum
final_momentum = 0.8;                               % value to which momentum is changed
mom_switch_iter = 250;                              % iteration at which momentum is changed
stop_lying_iter = 100;                              % iteration at which lying about P-values is stopped
max_iter = 200;                                    % maximum number of iterations
epsilon = 500;                                      % initial learning rate
min_gain = .01;                                     % minimum gain for delta-bar-delta

% Make sure P-vals are set properly
P(1:n + 1:end) = 0;                                 % set diagonal to zero
P = 0.5 * (P + P');                                 % symmetrize P-values
P = max(P ./ sum(P(:)), realmin);                   % make sure P-values sum to one
const = sum(P(:) .* log(P(:)));                     % constant in KL divergence
if ~initial_solution
    P = P * 4;                                      % lie about the P-vals to find better local minima
end

% Initialize the solution
if ~initial_solution
    ydata = .0001 * randn(n, no_dims);
end
y_incs  = zeros(size(ydata));
gains = ones(size(ydata));

% Run the iterations
for iter=1:max_iter
    
    % Compute joint probability that point i and j are neighbors
    sum_ydata = sum(ydata .^ 2, 2);
    num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * (ydata * ydata')))); % Student-t distribution
    num(1:n+1:end) = 0;                                                 % set diagonal to zero
    Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
    
    % Compute the gradients (faster implementation)
    L = (P - Q) .* num;
    y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
    
    % Update the solution
    gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...         % note that the y_grads are actually -y_grads
        + (gains * .8) .* (sign(y_grads) == sign(y_incs));
    gains(gains < min_gain) = min_gain;
    y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
    ydata = ydata + y_incs;
    ydata = bsxfun(@minus, ydata, mean(ydata, 1));
    
    % Update the momentum if necessary
    if iter == mom_switch_iter
        momentum = final_momentum;
    end
    if iter == stop_lying_iter && ~initial_solution
        P = P ./ 4;
    end
    
    % Print out progress
    if ~rem(iter, 10)
        cost = const - sum(P(:) .* log(Q(:)));
        %             disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
    end
    
    % Display scatter plot (maximally first three dimensions)
    if ~rem(iter, 10) && ~isempty(labels)
        if no_dims == 1
            scatter(ydata, ydata, 9, labels, 'filled');
        elseif no_dims == 2
            scatter(ydata(:,1), ydata(:,2), 9, labels, 'filled');
        else
            scatter3(ydata(:,1), ydata(:,2), ydata(:,3), 40, labels, 'filled');
        end
        axis tight
        axis off
        drawnow
    end
end
