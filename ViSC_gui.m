function varargout = ViSC_gui(varargin)
% VISC_GUI MATLAB code for ViSC_gui.fig
%      VISC_GUI, by itself, creates a new VISC_GUI or raises the existing
%      singleton*.
%
%      H = VISC_GUI returns the handle to a new VISC_GUI or the handle to
%      the existing singleton*.
%
%      VISC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISC_GUI.M with the given input arguments.
%
%      VISC_GUI('Property','Value',...) creates a new VISC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViSC_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViSC_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViSC_gui

% Last Modified by GUIDE v2.5 11-May-2018 11:22:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViSC_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @ViSC_gui_OutputFcn, ...
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


% --- Executes just before ViSC_gui is made visible.
function ViSC_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViSC_gui (see VARARGIN)

% Choose default command line output for ViSC_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ViSC_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ViSC_gui_OutputFcn(hObject, eventdata, handles) 
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
     end
     close(gcf);
end



function HYSCoutput_Callback(hObject, eventdata, handles)
% hObject    handle to HYSCoutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HYSCoutput as text
%        str2double(get(hObject,'String')) returns contents of HYSCoutput as a double
handles.HYSCoutput = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function HYSCoutput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HYSCoutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
set(gcf,'visible','off');
quit;

function geneSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to geneSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of geneSymbol as text
%        str2double(get(hObject,'String')) returns contents of geneSymbol as a double
handles.geneSymbol = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function geneSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geneSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
disp('ViSC Running...');
[~, ~, raw] = xlsread(handles.HYSCoutput.String, 2); % load tSNE scores
tSNE = cell2mat(raw);
[~, ~, raw] = xlsread(handles.HYSCoutput.String, 3); % load data
geneList = raw(:,1);
if ~iscellstr(geneList)
    A = cellfun(@ischar,geneList);
    idx = find(A==0);
    for i = 1:length(idx)
        geneList{idx(i)} = num2str(geneList{idx(i)});
    end
end
data = cell2mat(raw(:,2:end));
ViSC(tSNE(:,1:2), data, geneList, handles.geneSymbol, 0);

% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex]  = uigetfile('*.xlsx', 'Select Data file');
set(handles.HYSCoutput,'String',fullfile(PathName, FileName));
% handles.dataPathStr = fullfile(PathName, FileName);
guidata(hObject, handles);


function h = ViSC(ydata, X, geneList, geneID, colorOpt)

% colorOpt - choose color, 1 - blue, 2 - red, 0 - default

Step_size = 64;
[ia] = find(ismember(geneList, geneID));
temp = X(ia, :); %
temp = temp-min(temp);
temp = temp/max(temp);
Step=1/Step_size;
CC = [];
for i=1:Step_size-1
    CC{i}=find(temp<=(1)-(i-1)*Step);
end
CCSp = {NaN(size(ydata,1),size(ydata,2))};
CCSp_cell = repmat(CCSp, 1, Step_size-1);
for i=1:Step_size-1
    CCSp_cell{i}(CC{i},:)=ydata(CC{i},:);
end
h=figure;
switch colorOpt
    case 0 %% defualt color
        colormap(jet);
        % colormap(cool);
        Color=get(gcf,'Colormap');
    case 1 %% blue color
        Color = ones(64,3);
        Color(1,:) = [204, 204, 204]/255;
        stepWide = 204/Step_size;
        for i = 2:Step_size
            Color(i,:) = [204-(i-1)*stepWide, 204-(i-1)*stepWide, 204]/255;
        end
    case 2 %% red color
        Color = ones(64,3);
        Color(1,:) = [204, 204, 204]/255;
        stepWide = 204/Step_size;
        for i = 2:Step_size
            Color(i,:) = [204, 204-(i-1)*stepWide, 204-(i-1)*stepWide]/255;
        end
    case 3 %% magenta color
        Color = ones(64,3);
        Color(1,:) = [204, 204, 204]/255;
        stepWide = 204/Step_size;
        for i = 2:Step_size
            Color(i,:) = [204, 204-(i-1)*stepWide, 204]/255;
        end
end
%%

colormap(Color);
h = subplot('Position',[0.13 0.16 0.7 0.75]);
plot(CCSp_cell{1}(:,1), CCSp_cell{1}(:,2), 'o', 'markersize', 6, 'MarkerFaceColor', Color(64,:),  'MarkerEdgeColor', Color(64,:));
% plot3(CCSp_cell{1}(:,1), CCSp_cell{1}(:,2), CCSp_cell{1}(:,3), 'o', 'markersize', 4, 'MarkerFaceColor', Color(64,:),  'MarkerEdgeColor', Color(64,:));
hold on
for i=2:Step_size-1
    plot(CCSp_cell{i}(:,1), CCSp_cell{i}(:,2), 'o', 'markersize', 6, 'MarkerFaceColor', Color(64-i,:),  'MarkerEdgeColor', Color(64-i,:));
%     plot3(CCSp_cell{i}(:,1), CCSp_cell{i}(:,2), CCSp_cell{i}(:,3), 'o', 'markersize', 4, 'MarkerFaceColor', Color(64-i,:),  'MarkerEdgeColor', Color(64-i,:));
end

xlabel('tSNE_1', 'fontsize',14);
ylabel('tSNE_2', 'fontsize',14);
title(geneID);
set(gca,'FontSize',14);
subplot('Position',[0.92 0.16 0.05 0.75]);
imagesc(1,1:-Step:0,(1:56)');
set(gca,'xTickLabel','');
set(gca, 'yTick',[0,1]);
set(gcf,'color','w');
set(gca,'FontSize',14);
