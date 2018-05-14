function varargout = AdvancedSetting(varargin)
% ADVANCEDSETTING MATLAB code for AdvancedSetting.fig
%      ADVANCEDSETTING, by itself, creates a new ADVANCEDSETTING or raises the existing
%      singleton*.
%
%      H = ADVANCEDSETTING returns the handle to a new ADVANCEDSETTING or the handle to
%      the existing singleton*.
%
%      ADVANCEDSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCEDSETTING.M with the given input arguments.
%
%      ADVANCEDSETTING('Property','Value',...) creates a new ADVANCEDSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AdvancedSetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AdvancedSetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AdvancedSetting

% Last Modified by GUIDE v2.5 07-May-2018 15:32:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AdvancedSetting_OpeningFcn, ...
                   'gui_OutputFcn',  @AdvancedSetting_OutputFcn, ...
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


% --- Executes just before AdvancedSetting is made visible.
function AdvancedSetting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AdvancedSetting (see VARARGIN)

% Choose default command line output for AdvancedSetting
handles.output = hObject;
set(handles.minClusterSize_number,'Value',10);
set(handles.dimClustering,'Value',5);
set(handles.r2cutoff,'Value',0.5);
set(handles.pcutoff,'Value',0.05);
set(handles.Perplexity,'Value',30);
set(handles.cores,'Value',8);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AdvancedSetting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AdvancedSetting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function minClusterSize_number_Callback(hObject, eventdata, handles)
% hObject    handle to minClusterSize_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minClusterSize_number as text
%        str2double(get(hObject,'String')) returns contents of minClusterSize_number as a double
handles.minClusterSize_number.Value = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function minClusterSize_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minClusterSize_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dimClustering_Callback(hObject, eventdata, handles)
% hObject    handle to dimClustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dimClustering as text
%        str2double(get(hObject,'String')) returns contents of dimClustering as a double
handles.dimClustering.Value = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dimClustering_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dimClustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r2cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to r2cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r2cutoff as text
%        str2double(get(hObject,'String')) returns contents of r2cutoff as a double
handles.r2cutoff.Value = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function r2cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r2cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pcutoff_Callback(hObject, eventdata, handles)
% hObject    handle to pcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcutoff as text
%        str2double(get(hObject,'String')) returns contents of pcutoff as a double
handles.pcutoff.Value =  str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function pcutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Perplexity_Callback(hObject, eventdata, handles)
% hObject    handle to Perplexity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Perplexity as text
%        str2double(get(hObject,'String')) returns contents of Perplexity as a double
handles.Perplexity.Value =  str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Perplexity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Perplexity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cores_Callback(hObject, eventdata, handles)
% hObject    handle to cores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cores as text
%        str2double(get(hObject,'String')) returns contents of cores as a double
handles.cores.Value =  str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function cores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
handles.apply = 1;
guidata(hObject, handles);
set(gcf,'visible','off');

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cancel = 1;
guidata(hObject,handles);
set(gcf,'visible','off');



function filteringNo_Callback(hObject, eventdata, handles)
% hObject    handle to filteringNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filteringNo as text
%        str2double(get(hObject,'String')) returns contents of filteringNo as a double
handles.filteringNo.Value =  str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function filteringNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filteringNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filteringVar_Callback(hObject, eventdata, handles)
% hObject    handle to filteringVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filteringVar as text
%        str2double(get(hObject,'String')) returns contents of filteringVar as a double
handles.filteringVar.Value =  str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function filteringVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filteringVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
