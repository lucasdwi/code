function varargout = preProcessGUI(varargin)
% PREPROCESSGUI MATLAB code for preProcessGUI.fig
%      PREPROCESSGUI, by itself, creates a new PREPROCESSGUI or raises the existing
%      singleton*.
%
%      H = PREPROCESSGUI returns the handle to a new PREPROCESSGUI or the handle to
%      the existing singleton*.
%
%      PREPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSGUI.M with the given input arguments.
%
%      PREPROCESSGUI('Property','Value',...) creates a new PREPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preProcessGUI

% Last Modified by GUIDE v2.5 16-Jul-2018 11:12:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @preProcessGUI_OutputFcn, ...
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


% --- Executes just before preProcessGUI is made visible.
function preProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preProcessGUI (see VARARGIN)

% Choose default command line output for preProcess
handles.output = hObject;
% Prep zoom handles to pair with slider
handles.hZoom = zoom;
set(handles.hZoom,'ActionPostCallback',@zoomPostCallBack);
handles.eoi = cell(1,1);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes preProcessGUI wait for user response (see UIRESUME)
% uiwait(handles.preProcGUI);

function zoomPostCallBack(hObject,event)
handles = guidata(hObject);
set(handles.slider1,'Value',sum(xlim)/2);

% --- Outputs from this function are returned to the command line.
function varargout = preProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in channels.
function channels_Callback(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channels
% Grab data
preProc = guidata(findobj('Tag','preProcGUI'));
% Get current channel
currentChan = get(handles.channels,'Value');
if isfield(handles,'in')
    threshValue = handles.threshValue;
    in = handles.in;
    cla reset
    plotData
else
    cla reset
    % Plot
    plot(handles.chanView,preProc.LFPTs.tvec,preProc.LFPTs.data(currentChan,:));
end
% Set plot axes
xlabel(handles.chanView,'Time')
ylabel(handles.chanView,'V')


% --- Executes during object creation, after setting all properties.
function channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in applyThresh.
function applyThresh_Callback(hObject, eventdata, handles)
% hObject    handle to applyThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
preProc = guidata(findobj('Tag','preProcGUI'));
set(gcf,'pointer','watch')
drawnow
cla reset
% Grab on- and offset values and convert to samples
onsetValue = str2num(get(handles.onset,'String')).*preProc.adfreq;
offsetValue = str2num(get(handles.offset,'String')).*preProc.adfreq;
handles.onsetValue = onsetValue;
handles.offsetValue = offsetValue;
% Grab value in threshold
threshValue = str2num(get(handles.thresh,'String'));
handles.threshValue = threshValue;
% Get logical for threshold crosses
threshInds = abs(preProc.LFPTs.data)>=threshValue;
% Set threshold crosses to NaN
in = preProc.LFPTs.data;
in(threshInds) = NaN;

% Find values beyond threshold with on- and offset
for ii = 1:size(preProc.LFPTs.data,1)
    outInds = logicFind(1,threshInds(ii),'==');
    for jj = 1:length(outInds)
        if outInds(jj)-onsetValue < 0
            in(ii,1:outInds(jj)) = NaN;
        else
            in(ii,(outInds(jj)-onsetValue):outInds(jj)) = NaN;
        end
        if outInds(jj)+offsetValue > length(preProc.LFPTs.tvec)
            in(ii,outInds(jj):end) = NaN;
        else
            in(ii,outInds(jj):(outInds(jj)+offsetValue)) = NaN;
        end
    end
end
handles.in = in;
guidata(hObject,handles)
handles.out = preProc.LFPTs.data;
handles.out(~isnan(in)) = NaN;
% Propagate NaNs across channels
nanData = in;
nanData(:,logical(sum(isnan(in)))) = NaN;
% Determine which indices were propagated
propInd = nanData~=in;
propInd(threshInds) = 0;
% NaN values that were NOT propagated
handles.prop = preProc.LFPTs.data;
handles.prop(~propInd) = NaN;
handles.prop(~isnan(handles.out)) = NaN;
% Apply NaNs to in data
in = nanData;
% Plot 
plotData
set(gcf,'pointer','arrow')
% Update handles
guidata(hObject, handles);

% --- Executes on button press in deleteChan.
function deleteChan_Callback(hObject, eventdata, handles)
% hObject    handle to deleteChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Grab data
preProc = guidata(findobj('Tag','preProcGUI'));

currentChan = get(handles.channels,'Value');
opts.Interpreter = 'tex';
opts.Default = 'No';
resp = questdlg(['Warning: You are about to delete all data in ',...
    preProc.LFPTs.label{currentChan},newline,...
    'Are you sure you want to do this?'],'Warning','Yes','No',opts);
switch resp
    case 'Yes'
        preProc.LFPTs.data(currentChan,:) = [];
        preProc.LFPTs.label{currentChan} = [];
        set(handles.channels,'String',preProc.LFPTs.label)
    case 'No'
end

% --- Executes on button press in editName.
function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Grab data
preProc = guidata(findobj('Tag','preProcGUI'));
% Prompt user to input new name
newName = inputdlg('Enter new channel name','Edit Channel Name');
% If canceled, then newName will be empty
if ~isempty(newName)
    % Check if newName is already being used
    if ~isempty(find(contains(preProc.LFPTs.label,newName)))
        warndlg([newName,' is already being used as a name for a channel'])
    else
    % Get current channel number
    currentChan = get(handles.channels,'Value');
    % Replace channel name
    preProc.LFPTs.label(currentChan) = newName;
    % Reset popup menu
    set(handles.channels,'String',preProc.LFPTs.label)
    end
end


function onset_Callback(hObject, eventdata, handles)
% hObject    handle to onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of onset as text
%        str2double(get(hObject,'String')) returns contents of onset as a double


% --- Executes during object creation, after setting all properties.
function onset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in runTrializeGUI.
function runTrializeGUI_Callback(hObject, eventdata, handles)
% hObject    handle to runTrializeGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run('trializeGUI.m')

% --- Executes on button press in runTrialize.
function runTrialize_Callback(hObject, eventdata, handles)
% hObject    handle to runTrialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Plots data
function plotData(hObject,eventdata,handles)

preProc = guidata(findobj('Tag','preProcGUI'));
currentChan = preProc.channels.Value;
hold on
% Plot positive and negative
plot(preProc.chanView,[preProc.LFPTs.tvec(1) preProc.LFPTs.tvec(end)],[preProc.threshValue preProc.threshValue],'--k');
plot(preProc.chanView,[preProc.LFPTs.tvec(1) preProc.LFPTs.tvec(end)],[-1*preProc.threshValue -1*preProc.threshValue],'--k');
% Plot
plot(preProc.LFPTs.tvec,preProc.in(currentChan,:),'-b')
plot(preProc.LFPTs.tvec,preProc.out(currentChan,:),'-r')
plot(preProc.LFPTs.tvec,preProc.prop(currentChan,:),'color',[.75 0 .75])
% Set plot axes
xlabel(preProc.chanView,'Time')
ylabel(preProc.chanView,'V')

% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load file
uiopen('*.mat')
handles.eventTs = eventTs;
handles.LFPTs = LFPTs;
handles.adfreq = adfreq;
guidata(hObject,handles);
run('filtDwn.m')
% Change arrow to watch to indicate busy
set(handles.preProcGUI,'pointer','watch')
drawnow
% Populate popup menu with channel labels
set(handles.channels,'String',LFPTs.label)
% Get current channl
currentChan = get(handles.channels,'Value');
% Plot first channel
plot(handles.chanView,LFPTs.tvec,LFPTs.data(currentChan,:));
% Set plot axes
xlabel(handles.chanView,'Time')
ylabel(handles.chanView,'V')
% Set slider
set(handles.slider1,'Min',min(LFPTs.tvec),'Max',max(LFPTs.tvec),'Value',1);
% Change back to arrow
set(handles.preProcGUI,'pointer','arrow')
