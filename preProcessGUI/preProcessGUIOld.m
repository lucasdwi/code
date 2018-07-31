function varargout = preProcessGUI(varargin)
% PREPROCESS MATLAB code for preProcess.fig
%      PREPROCESS, by itself, creates a new PREPROCCESS or raises the existing
%      singleton*.
%
%      H = PREPROCESS returns the handle to a new PREPROCCESS or the handle to
%      the existing singleton*.
%
%      PREPROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCCESS.M with the given input arguments.
%
%      PREPROCESS('Property','Value',...) creates a new PREPROCCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preProccess_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preProccess_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preProccess

% Last Modified by GUIDE v2.5 13-Jul-2018 13:59:00

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

% --- Executes just before preProcess is made visible.
function preProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preProcess (see VARARGIN)
global eoi eventTs
% Choose default command line output for preProcess
handles.output = hObject;
% Prep zoom handles to pair with slider
handles.hZoom = zoom;
set(handles.hZoom,'ActionPostCallback',@zoomPostCallBack);
% Set up eoi
eoi = {'None'};
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes preProcess wait for user response (see UIRESUME)
% uiwait(handles.gui1);

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
global LFPTs out prop threshValue
% Get current channel
currentChan = get(handles.channels,'Value');

if isfield(handles,'in')
    in = handles.in;
    cla reset
    hold on
    % Plot positive and negative
    h1 = plot(handles.chanView,[LFPTs.tvec(1) LFPTs.tvec(end)],[threshValue threshValue],'--k');
    h2 = plot(handles.chanView,[LFPTs.tvec(1) LFPTs.tvec(end)],[-1*threshValue -1*threshValue],'--k');
    % Plot
    plot(LFPTs.tvec,in(currentChan,:),'-b')
    plot(LFPTs.tvec,out(currentChan,:),'-r');
    plot(LFPTs.tvec,prop(currentChan,:),'color',[.75 0 .75]);
else
    cla reset
    % Plot
    plot(handles.chanView,LFPTs.tvec,LFPTs.data(currentChan,:));
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
global LFPTs adfreq out threshValue prop
handles.eventTs
set(gcf,'pointer','watch')
drawnow
cla reset
% Grab on- and offset values and convert to samples
onsetValue = str2num(get(handles.onset,'String')).*adfreq;
offsetValue = str2num(get(handles.offset,'String')).*adfreq;
% Grab value in threshold
threshValue = str2num(get(handles.thresh,'String'));
% Get logical for threshold crosses
threshInds = abs(LFPTs.data)>=threshValue;
% Set threshold crosses to NaN
in = LFPTs.data;
in(threshInds) = NaN;

% Find values beyond threshold with on- and offset
for ii = 1:size(LFPTs.data,1)
    outInds = logicFind(1,threshInds(ii),'==');
    for jj = 1:length(outInds)
        if outInds(jj)-onsetValue < 0
            in(ii,1:outInds(jj)) = NaN;
        else
            in(ii,(outInds(jj)-onsetValue):outInds(jj)) = NaN;
        end
        if outInds(jj)+offsetValue > length(LFPTs.tvec)
            in(ii,outInds(jj):end) = NaN;
        else
            in(ii,outInds(jj):(outInds(jj)+offsetValue)) = NaN;
        end
    end
end
handles.in = in;
guidata(hObject,handles)
out = LFPTs.data;
out(~isnan(in)) = NaN;
% Propagate NaNs across channels
nanData = in;
nanData(:,logical(sum(isnan(in)))) = NaN;
% Determine which indices were propagated
propInd = nanData~=in;
propInd(threshInds) = 0;
% NaN values that were NOT propagated
prop = LFPTs.data;
prop(~propInd) = NaN;
prop(~isnan(out)) = NaN;
% Apply NaNs to in data
in = nanData;
% Get current channl
currentChan = get(handles.channels,'Value');
hold on
% Plot positive and negative
h1 = plot(handles.chanView,[LFPTs.tvec(1) LFPTs.tvec(end)],[threshValue threshValue],'--k');
h2 = plot(handles.chanView,[LFPTs.tvec(1) LFPTs.tvec(end)],[-1*threshValue -1*threshValue],'--k');
% Plot
plot(LFPTs.tvec,in(currentChan,:),'-b')
plot(LFPTs.tvec,out(currentChan,:),'-r')
plot(LFPTs.tvec,prop(currentChan,:),'color',[.75 0 .75])
% Set plot axes
xlabel(handles.chanView,'Time')
ylabel(handles.chanView,'V')
set(gcf,'pointer','arrow')

function channelName_Callback(hObject, eventdata, handles)
% hObject    handle to channelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channelName as text
%        str2double(get(hObject,'String')) returns contents of channelName as a double


% --- Executes during object creation, after setting all properties.
function channelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deleteChan.
function deleteChan_Callback(hObject, eventdata, handles)
% hObject    handle to deleteChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LFPTs
currentChan = get(handles.channels,'Value');
opts.Interpreter = 'tex';
opts.Default = 'No';
resp = questdlg(['Warning: You are about to delete all data in ',...
    LFPTs.label{currentChan},newline,...
    'Are you sure you want to do this?'],'Warning','Yes','No',opts);
switch resp
    case 'Yes'
        LFPTs.data(currentChan,:) = [];
        LFPTs.label{currentChan} = [];
        set(handles.channels,'String',LFPTs.label)
    case 'No'
end
% --- Executes on button press in saveFile.
function saveFile_Callback(hObject, eventdata, handles)
% hObject    handle to saveFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadFile.
function loadFile_Callback(hObject, eventdata, handles)
% hObject    handle to loadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in editName.
function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LFPTs
% Prompt user to input new name
newName = inputdlg('Enter new channel name','Edit Channel Name');
% If canceled, then newName will be empty
if ~isempty(newName)
    % Check if newName is already being used
    if ~isempty(find(contains(LFPTs.label,newName)))
        warndlg([newName,' is already being used as a name for a channel'])
    else
    % Get current channel number
    currentChan = get(handles.channels,'Value');
    % Replace channel name
    LFPTs.label(currentChan) = newName;
    % Reset popup menu
    set(handles.channels,'String',LFPTs.label)
    end
end

% --------------------------------------------------------------------
function open_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LFPTs currentPlot adfreq 

% Load file
uiopen('*.mat')
handles.eventTs = eventTs;
guidata(hObject,handles);
run('filtDwn.m')
% Change arrow to watch to indicate busy
set(handles.gui1,'pointer','watch')
% Populate popup menu with channel labels
set(handles.channels,'String',LFPTs.label)
% Get current channl
currentChan = get(handles.channels,'Value');
% Plot first channel
currentPlot = plot(handles.chanView,LFPTs.tvec,LFPTs.data(currentChan,:));
% Set plot axes
xlabel(handles.chanView,'Time')
ylabel(handles.chanView,'V')
% Set slider
set(handles.slider1,'Min',min(LFPTs.tvec),'Max',max(LFPTs.tvec),'Value',1);
% Change back to arrow
set(handles.gui1,'pointer','arrow')



% --------------------------------------------------------------------
function save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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

function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

adjustAxisLimits(handles)

function adjustAxisLimits(handles)
minVal = get(handles.slider1,'Value');
maxVal = minVal + diff(xlim);
set(handles.chanView,'XLim',[minVal maxVal])

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
