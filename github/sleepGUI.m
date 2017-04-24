function varargout = sleepGUI(varargin)
%SLEEPGUI M-file for sleepGUI.fig
%      SLEEPGUI, by itself, creates a new SLEEPGUI or raises the existing
%      singleton*.
%
%      H = SLEEPGUI returns the handle to a new SLEEPGUI or the handle to
%      the existing singleton*.
%
%      SLEEPGUI('Property','Value',...) creates a new SLEEPGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to sleepGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SLEEPGUI('CALLBACK') and SLEEPGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SLEEPGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sleepGUI

% Last Modified by GUIDE v2.5 02-Oct-2016 21:08:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sleepGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @sleepGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before sleepGUI is made visible.
function sleepGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for sleepGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.debugbox,'String',[{'1. Use browse and then load to open a file.'};{'2. Using the instantaneous amplitude plot find a channel with low signal/noise (i.e. sleep events are obvious)'};{'3. Adjust percentile and use "Test" to find a good threshold.'};{'4. With the cleanest channel selected, use "Set Cut-Off" to find intervals.'};{'5. Go through channels and "Plot" to review, if needed tune onset, offset, and percentile until all sleep events are selected across channels, always reselecting the cleanest channel before re-"Set Cut-Offs".'};{'6. Save!'}])

% UIWAIT makes sleepGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sleepGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushBrowse.
function pushBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filepath filename
% Clears plots
axes(handles.axes1)
cla
axes(handles.axes2)
cla
[pathname,filename] = uigetfile(('*.mat'),'File Selector');
filepath = [filename,pathname];
set(handles.fileShow,'String',filepath);

% --- Executes on button press in pushLoad.
function pushLoad_Callback(hObject, eventdata, handles)
% hObject    handle to pushLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filepath data
% Start busy indicator
set(gcf,'Pointer','watch')
% Load file
load(filepath);
% Add LFPTs and smInstAmp to global 'data'
data.LFPTs = LFPTs;
data.smInstAmp = smInstAmp;
% Find good (non-nan) indices
data.goodInds = find(~isnan(data.LFPTs.data(1,:)));
% Populate popup menu with channel names
set(handles.popupmenu1,'String',data.LFPTs.label)
% Finds default percentile (95%) of channel one
defaultP = prctile(data.smInstAmp(1,data.goodInds),95);
% Plot channel one instantaneous amplitude by default
axes(handles.axes2)
plot(data.LFPTs.tvec,data.smInstAmp(1,:),'r')
hline(defaultP,'--','95%')
title(['Instantaneous Amplitude of ',data.LFPTs.label{1}])
% Stop busy indicator
set(gcf,'Pointer','arrow')
% Grab eventTs
if exist('eventTs')
    data.eventTs = eventTs;
end

function editPercent_Callback(hObject, eventdata, handles)
% hObject    handle to editPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPercent as text
%        str2double(get(hObject,'String')) returns contents of editPercent as a double


% --- Executes during object creation, after setting all properties.
function editPercent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushPlot.
function pushPlot_Callback(hObject, eventdata, handles)
% hObject    handle to pushPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data
% Start busy indicator
set(gcf,'Pointer','watch')
% Get channel number from popupmenu1
chan = get(handles.popupmenu1,'Value');

% Use axes1 to plot raw signal
% Plot raw signal with intervals in axes1
axes(handles.axes1)
cla
plot(data.LFPTs.tvec,data.LFPTs.data(chan,:)); 
title(data.LFPTs.label{chan}); xlabel('Time (secs)'); ylabel('Amplitude (V)');
% Plot sleep intervals below
if ~isempty(data.thresh)
    for int = 1:size(data.sleepInt,2)
        hold on;
        plot([data.sleepInt(1,int) data.sleepInt(2,int)],[-.9 -.9],'-k')
    end
end

% Get threshold label
threshLabel = [get(handles.editPercent,'String'),'%'];
% Plot instantaneous amplitude with threshold bar in axes2
axes(handles.axes2)
plot(data.LFPTs.tvec,data.smInstAmp(chan,:),'r')
if ~isempty(data.thresh)
    hline(data.thresh(chan),'--',threshLabel)
end
title(['Instantaneous Amplitude of ',data.LFPTs.label{chan}])
% End busy indicator
set(gcf,'Pointer','arrow')
% --- Executes on button press in pushSave.
function pushSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data filepath

% Grab parameters for 'hist.sleep'
hist.sleep.thresh = data.thresh;
hist.sleep.percentile = data.percentile;
hist.sleep.onset = data.onset;
hist.sleep.offset = data.offset;

% Setup sleep event structure to be merged with eventTs
sleep.label = {'Sleep (Start)','Sleep (End)'};
sleep.t = {data.sleepInt(1,:)',data.sleepInt(2,:)'};
set(handles.debugbox,'string',filepath)
% Remove .mat ending
thisSave = filepath(1:end-4);
% Pick directory and save
uisave({'hist','sleep'},[thisSave,'_smoothed'])


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in setInts.
function setInts_Callback(hObject, eventdata, handles)
% hObject    handle to setInts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data

% Start 'busy' indicator
set(gcf,'Pointer','watch')
% Grab percentile, onset, and offset values
percent = str2double(get(handles.editPercent,'String'));
onset = str2double(get(handles.editOnset,'String'));
offset = str2double(get(handles.editOffset,'String'));

% for c = 1:size(data.LFPTs.data,1) 
%    % Get thresholds from percentile
%    thresh(c) = prctile(data.smInstAmp(c,data.goodInds),percent);
%    % Find indices above threshold
%    sleepInds{1,c} = data.LFPTs.tvec(data.smInstAmp(c,:) >= thresh(c));
%    % Skip to first index larger than onset to avoid negative indexing
%    firstIndex = find(sleepInds{1,c} > onset);
%    % Apply onset and offset to creat intervals
%    sleepInds{2,c}(1,:) = sleepInds{1,c}(firstIndex) - onset;
%    sleepInds{2,c}(2,:) = sleepInds{1,c}(firstIndex) + offset;
% end
% Only create intervals for percentiles below 100
if percent < 100
    for ct = 1:size(data.LFPTs.data,1)
        % Get thresholds from percentile for all channels
        thresh(ct) = prctile(data.smInstAmp(ct,data.goodInds),percent);
    end
    for c = get(handles.popupmenu1,'Value')
       % Find indices above threshold
       sleepInds{1,1} = data.LFPTs.tvec(data.smInstAmp(c,:) >= thresh(get(handles.popupmenu1,'Value')));
       % Skip to first index larger than onset to avoid negative indexing
       firstIndex = find(sleepInds{1} > onset);
       % Apply onset and offset to creat intervals
       sleepInds{2,1}(1,:) = sleepInds{1,1}(firstIndex) - onset;
       sleepInds{2,1}(2,:) = sleepInds{1,1}(firstIndex) + offset;
    end
    % Combine overlapping intervals
    [overInt] = intMaker(sleepInds(2,:));

    % Combinge all channels
    allChan = [];
    if size(overInt,1) >1
        for c = 1:size(overInt,1)
            allChan = horzcat(allChan,overInt{1,c}); 
        end
        % Sort by ascending order
        allChanSort{1,1} = sortrows(allChan',1)';
        % Re-merge intervals across all channels
        [overIntAll] = intMaker(allChanSort);
    else
        overIntAll = overInt;
    end
    % Check if last interval goes past end of recording and truncate
    if overIntAll{1,1}(2,end) > data.LFPTs.tvec(1,end)
        overIntAll{1,1}(2,end) = data.LFPTs.tvec(1,end);
    end
    data.sleepInt = overIntAll{1,1};
else
    thresh = [];
end
% Plot
% Get channel number from popupmenu1
chan = get(handles.popupmenu1,'Value');
% Use axes1 to plot raw signal
axes(handles.axes1)
cla
plot(data.LFPTs.tvec,data.LFPTs.data(chan,:)); 
title(data.LFPTs.label{chan}); xlabel('Time (secs)'); ylabel('Amplitude (V)');
% Plot rest and binge intervals if exist
% Find Binge Start and End indices
if isfield(data,'eventTs')
    indBs = find(not(cellfun('isempty',strfind(data.eventTs.label,'Binge (Start)'))));
    indBe = find(not(cellfun('isempty',strfind(data.eventTs.label,'Binge (End)'))));
    if ~isempty(indBs) && ~isempty(indBe)
        for ii = 1:size(data.eventTs.t{1,indBs},1)
            hold on;
            plot([data.eventTs.t{indBs,ii} data.eventTs.t{indBe,ii}],[-.9 -.9],'-r')
        end
    end
end

% Plot sleep intervals in blue below if they exist
if ~isempty(data.thresh)
    for int = 1:size(data.sleepInt,2)
        hold on;
        plot([data.sleepInt(1,int) data.sleepInt(2,int)],[-.88 -.88],'-bl')
    end
end

data.thresh = thresh;
data.onset = onset;
data.offset = offset;

% End busy indicator
set(gcf,'Pointer','arrow')

function editOnset_Callback(hObject, eventdata, handles)
% hObject    handle to editOnset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOnset as text
%        str2double(get(hObject,'String')) returns contents of editOnset as a double


% --- Executes during object creation, after setting all properties.
function editOnset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOnset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editOffset_Callback(hObject, eventdata, handles)
% hObject    handle to editOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOffset as text
%        str2double(get(hObject,'String')) returns contents of editOffset as a double


% --- Executes during object creation, after setting all properties.
function editOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushTest.
function pushTest_Callback(hObject, eventdata, handles)
% hObject    handle to pushTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data

% Grab percentile and channel
percent = str2double(get(handles.editPercent,'String'));
chan = get(handles.popupmenu1,'Value'); 
% Use percentile to get threshold
thresh = prctile(data.smInstAmp(chan,data.goodInds),percent);
threshLabel = [num2str(percent),'%'];
% Plot channel with percentile
axes(handles.axes2)
plot(data.LFPTs.tvec,data.smInstAmp(chan,:),'r')
hline(thresh,'--',threshLabel)
title(['Instantaneous Amplitude of ',data.LFPTs.label{chan}])
data.percentile = percent;
data.thresh = thresh;
set(handles.debugbox,'String',threshLabel)
