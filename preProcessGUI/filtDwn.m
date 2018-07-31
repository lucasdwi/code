function varargout = filtDwn(varargin)
% FILTGUI MATLAB code for filtGUI.fig
%      FILTGUI, by itself, creates a new FILTGUI or raises the existing
%      singleton*.
%
%      H = FILTGUI returns the handle to a new FILTGUI or the handle to
%      the existing singleton*.
%
%      FILTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILTGUI.M with the given input arguments.
%
%      FILTGUI('Property','Value',...) creates a new FILTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before filtDwn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to filtDwn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help filtGUI

% Last Modified by GUIDE v2.5 13-Jul-2018 15:18:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @filtDwn_OpeningFcn, ...
                   'gui_OutputFcn',  @filtDwn_OutputFcn, ...
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


% --- Executes just before filtGUI is made visible.
function filtDwn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to filtGUI (see VARARGIN)

% Choose default command line output for filtGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes filtGUI wait for user response (see UIRESUME)
% uiwait(handles.filtGUI);

% Set custom filter and dsf to uneditable
set(handles.lowFreq,'enable','off')
set(handles.highFreq,'enable','off')
set(handles.filterOrder,'enable','off')
set(handles.filterRipple,'enable','off')
set(handles.visFilter,'enable','off')
set(handles.dsf,'enable','off')

% Grab data
preProc = guidata(findobj('Tag','preProcGUI'));
% Auto populate channels
set(handles.channels,'String',preProc.LFPTs.label);

% --- Outputs from this function are returned to the command line.
function varargout = filtDwn_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function lowFreq_Callback(hObject, eventdata, handles)
% hObject    handle to lowFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowFreq as text
%        str2double(get(hObject,'String')) returns contents of lowFreq as a double


% --- Executes during object creation, after setting all properties.
function lowFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function highFreq_Callback(hObject, eventdata, handles)
% hObject    handle to highFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highFreq as text
%        str2double(get(hObject,'String')) returns contents of highFreq as a double


% --- Executes during object creation, after setting all properties.
function highFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filterOrder_Callback(hObject, eventdata, handles)
% hObject    handle to filterOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterOrder as text
%        str2double(get(hObject,'String')) returns contents of filterOrder as a double


% --- Executes during object creation, after setting all properties.
function filterOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filterRipple_Callback(hObject, eventdata, handles)
% hObject    handle to filterRipple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterRipple as text
%        str2double(get(hObject,'String')) returns contents of filterRipple as a double


% --- Executes during object creation, after setting all properties.
function filterRipple_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterRipple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in visFilter.
function visFilter_Callback(hObject, eventdata, handles)
% hObject    handle to visFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Grab data
preProc = guidata(findobj('Tag','preProcGUI'));
set(handles.filtGUI,'pointer','watch')
drawnow
low = str2num(get(handles.lowFreq,'String'));
high = str2num(get(handles.highFreq,'String'));
[b,a] = cheby1(2,.5,[low high]*2/preProc.adfreq,'stop');
fvtool(b,a,'Fs',preProc.adfreq)
set(handles.filtGUI,'pointer','arrow')
% --- Executes on button press in customFilter.
function customFilter_Callback(hObject, eventdata, handles)
% hObject    handle to customFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of customFilter
set(handles.lowFreq,'enable','on')
set(handles.highFreq,'enable','on')
set(handles.filterOrder,'enable','on')
set(handles.filterRipple,'enable','on')
set(handles.visFilter,'enable','on')

% --- Executes on button press in noFilter.
function noFilter_Callback(hObject, eventdata, handles)
% hObject    handle to noFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noFilter
set(handles.lowFreq,'enable','off')
set(handles.highFreq,'enable','off')
set(handles.filterOrder,'enable','off')
set(handles.filterRipple,'enable','off')
set(handles.visFilter,'enable','off')

% --- Executes on button press in defaultFilter.
function defaultFilter_Callback(hObject, eventdata, handles)
% hObject    handle to defaultFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of defaultFilter
set(handles.lowFreq,'enable','off')
set(handles.highFreq,'enable','off')
set(handles.filterOrder,'enable','off')
set(handles.filterRipple,'enable','off')
set(handles.visFilter,'enable','off')


% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Grab data
preProc = guidata(findobj('Tag','preProcGUI'));
% Change arrow to watch to indicate busy
set(handles.filtGUI,'pointer','watch')
drawnow
dsf = str2double(get(handles.dsf,'String'));
low = str2double(get(handles.lowFreq,'String'));
high = str2double(get(handles.highFreq,'String'));
ripple = str2double(get(handles.filterRipple,'String'));
order = str2double(get(handles.filterOrder,'String'));

% Run filter60.m
LFPTs.data = filter60(preProc.LFPTs.data,preProc.adfreq,0,'low',low,'high',high,'ripple',ripple,'order',order);
% Run dwnSample.m
LFPTs = dwnSample(preProc.LFPTs,dsf,preProc.adfreq);
% Change arrow to watch to indicate busy
set(handles.filtGUI,'pointer','arrow')
% Exit
close(filtDwn)
function dsf_Callback(hObject, eventdata, handles)
% hObject    handle to dsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dsf as text
%        str2double(get(hObject,'String')) returns contents of dsf as a double


% --- Executes during object creation, after setting all properties.
function dsf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultDwn.
function defaultDwn_Callback(hObject, eventdata, handles)
% hObject    handle to defaultDwn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of defaultDwn

set(handles.dsf,'enable','off')
% --- Executes on button press in customDwn.
function customDwn_Callback(hObject, eventdata, handles)
% hObject    handle to customDwn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of customDwn

set(handles.dsf,'enable','on')


function defaultDwn_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in psdPlot.
function psdPlot_Callback(hObject, eventdata, handles)
% hObject    handle to psdPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Grab data
preProc = guidata(findobj('Tag','preProcGUI'));
% Change arrow to watch to indicate busy
set(handles.filtGUI,'pointer','watch')
drawnow
% Get current channl
currentChan = get(handles.channels,'Value');
% Run pwelch
[Pxx,f] = pwelch(preProc.LFPTs.data(currentChan,:),2048,1024,1:100,preProc.adfreq);
% Plot
figure
plot(f,10*log10(Pxx))
xlabel('Frequency (Hz)')
ylabel('dB')
title(preProc.LFPTs.label{currentChan})
box off
% Change watch back to arrow
set(handles.filtGUI,'pointer','arrow')


% --- Executes on selection change in channels.
function channels_Callback(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channels


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
