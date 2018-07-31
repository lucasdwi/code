function varargout = trializeGUI(varargin)
% TRIALIZEGUI MATLAB code for trializeGUI.fig
%      TRIALIZEGUI, by itself, creates a new TRIALIZEGUI or raises the existing
%      singleton*.
%
%      H = TRIALIZEGUI returns the handle to a new TRIALIZEGUI or the handle to
%      the existing singleton*.
%
%      TRIALIZEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIALIZEGUI.M with the given input arguments.
%
%      TRIALIZEGUI('Property','Value',...) creates a new TRIALIZEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trializeGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trializeGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trializeGUI

% Last Modified by GUIDE v2.5 16-Jul-2018 12:57:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trializeGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @trializeGUI_OutputFcn, ...
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

% --- Executes just before trializeGUI is made visible.
function trializeGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trializeGUI (see VARARGIN)

% Choose default command line output for trializeGUI
handles.output = hObject;
% Load gui1 data
preProc = guidata(findobj('Tag','preProcGUI'));
if ~isempty(preProc.eoi{1,1})
    set(handles.eventMenu,'String',preProc.eoi{:,1})
    handles.eoi = preProc.eoi;
else
    set(handles.eventMenu,'String',' ');
    handles.eoi = cell(1,1);
end

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes trializeGUI wait for user response (see UIRESUME)
% uiwait(handles.trialGUI);


% --- Outputs from this function are returned to the command line.
function varargout = trializeGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in eventMenu.
function eventMenu_Callback(hObject, eventdata, handles)
% hObject    handle to eventMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns eventMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eventMenu


% --- Executes during object creation, after setting all properties.
function eventMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in addEvent.
function addEvent_Callback(hObject, eventdata, handles)
% hObject    handle to addEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

run('trializeCfgGUI')


% --- Executes on button press in editEvent.
function editEvent_Callback(hObject, eventdata, handles)
% hObject    handle to editEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
