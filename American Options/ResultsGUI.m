function varargout = ResultsGUI(varargin)
% RESULTSGUI M-file for ResultsGUI.fig
%      RESULTSGUI, by itself, creates a new RESULTSGUI or raises the existing
%      singleton*.
%
%      H = RESULTSGUI returns the handle to a new RESULTSGUI or the handle to
%      the existing singleton*.
%
%      RESULTSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTSGUI.M with the given input arguments.
%
%      RESULTSGUI('Property','Value',...) creates a new RESULTSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ResultsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ResultsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ResultsGUI

% Last Modified by GUIDE v2.5 15-Aug-2007 18:14:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ResultsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ResultsGUI_OutputFcn, ...
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


% --- Executes just before ResultsGUI is made visible.
function ResultsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ResultsGUI (see VARARGIN)

% Choose default command line output for ResultsGUI
handles.output = hObject;

load AMERICAN_OPTION_DATA;
handles.data = data;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ResultsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

iPlot(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = ResultsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in CRR.
function CRR_Callback(hObject, eventdata, handles)
% hObject    handle to CRR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CRR

iPlot(hObject);

% --- Executes on button press in FD.
function FD_Callback(hObject, eventdata, handles)
% hObject    handle to FD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FD

iPlot(hObject);

% --- Executes on button press in LSM.
function LSM_Callback(hObject, eventdata, handles)
% hObject    handle to LSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LSM

iPlot(hObject);

function iPlot(h)

d = guidata(h);
cla;

v = view;
if v == view(2);
    v = 3;
end
if get(d.CRR,'value') == 1
    line(d.data.CRR.S,d.data.CRR.T,d.data.CRR.P,'color','b','linestyle','none','marker','.');
end

if get(d.FD,'value') == 1
    hold on
    surf(d.data.FD.S,d.data.FD.T,d.data.FD.P); shading interp;
end

if get(d.LSM,'value') == 1
    line(d.data.LSM.S,d.data.LSM.T,d.data.LSM.P,'color','r','linestyle','none','marker','.');
end
grid on
view(v);
title('American Option Price','fontsize',18);
xlabel('Share price','fontsize',18);
ylabel('Time','fontsize',18);
zlabel('Price of option','fontsize',18);
rotate3d('on');

