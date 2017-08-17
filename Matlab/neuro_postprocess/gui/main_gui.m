function varargout = main_gui(varargin)
% MAIN_GUI MATLAB code for main_gui.fig
%      MAIN_GUI, by itself, creates a new MAIN_GUI or raises the existing
%      singleton*.
%
%      H = MAIN_GUI returns the handle to a new MAIN_GUI or the handle to
%      the existing singleton*.
%
%      MAIN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_GUI.M with the given input arguments.
%
%      MAIN_GUI('Property','Value',...) creates a new MAIN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_gui

% Last Modified by GUIDE v2.5 05-Dec-2016 21:52:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_gui_OpeningFcn, ...
    'gui_OutputFcn',  @main_gui_OutputFcn, ...
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
end

% --- Executes just before main_gui is made visible.
function main_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_gui (see VARARGIN)

% Choose default command line output for main_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global gui_data

% Load brain mesh
load('data/brain_coordinates_simple.mat','fv');
gui_data.fv = fv;

% Load 3D coordinates of electrodes
load('data/cap_coordinates_3D.mat','cap_coordinates_3D');
gui_data.cap_coordinates_3D = cap_coordinates_3D;

% Load data of time-series
load('data/input_for_clustering_problem_from_Susanne.mat','x');

R = size(cap_coordinates_3D,1); % = 64, number of electrodes
gui_data.R = R;

T = length(x)/R;
gui_data.T = T;

% add some data to electrodes, each timestep stored in column
electrode_data = reshape(x,R,T);
gui_data.electrode_data = electrode_data;

%% plot axes_3D
axes(handles.axes_3D);

% initial plot, returns pointer to brain object
brain_patch = plot_brain(fv.vertices,fv.faces,zeros(size(fv.vertices)));
gui_data.brain_patch = brain_patch;

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

colorbar

axis off
view([-135 35]);

set(gcf,'Renderer','OpenGL');
set(gcf,'RendererMode','manual');

% plot t=1
t = 1;
gui_data.t = t;

% map electrode data to brain vertex data
brain_data = map_electrode_to_brain( fv.vertices, cap_coordinates_3D, electrode_data(:,t) );

% update brain coloring
update_brain(brain_patch, brain_data);
    
% set color axis
caxis([-100 100])
    
hold off

%% set timer text
set(handles.text3,'String',['T = ' num2str(t) ' (' num2str(T) ')' ]);

%% set slider properties
set(handles.time_slider,'Value',t);
set(handles.time_slider,'Min',1);
set(handles.time_slider,'Max',T);
set(handles.time_slider,'SliderStep',[1/(T-1) 1/(T-1)]);

end

% --- Outputs from this function are returned to the command line.
function varargout = main_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

% --- Executes on slider movement.
function time_slider_Callback(hObject, eventdata, handles)
% hObject    handle to time_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global gui_data

% round slider value
t = round(get(hObject,'Value'));
set(hObject,'Value',t);
gui_data.t = t;

%% set timer text
set(handles.text3,'String',['T = ' num2str(t) ' (' num2str(gui_data.T) ')' ]);

% update axis_3D
% map electrode data to brain vertex data
brain_data = map_electrode_to_brain( gui_data.fv.vertices, gui_data.cap_coordinates_3D, gui_data.electrode_data(:,t) );

% update brain coloring
update_brain(gui_data.brain_patch, brain_data);

end

% --- Executes during object creation, after setting all properties.
function time_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


function pause_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pause_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pause_edit as text
%        str2double(get(hObject,'String')) returns contents of pause_edit as a double

end

% --- Executes during object creation, after setting all properties.
function pause_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pause_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
