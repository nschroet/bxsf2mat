function varargout = quick_bxsf2mat(varargin)
% QUICK_BXSF2MAT MATLAB code for quick_bxsf2mat.fig
%      QUICK_BXSF2MAT, by itself, creates a new QUICK_BXSF2MAT or raises the existing
%      singleton*.
%
%      H = QUICK_BXSF2MAT returns the handle to a new QUICK_BXSF2MAT or the handle to
%      the existing singleton*.
%
%      QUICK_BXSF2MAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUICK_BXSF2MAT.M with the given input arguments.
%
%      QUICK_BXSF2MAT('Property','Value',...) creates a new QUICK_BXSF2MAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quick_bxsf2mat_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quick_bxsf2mat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quick_bxsf2mat

% Last Modified by GUIDE v2.5 30-Oct-2016 14:07:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @quick_bxsf2mat_OpeningFcn, ...
                   'gui_OutputFcn',  @quick_bxsf2mat_OutputFcn, ...
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


% --- Executes just before quick_bxsf2mat is made visible.
function quick_bxsf2mat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quick_bxsf2mat (see VARARGIN)

% Choose default command line output for quick_bxsf2mat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes quick_bxsf2mat wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = quick_bxsf2mat_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
interp_points=str2num(get(handles.edit_no_interp_points,'String'));
length_interp_vect=str2num(get(handles.edit_length_interp_vect,'String'));

% load and convert bxsf
rawdata_converted=bxsf2mat(load_bxsf,interp_points,length_interp_vect);

% write loaded data to workspace for other functions to access
%set(handles.pushbutton_load,'UserData',{rawdata_converted});
assignin('base', 'bxsf_data', rawdata_converted);

set(handles.pushbutton_load,'BackgroundColor','green'); 

% list bands in listbox
list_band_numbers=num2cell(1:rawdata_converted.N_band);
list_band_numbers_crossing_Ef=num2cell(rawdata_converted.band_numbers_crossing_Ef);

if get(handles.radiobutton_band_crossing,'Value')
    set(handles.listbox_select_bands,'String',list_band_numbers_crossing_Ef);
else
    set(handles.listbox_select_bands,'String',list_band_numbers);
end;




% --- Executes on button press in pushbutton_cut_kz.
function pushbutton_cut_kz_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cut_kz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load 4D data
%raw_data=get(handles.pushbutton_load,'UserData');
raw_data=evalin('base','bxsf_data');
%raw_data=raw_data{1};
kz_value=str2num(get(handles.edit_kz_value,'String'));
[~,kz_index]=min(abs(raw_data.kz-kz_value));
kz_cut_data=raw_data;
kz_cut_data.kz=kz_cut_data.kz(kz_index);

%cut 2D slice out of 3D energy data
for ii=1:kz_cut_data.N_band
    kz_cut_data.E{ii}=squeeze(kz_cut_data.E{ii}(:,:,kz_index));
 
end;

%write 3D data to UserData
%set(handles.pushbutton_cut_kz,'UserData',{kz_cut_data});
assignin('base', 'bxsf_kzcut_data', kz_cut_data);

set(handles.pushbutton_cut_kz,'BackgroundColor','green'); 

% --- Executes on selection change in listbox_select_bands.
function listbox_select_bands_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_select_bands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_select_bands contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_select_bands


% --- Executes during object creation, after setting all properties.
function listbox_select_bands_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_select_bands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kz_value_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kz_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kz_value as text
%        str2double(get(hObject,'String')) returns contents of edit_kz_value as a double


% --- Executes during object creation, after setting all properties.
function edit_kz_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kz_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_contour_plotting.
function pushbutton_contour_plotting_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_contour_plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_plot=str2num(get(handles.edit_plot_figure,'String'));

%load kz cut data from workspace
kz_cut_data=evalin('base','bxsf_kzcut_data');
% kz_cut_data=get(handles.pushbutton_cut_kz,'UserData');
% kz_cut_data=kz_cut_data{1};


%extract selected bands
band_list_plotting_index=get(handles.listbox_select_bands,'Value');
band_list_plotting=cellfun(@str2num,get(handles.listbox_select_bands,'String'),'un',0);
band_list_plotting=cell2mat(band_list_plotting(band_list_plotting_index));



%extract contour energies
contour_energies=sort(str2num(get(handles.edit_plot_energies, 'String')));
if length(contour_energies)==1
    contour_energies=[contour_energies contour_energies];
end;

%plot contours

figure(figure_plot)
hold on
for ii=1:length(band_list_plotting)
    contour(kz_cut_data.kx,kz_cut_data.ky,kz_cut_data.E{band_list_plotting(ii)},...
        contour_energies,'ShowText','on','Color','red')
end;





function edit_plot_energies_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_energies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_energies as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_energies as a double


% --- Executes during object creation, after setting all properties.
function edit_plot_energies_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_energies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_band_crossing.
function radiobutton_band_crossing_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_band_crossing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_band_crossing



function edit_plot_figure_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_figure as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_figure as a double


% --- Executes during object creation, after setting all properties.
function edit_plot_figure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_symmetrize.
function pushbutton_symmetrize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_symmetrize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load rawdata
raw_data=evalin('base','bxsf_data');
% raw_data=get(handles.pushbutton_load,'UserData');
% raw_data=raw_data{1}

% use symmetrize function
[ raw_data ] = symmetrize_mat( raw_data );

% write symmetrized result into workspace
%set(handles.pushbutton_load,'UserData',{raw_data});
assignin('base', 'bxsf_data', raw_data);

set(handles.pushbutton_symmetrize,'BackgroundColor','green'); 

a=5;
% do symmetrization


% --- Executes during object creation, after setting all properties.
function pushbutton_contour_plotting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_contour_plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_plot_3D_isosurface.
function pushbutton_plot_3D_isosurface_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_3D_isosurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_plot=str2num(get(handles.edit_fig_isosurface,'String'));
isosurface_energy=str2num(get(handles.edit_energy_isosurface, 'String'));

%load bxsf data
% isosurface_data=get(handles.pushbutton_load,'UserData');
% isosurface_data=isosurface_data{1};
isosurface_data=evalin('base','bxsf_data');

%extract selected bands
band_list_plotting_index=get(handles.listbox_select_bands,'Value');
band_list_plotting=cellfun(@str2num,get(handles.listbox_select_bands,'String'),'un',0);
band_list_plotting=cell2mat(band_list_plotting(band_list_plotting_index));
[X,Y,Z]=meshgrid(isosurface_data.ky, ...
    isosurface_data.kx, ...
    isosurface_data.kz);
figure(figure_plot)
hold on
color_list=['y','m','cyan', 'red', 'green', 'blue'];
for ii=1:length(band_list_plotting)
    fv = isosurface(X,Y,Z,isosurface_data.E{band_list_plotting(ii)},isosurface_energy);
    patch('Faces',fv.faces,'Vertices',fv.vertices,'FaceColor',color_list(ii))
end;


%extract contour energies


function edit_energy_isosurface_Callback(hObject, eventdata, handles)
% hObject    handle to edit_energy_isosurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_energy_isosurface as text
%        str2double(get(hObject,'String')) returns contents of edit_energy_isosurface as a double


% --- Executes during object creation, after setting all properties.
function edit_energy_isosurface_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_energy_isosurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fig_isosurface_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fig_isosurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fig_isosurface as text
%        str2double(get(hObject,'String')) returns contents of edit_fig_isosurface as a double


% --- Executes during object creation, after setting all properties.
function edit_fig_isosurface_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fig_isosurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_no_interp_points_Callback(hObject, eventdata, handles)
% hObject    handle to edit_no_interp_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_no_interp_points as text
%        str2double(get(hObject,'String')) returns contents of edit_no_interp_points as a double


% --- Executes during object creation, after setting all properties.
function edit_no_interp_points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_no_interp_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_length_interp_vect_Callback(hObject, eventdata, handles)
% hObject    handle to edit_length_interp_vect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_length_interp_vect as text
%        str2double(get(hObject,'String')) returns contents of edit_length_interp_vect as a double


% --- Executes during object creation, after setting all properties.
function edit_length_interp_vect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_length_interp_vect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
