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

% Last Modified by GUIDE v2.5 07-Dec-2016 14:50:05

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

%load kz direction vector
kz_direction=str2num(get(handles.edit_kz_direction, 'String'));
kz_direction=kz_direction./norm(kz_direction);
kz_length=str2num(get(handles.edit_kz_value, 'String'));

if kz_direction==[0 0 1]
    [~,kz_index]=min(abs(raw_data.kz-kz_length));
    raw_data.kz=raw_data.kz(kz_index);
%     kz_cut_data=raw_data;
%     kz_cut_data.kz=kz_cut_data.kz(kz_index);
    for ii=1:raw_data.N_band
    raw_data.E{ii}=squeeze(raw_data.E{ii}(:,:,kz_index));
    end
    plane=createPlane(kz_length*kz_direction,kz_direction);
else
    % build list of 3D points forming 2D plane by defining
    % orthogonal vectors to g_hkl vector. Choose first vector as projection of
    % the z-axis to the new plane (only seems to work for other planes than 
    % 1 0 0, 0 1 0, 0 0 1 
    origin_plane=kz_length*kz_direction;
    plane=createPlane(origin_plane,kz_direction);
    kz_point=projPointOnPlane([0 0 1],plane);
    ky_2D_unit=round(kz_point-origin_plane,5);%needs better normalization
    ky_2D_unit=ky_2D_unit./norm(ky_2D_unit);
    kx_2D_unit=cross(kz_direction,ky_2D_unit);%needs normalization

    % generate meshgrid of coordinate vectors for new 2D plane system
    length_kz_cut_plane_side=str2num(get(handles.edit_length_kz_cut_plane_side, 'String'));
    resolution_cut=str2num(get(handles.edit_points_kz_cut_plane, 'String'));
    kx_ky_vectors=linspace(-length_kz_cut_plane_side,length_kz_cut_plane_side,resolution_cut);
    [X,Y]=meshgrid(kx_ky_vectors);
    origin_offset=ones(numel(X),1);

    % with the coordinate vectors and 3D basis vectors of plane, construct set
    % of 3D points that are evenly spaced on the plane
    temp=X(:)*kx_2D_unit+Y(:)*ky_2D_unit+origin_offset(:)*origin_plane;
    kx_2D=kx_ky_vectors*norm(kx_2D_unit);
    ky_2D=kx_ky_vectors*norm(ky_2D_unit);

    [KX,KY,KZ]=meshgrid(raw_data.kx,raw_data.ky,raw_data.kz);
    for ii=1:raw_data.N_band
        data_2D_plane=interp3(KX,KY,KZ, raw_data.E{ii},temp(:,1),temp(:,2),temp(:,3));  
        raw_data.E{ii}=reshape(data_2D_plane,length(kx_2D), length(ky_2D));
    end;
    raw_data.kx=kx_2D;
    raw_data.ky=ky_2D;
    raw_data.kz=kz_length;
end;
raw_data.kz_plane=plane;
%write 3D data to UserData
%set(handles.pushbutton_cut_kz,'UserData',{kz_cut_data});
assignin('base', 'bxsf_kzcut_data', raw_data);

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
color_switch_value=get(handles.popupmenu_contour_color,'Value');
color_switch_string_list=get(handles.popupmenu_contour_color,'String');
color=color_switch_string_list(color_switch_value);

show_text_switch=get(handles.radiobutton_contour_text_labels,'Value');
if show_text_switch
    show_text_switch='on';
else
    show_text_switch='off';
end

for ii=1:length(band_list_plotting)
    if color_switch_value==1
        contour(kz_cut_data.kx,kz_cut_data.ky,kz_cut_data.E{band_list_plotting(ii)},...
        contour_energies,'ShowText',show_text_switch);
    else
        contour(kz_cut_data.kx,kz_cut_data.ky,kz_cut_data.E{band_list_plotting(ii)},...
        contour_energies,'ShowText',show_text_switch,'Color',color{:});
    end
end;
axis equal





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


% --- Executes on button press in pushbutton_plot_cut.
function pushbutton_plot_cut_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_plot_cut=str2num(get(handles.edit_fig_plot_cuts,'String'));

% load 2D data from kz cut
bxsf_kzcut_data=evalin('base','bxsf_kzcut_data');

% select bands for plot
band_list_plotting_index=get(handles.listbox_select_bands,'Value');
band_list_plotting=cellfun(@str2num,get(handles.listbox_select_bands,'String'),'un',0);
band_list_plotting=cell2mat(band_list_plotting(band_list_plotting_index));

% load high sym points for cut
k_path = get(handles.uitable_k_path, 'data');
k_length=size(k_path);
k_length=k_length(1);

% interpolate cutting path
[X,Y]=meshgrid(bxsf_kzcut_data.kx,bxsf_kzcut_data.ky);

k_path_interp={};
interpolated_energy={};
s=0; %sets starting value for k-path length
l=1; %initializes running index
for ii=2:k_length
    if ~isnan(k_path(ii,1))
        x=linspace(k_path(ii-1,1),k_path(ii,1),51); %interpolate 100 point path between kx cooridnates
        y=linspace(k_path(ii-1,2),k_path(ii,2),51); %interpolate 100 point path between ky cooridnates
        k_path_interp_length{l}=norm([k_path(ii-1,1)-k_path(ii,1);k_path(ii-1,2)-k_path(ii,2)]); % measure length between points
        k_path_coordinates{l}=linspace(s,s+k_path_interp_length{l},51);
        s=s+k_path_interp_length{l};
        interpolated_energy{l}=smooth(interp2(X,Y,bxsf_kzcut_data.E{band_list_plotting},x,y),5);
        l=l+1;
    end
end;

no_high_sym_paths=length(interpolated_energy);
figure(fig_plot_cut)
hold on
for ii=1:no_high_sym_paths
%     subplot(1,no_high_sym_paths,ii)
    plot(k_path_coordinates{ii},interpolated_energy{ii})
end

a=5;


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fig_plot_cuts_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fig_plot_cuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fig_plot_cuts as text
%        str2double(get(hObject,'String')) returns contents of edit_fig_plot_cuts as a double


% --- Executes during object creation, after setting all properties.
function edit_fig_plot_cuts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fig_plot_cuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kz_direction_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kz_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kz_direction as text
%        str2double(get(hObject,'String')) returns contents of edit_kz_direction as a double


% --- Executes during object creation, after setting all properties.
function edit_kz_direction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kz_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_kz_plane.
function pushbutton_plot_kz_plane_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_kz_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bxsf_kzcut_plane=evalin('base','bxsf_kzcut_data.kz_plane');
drawPlane3d(bxsf_kzcut_plane);
alpha(0.2)



function edit_length_kz_cut_plane_side_Callback(hObject, eventdata, handles)
% hObject    handle to edit_length_kz_cut_plane_side (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_length_kz_cut_plane_side as text
%        str2double(get(hObject,'String')) returns contents of edit_length_kz_cut_plane_side as a double


% --- Executes during object creation, after setting all properties.
function edit_length_kz_cut_plane_side_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_length_kz_cut_plane_side (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_points_kz_cut_plane_Callback(hObject, eventdata, handles)
% hObject    handle to edit_points_kz_cut_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_points_kz_cut_plane as text
%        str2double(get(hObject,'String')) returns contents of edit_points_kz_cut_plane as a double


% --- Executes during object creation, after setting all properties.
function edit_points_kz_cut_plane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_points_kz_cut_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_2D_energy_surface.
function pushbutton_2D_energy_surface_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_2D_energy_surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_plot=str2num(get(handles.edit_plot_figure,'String'));

%load kz cut data from workspace
kz_cut_data=evalin('base','bxsf_kzcut_data');


%extract selected bands
band_list_plotting_index=get(handles.listbox_select_bands,'Value');
band_list_plotting=cellfun(@str2num,get(handles.listbox_select_bands,'String'),'un',0);
band_list_plotting=cell2mat(band_list_plotting(band_list_plotting_index));

%plot contours
figure(figure_plot)
hold on
for ii=1:length(band_list_plotting)

    surf(kz_cut_data.kx,kz_cut_data.ky,kz_cut_data.E{band_list_plotting(ii)})
end;
axis equal


% --- Executes on button press in pushbutton_translate_BZ.
function pushbutton_translate_BZ_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_translate_BZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load rawdata
raw_data=evalin('base','bxsf_data');

% find BZ boundary index
BZ_boundary=str2num(get(handles.edit_BZ_boundary,'String'));

[~,BZ_boundary_index]=min(abs(abs(raw_data.kx)-BZ_boundary));
BZ_boundary_index2=length(raw_data.kz)-BZ_boundary_index;

% cut index vector kx at boundary and then replicate twice along direction
raw_data.kx=raw_data.kx(BZ_boundary_index:BZ_boundary_index2);
raw_data.kx=[raw_data.kx-2*BZ_boundary,raw_data.kx,raw_data.kx+2*BZ_boundary];
raw_data.ky=raw_data.kx;
% raw_data.kz=raw_data.kx;

for ii=1:raw_data.N_band
    temp=raw_data.E{ii};
    temp=temp(BZ_boundary_index:BZ_boundary_index2,BZ_boundary_index:BZ_boundary_index2,...
        :); %cut ends of cube at BZ boundary
    temp=cat(1,temp,temp,temp);
    temp=cat(2,temp,temp,temp);
%     temp=cat(3,temp,temp,temp);
    raw_data.E{ii}=temp;
end;
% write symmetrized result into workspace
%set(handles.pushbutton_load,'UserData',{raw_data});
assignin('base', 'bxsf_data', raw_data);

set(handles.pushbutton_translate_BZ,'BackgroundColor','green'); 




function edit_BZ_boundary_Callback(hObject, eventdata, handles)
% hObject    handle to edit_BZ_boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_BZ_boundary as text
%        str2double(get(hObject,'String')) returns contents of edit_BZ_boundary as a double


% --- Executes during object creation, after setting all properties.
function edit_BZ_boundary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BZ_boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_translation_direction.
function popupmenu_translation_direction_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_translation_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_translation_direction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_translation_direction


% --- Executes during object creation, after setting all properties.
function popupmenu_translation_direction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_translation_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_delete_gca_contours.
function pushbutton_delete_gca_contours_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete_gca_contours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj(gca,'Type','contour');
delete(h(:));


% --- Executes on button press in pushbutton_delete_gca_surface.
function pushbutton_delete_gca_surface_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete_gca_surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj(gca,'Type','surface');
delete(h(:));


% --- Executes on selection change in popupmenu_contour_color.
function popupmenu_contour_color_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_contour_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_contour_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_contour_color


% --- Executes during object creation, after setting all properties.
function popupmenu_contour_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_contour_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_contour_text_labels.
function radiobutton_contour_text_labels_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_contour_text_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_contour_text_labels
