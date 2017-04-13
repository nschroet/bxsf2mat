function [mat_data]=bxsf2mat(bxsf_rawdata,no_interpolation_points,interpolation_length)
%bxsf2mat interpolates bxsf data (which is generally in a non-orthoginal)
%into orthogonal cartesian coordinates in inverse Angstroms
%   bxsf_rawdata is in the format specified by Teng's loader, as given by
%   load_bxsf

% build transformation matrix from bxsf vectors to cartesian vectors 
% [1 0 0] [0 1 0] [0 0 1]

matrix_bxsf_vect=[bxsf_rawdata.v1' bxsf_rawdata.v2' bxsf_rawdata.v3'];

%invert the matrix to get trafo from cartesian to bxfs vectos
trafo_matrix_cart_2_bxsf_vect_space=inv(matrix_bxsf_vect);

% create cartesian meshgrid
cartesian_length_vect=linspace(-interpolation_length,interpolation_length,no_interpolation_points); %here we use cube with side length 2
[X_1,Y_1,Z_1] = meshgrid(cartesian_length_vect,...
    cartesian_length_vect,...
    cartesian_length_vect);

% Create an N x 3 matrix of coordinates
points_cartesian = [X_1(:), Y_1(:), Z_1(:)];

% Transform the coordinates
points_transformed = points_cartesian * trafo_matrix_cart_2_bxsf_vect_space;

mat_data=bxsf_rawdata;
mat_data.kx=cartesian_length_vect;
mat_data.ky=cartesian_length_vect;
mat_data.kz=cartesian_length_vect;
mat_data.N_band_E_range=num2cell(mat_data.E_range,2);

% compile list of bands crossing Ef
band_numbers_crossing_Ef=[];
for ii=1:mat_data.N_band
    if (mat_data.Ef<max(mat_data.N_band_E_range{ii})) && (mat_data.Ef>min(mat_data.N_band_E_range{ii}))
    band_numbers_crossing_Ef=[band_numbers_crossing_Ef,ii];
    end;
end;
mat_data.band_numbers_crossing_Ef=band_numbers_crossing_Ef;


% clean fields from Teng's data 
mat_data=rmfield(mat_data,{'Nx','Ny','Nz','G0','E_range'});


%generate normalized meshgrid for interpolation (since we are in the
%non-orthogonal space, which is now treated as a cartesian cube)
cartesian_length_vect1=linspace(-1,1,2*bxsf_rawdata.Nx-1); %range goes from -1 to 1 since we extend to larger cube during interpolation
cartesian_length_vect2=linspace(-1,1,2*bxsf_rawdata.Ny-1); %note that -1 originates from not counting the boundary between cubes twice
cartesian_length_vect3=linspace(-1,1,2*bxsf_rawdata.Nz-1);

[X,Y,Z] = meshgrid(cartesian_length_vect1,...
    cartesian_length_vect2,...
    cartesian_length_vect3);
for ii=1:bxsf_rawdata.N_band
    %// interpolate the energies of the transformed cartesian coordinates.
    
    %copy the raw data that is defined over positive quadrant kx=[0
    %1],ky=[0 1],kz=[0 1] to full quadrant kx=[-1
    %1],ky=[-1 1],kz=[-1 1]. since bxsf data is using general grids, not
    %periodic ones (see section "bandgrids" in http://www.xcrysden.org/doc/XSF.html
    % we will need to avoid copying the redundant boundary of the data grid
    % twice, hence the subsequent 2:end parts of the code
    tmp=bxsf_rawdata.E{ii};
    tmp=cat(1,tmp,tmp(2:end,:,:));
    tmp=cat(2,tmp,tmp(:,2:end,:));
    tmp=cat(3,tmp,tmp(:,:,2:end));
    
    % now interpolation happens on enlarged grid
    data_cartesian=interp3(X,Y,Z, tmp, points_transformed(:,1), points_transformed(:,2), points_transformed(:,3));
    mat_data.E{ii}=reshape(data_cartesian,[no_interpolation_points,no_interpolation_points,no_interpolation_points]);
    mat_data.E{ii}=mat_data.E{ii}-mat_data.Ef;
end;
a=5;
end
