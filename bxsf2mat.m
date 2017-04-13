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
cartesian_length_vect=linspace(0,interpolation_length,no_interpolation_points); %here we use cube with side length 2
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
mat_data=rmfield(mat_data,{'Nx','Ny','Nz','G0','v1','v2','v3','E_range'});


%generate normalized meshgrid for interpolation (since we are in the
%non-orthogonal space, which is now treated as a cartesian cube)
cartesian_length_vect1=linspace(0,1,bxsf_rawdata.Nx);
cartesian_length_vect2=linspace(0,1,bxsf_rawdata.Ny);
cartesian_length_vect3=linspace(0,1,bxsf_rawdata.Nz);

[X,Y,Z] = meshgrid(cartesian_length_vect1,...
    cartesian_length_vect2,...
    cartesian_length_vect3);
for ii=1:bxsf_rawdata.N_band
    %// interpolate the energies of the transformed cartesian coordinates
    data_cartesian=interp3(X,Y,Z, bxsf_rawdata.E{ii}, points_transformed(:,1), points_transformed(:,2), points_transformed(:,3));
    mat_data.E{ii}=reshape(data_cartesian,[no_interpolation_points,no_interpolation_points,no_interpolation_points]);
    mat_data.E{ii}=mat_data.E{ii}-mat_data.Ef;
end;
a=5;
end
