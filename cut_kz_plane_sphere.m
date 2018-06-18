function [ kz_plane_data ] = cut_kz_plane_sphere( raw4D_data,kz_direction, EplusV0,length_kz_cut_plane_side,resolution_cut )
%CUT_KZ_PLANE Summary of this function goes here
%   Detailed explanation goes here

% load 4D data
% 2*m_e/h_bar=constant_prefactor 
% A=constant_prefactor * 1 eV / 1e-20 for scaling to input E+V_0 in units
% of eV, and to get out k in units of inverse Angstrom
% A=2*9.10938291e-31*1.6e-19/1.054571726e-34^2*10^-20;
% kz_plane_intercept=sqrt(A.*EplusV0);
A=2*9.10938291e-31*1.6e-19/1.054571726e-34^2*10^-20;
kz_plane_intercept=sqrt(A*EplusV0);

% if norm(kz_direction-[0 0 1])<0.001
%     [~,kz_index]=min(abs(raw4D_data.kz-kz_plane_intercept));
%     raw4D_data.kz=raw4D_data.kz(kz_index);
% %     kz_cut_data=raw_data;
% %     kz_cut_data.kz=kz_cut_data.kz(kz_index);
%     for ii=1:raw4D_data.N_band
%     raw4D_data.E{ii}=squeeze(raw4D_data.E{ii}(:,:,kz_index));
%     end
%     plane=createPlane(kz_plane_intercept*kz_direction,kz_direction);
% else
    % build list of 3D points forming 2D plane by defining
    % orthogonal vectors to g_hkl vector. Choose first vector as projection of
    % the z-axis to the new plane (only seems to work for other planes than 
    % 1 0 0, 0 1 0, 0 0 1 
if norm(kz_direction-[0 0 1])<0.001
    ky_2D_unit=[0 1 0];
    kx_2D_unit=[1 0 0];
    kz_unit=kz_direction./norm(kz_direction);
else
    origin_plane=kz_plane_intercept*kz_direction;
    plane=createPlane(origin_plane,kz_direction);
    kz_unit=kz_direction./norm(kz_direction);
    kz_point=projPointOnPlane([0 0 1],plane);
    ky_2D_unit=round(kz_point-origin_plane,5);%needs better normalization
    ky_2D_unit=ky_2D_unit./norm(ky_2D_unit);
    kx_2D_unit=cross(kz_direction,ky_2D_unit);%needs normalization
end
    % generate meshgrid of coordinate vectors for new 2D plane system
    kx_vector=linspace(min(raw4D_data.kx),max(raw4D_data.kx),resolution_cut);
    ky_vector=linspace(min(raw4D_data.ky),max(raw4D_data.ky),resolution_cut);
    [X,Y]=meshgrid(kx_vector,ky_vector);
    k_X=X(:).*kx_2D_unit;
    k_Y=Y(:).*ky_2D_unit;
    
    k_X_norm=X(:).*norm(kx_2D_unit);
    k_Y_norm=Y(:).*norm(ky_2D_unit);
    
    Z=sqrt(A.*EplusV0-k_X_norm.^2-k_Y_norm.^2); %define z value such that points
    % form a sphere
    Z(abs(imag(Z))>0)=NaN;
    k_Z=Z(:).*kz_unit;

    % with the coordinate vectors and 3D basis vectors of plane, construct set
    % of 3D points that are evenly spaced on the plane
    temp=k_X+k_Y+k_Z;
%     kx_2D=kx_ky_vectors*norm(kx_2D_unit);
%     ky_2D=kx_ky_vectors*norm(ky_2D_unit);

    [KX,KY,KZ]=meshgrid(raw4D_data.kx,raw4D_data.ky,raw4D_data.kz);
    for ii=1:raw4D_data.N_band
        data_2D_plane=interp3(KX,KY,KZ, raw4D_data.E{ii},temp(:,1),temp(:,2),temp(:,3));  
        raw4D_data.E{ii}=reshape(data_2D_plane,length(kx_vector), length(ky_vector));
    end
%     raw4D_data.kx=kx_2D;
%     raw4D_data.ky=ky_2D;
    raw4D_data.kz=kz_plane_intercept;

kz_plane_data.kx=raw4D_data.kx;
kz_plane_data.ky=raw4D_data.ky;
kz_plane_data.kz_radius=raw4D_data.kz;
kz_plane_data.kz_sphere=[0 0 0 kz_plane_intercept];
kz_plane_data.E=raw4D_data.E;

end

