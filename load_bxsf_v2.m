function [Rawdata]=load_bxsf_v2()
%This function loads bxsf files according to Teng's conventions
% This file was copied from Teng's original program
[filename,pathname]=uigetfile('*.bxsf');
if isequal(filename,0)
    return;
end
fullpath=fullfile(pathname,filename);

%bxsf2mat
fid=fopen(fullpath,'r');
tline='t'; %sets tline to arbitrary char
Ef=[];
bandgrid_counter=0;
N=[];

while ischar(tline)
    tline = fgetl(fid);
    
    if ~ischar(tline)
        continue;
    else
        tline=strtrim(tline);
    end
        
    if isempty(tline) || tline(1)=='#' 
        continue
    end
    
    if isempty(Ef)
        k = strfind(tline, 'Fermi Energy:');
        if ~isempty(k)
            Ef = textscan(tline,'%s %s %f');
            Ef=Ef{1,3};
        end
    elseif isempty(N)
        
        k = strfind(tline, 'BANDGRID_3D');
        
        if ~isempty(k)
            bandgrid_counter=bandgrid_counter+1;
        end

        if bandgrid_counter==2 
        counter=1; %sets counter to read the next 6 meaningful lines
            while counter <= 6
                tline = fgetl(fid);
                if ~ischar(tline)
                    continue;
                else
                    tline=strtrim(tline);
                end

                if isempty(tline) || tline(1)=='#'
                    continue
                end

            switch counter
                case 1
                    N_band=str2num(tline);
                case 2
                    dummy=textscan(tline, '%f %f %f');
                    Nx=dummy{1};
                    Ny=dummy{2};
                    Nz=dummy{3};
                    N=Nx*Ny*Nz;
                case 3
                    dummy=textscan(tline, '%f %f %f');
                    G0=cell2mat(dummy);
                case 4
                    dummy=textscan(tline, '%f %f %f');
                    v1=cell2mat(dummy);
                case 5
                    dummy=textscan(tline, '%f %f %f');
                    v2=cell2mat(dummy);
                case 6
                    dummy=textscan(tline, '%f %f %f');
                    v3=cell2mat(dummy);

            end   

            counter=counter+1;

            end
            break
            
        end
    
    end
    

    
end

for j=1:N_band
    fgetl(fid);
    M=fscanf(fid,'%f',[1,N]);
    value{j,1}=permute(reshape(M,[Nz Ny Nx]),[3 2 1]);
    E_range(j,1)=min(M(:));
    E_range(j,2)=max(M(:));
    fgetl(fid);
    display(['loaded band number ',num2str(j)])
end

%  a=5;   
%   disp(tline)
  fclose(fid);

Rawdata.E=value;
Rawdata.E_range=E_range;
Rawdata.Ef=Ef;
Rawdata.N_band=N_band;
Rawdata.Nx=double(Nx);
Rawdata.Ny=double(Ny);
Rawdata.Nz=double(Nz);
Rawdata.G0=G0;
Rawdata.v1=v1;
Rawdata.v2=v2;
Rawdata.v3=v3;
end




% 
% for j=1:9
%     temp=fgetl(fid);
% end
% cell_tmp=textscan(temp,'%s%s%f','delimiter',' ','MultipleDelimsAsOne',1);
% Ef=cell_tmp{3};
% for j=1:6
%     temp=fgetl(fid);
% end
% cell_tmp=textscan(temp,'%d','delimiter',' ','MultipleDelimsAsOne',1);
% N_band=cell_tmp{1};
% temp=fgetl(fid);
% cell_tmp=textscan(temp,'%d%d%d','delimiter',' ','MultipleDelimsAsOne',1);
% Nx=cell_tmp{1};
% Ny=cell_tmp{2};
% Nz=cell_tmp{3};
% N=Nx*Ny*Nz;
% temp=fgetl(fid);
% cell_tmp=textscan(temp,'%f%f%f','delimiter',' ','MultipleDelimsAsOne',1);
% G0=cell2mat(cell_tmp);
% for j=1:3
%     temp=fgetl(fid);
%     cell_tmp=textscan(temp,'%f%f%f','delimiter',' ','MultipleDelimsAsOne',1);
%     v(j,1:3)=cell2mat(cell_tmp);
% end
% 
% tic;
% for j=1:N_band
%     fgetl(fid);
%     M=fscanf(fid,'%f',[1,N]);
%     value{j,1}=permute(reshape(M,[Nz Ny Nx]),[3 2 1]);
%     E_range(j,1)=min(M(:));
%     E_range(j,2)=max(M(:));
%     fgetl(fid);
% end
% fclose(fid);
% toc;
% 
% Rawdata.E=value;
% Rawdata.E_range=E_range;
% Rawdata.Ef=Ef;
% Rawdata.N_band=N_band;
% Rawdata.Nx=double(Nx);
% Rawdata.Ny=double(Ny);
% Rawdata.Nz=double(Nz);
% Rawdata.G0=G0;
% Rawdata.v1=v(1,:);
% Rawdata.v2=v(2,:);
% Rawdata.v3=v(3,:);
% end
