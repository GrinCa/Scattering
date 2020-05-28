function genfolders(mesh,param)

Filename = mesh.file;

%--------------------------------------------------------------------------
% Matrices
%--------------------------------------------------------------------------


if exist(['Matrices/',Filename]) == 0
    command = ['Matrices/',Filename];
    system(['mkdir ' command]);
end

if exist(['Matrices/',Filename,'/',param.idData]) == 0
    command = ['Matrices/',Filename,'/',param.idData];
    system(['mkdir ' command]);
end

%--------------------------------------------------------------------------
% Geometry : copy .msh from Geometry/filename into the main folder
%--------------------------------------------------------------------------

try
    command = ['cp ','Geometry/',Filename,'/',Filename,'.msh',' ',Filename,'.msh'];
    system(command);
catch
    disp(['[genfolders] No .msh file in the folder : /Geometry/',Filename]);
    return;
end
%--------------------------------------------------------------------------
% DataMap
%--------------------------------------------------------------------------

if exist(['DataMap/',Filename]) == 0
    command = ['DataMap/',Filename];
    system(['mkdir ' command]);
end

end