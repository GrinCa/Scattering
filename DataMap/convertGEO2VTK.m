function convertGEO2VTK(FEmatrices,mesh,sizemesh,SOL,param,index)

%##########################################################################
%Please read the following lines for further informations
%##########################################################################
%this function aims to convert "SOL" array calculated with FE or WCAWE
%method into a .vtk file which Paraview can read. To understand how this
%scripts work, it is important to know what is the structure of .vtk file.
%You may find good info in the "vtk_fil_documentation.pdf" in the
%Documentation folder.
%We firstly need the Nodes file in order to get the coordinates of each
%node. Then we need SOL array (size = (i*ndof,nfreq), i=number of function). 
%It is possible thatwith the time, version compatibilities are no longer 
%working. The main probleme of this type of file(.vtk) is that Paraview 
%won't give you accurate information if it fails to read.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To sum up .vtk file contains(in this particular case):
% -header: gives the version of the vtk file
% -DATA_SET type : in our case type=UNSCTRUCTURED_GRID which enables us to
% choose as we want triangle elements for 2D or tetrahedral elements for 3D
% -POINTS : content of the Nodes file (coordinates)
% -CELLS : which refers to the .msh file. It contains the ID of the nodes
% of each elements
% -CELLS_TYPE : contains the ID of the element to use, according .vtk
% files. For instance triangles ID=5, and tetra ID=10.
% -POINT_DATA : can take several data, from scalar for pressure to tensor
% for stress field. In our case it is SCALAR id.


% FILENAME refers to .msh file and it enables us to create the path to
% store .vtk files
FILENAME = mesh.file;
connectivity_table = load(['Matrices/',FILENAME,'/connectivity_table.txt']);
connectivity_table = connectivity_table(:,[1 2 3 4 5 8 6 7 9 10]);
text_field = get_text_field(FEmatrices.Nodes,connectivity_table,FILENAME);
convertVTK3(FEmatrices,text_field,SOL,FILENAME,param,index,sizemesh)

end


function convertVTK3(FEmatrices,text_field,SOL,FILENAME,param,index,sizemesh)
% -This function create as much .vtk file as we have frequencies.
% -text field contains all the text data that every .vtk file needs,
% and which is the same no matter the frequency considerate. It was
% therefore useful to get it once and for all, print it in each
% file, instead of "recalculating" it every loop iteration...
% -index is given by the user at the call of convertGEO2VTK and contains
% the IDs of the frequencies (regarding freq array) that we want to write 
% in memory.
ndof = size(FEmatrices.Nodes,1);

indexfreq = index{1};
indextheta = index{2};

SOLUTION = zeros(ndof,length(indexfreq),length(indextheta));
SOLUTION(:,:,:) = real(SOL(:,:,:));



for ii=indexfreq
    for jj=indextheta
        disp(['***Converting [FREQ,THETA] = [',num2str(param.freq(ii)),',',num2str(180*param.theta(jj)/pi),']***']);
        Scattered_field = zeros(ndof,length(indexfreq),length(indextheta));
        Scattered_field(:,ii,jj) = SOLUTION(:,ii,jj);
        
        file_name = strcat('DataMap/',FILENAME,'/',FILENAME,'_sizemesh_',num2str(sizemesh),'_freq_',num2str(param.freq(ii)),'_theta_',num2str(180*param.theta(jj)/pi),'.vtk');

        fileID = fopen(file_name,'wt');
        fprintf(fileID,text_field);
        text_data = [];
        for kk=1:ndof
            text_data = [text_data [num2str(real(Scattered_field(kk,ii,jj))) '\n']];%Cavity_pressure(jj,ii)
        end
        text_data = [text_data 'SCALARS BACKGROUND_PRESSURE float 1\n'];
        text_data = [text_data 'LOOKUP_TABLE default\n'];
        for kk=1:ndof
            text_data = [text_data [num2str(real(FEmatrices.BG_pressure(kk,ii,jj))) '\n']];
        end
        fprintf(fileID,text_data);
        fclose(fileID);
    end
end
end


function text_field = get_text_field(Nodes,connectivity_table,FILENAME)
ndof = size(Nodes,1);
text_field = [];
% -text field contains all the text data that every .vtk file needs,
% and which is the same no matter the frequency considerate. It was
% therefore useful to get it once and for all, print it in each
% file, instead of "recalculating" it every loop iteration...
text_field = [text_field ['# vtk DataFile Version 2.0\n',FILENAME,'\nASCII\n']];
% above is the header refering to the version of vtk. It might have
% changed...
text_field = [text_field 'DATASET UNSTRUCTURED_GRID\n'];       %refer to doc
text_field = [text_field ['POINTS ',num2str(ndof),' float\n']];%refer to doc
% wrinting coordinates of each node
for ii=1:ndof
    text_field = [text_field [num2str(Nodes(ii,1)),' ',...
                              num2str(Nodes(ii,2)),' ',...
                              num2str(Nodes(ii,3)),'\n']];
end

disp('***Initialize conversion 3D***');
text_field = [text_field ['CELLS ',num2str(size(connectivity_table,1)),' ',...
                                   num2str(11*size(connectivity_table,1)),'\n']];
for ii=1:size(connectivity_table,1)
    text_field = [text_field ['10 ',num2str(connectivity_table(ii,:))],'\n'];
end 
text_field = [text_field ['CELL_TYPES ',num2str(size(connectivity_table,1)),'\n']];
for ii=1:size(connectivity_table,1)
    text_field = [text_field '24\n'];
end

text_field = [text_field ['POINT_DATA ' num2str(ndof) '\n']];
text_field = [text_field 'SCALARS SCATTERED_PRESSURE float 1\n'];
text_field = [text_field 'LOOKUP_TABLE default\n'];
end


