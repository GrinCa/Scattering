function [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param)

isupdated = 0; % isupdated=1 : recalculate freefrem++ script

for ii=1:length(matrix_names)
    matrix_names{ii} = strcat(mesh.file,'/',matrix_names{ii});
    if exist(matrix_names(ii),"file")
    else
        isupdated=1;
    end
end


%--------------------------------------------------------------------------
%Run FreeFem++ script, IF NEEDED
%--------------------------------------------------------------------------

if (isupdated||flag.rerun) % EDP updated and/or mesh updated
   t_0 = cputime;
   disp('************************');
   disp('*Rerun FreeFem++ script*');
   disp('************************');
   edpcommand = strcat('FreeFem++'," ",mesh.file,'.edp');
   system(edpcommand);
   timing.freefem = cputime-t_0;
   disp('*********************************************************');
   output = sprintf('[Get_matrices:infos] CPUtime for building of matrices %.4f s',timing.freefem);
   disp(output);
   disp('*********************************************************');
   flag.recalculated = 1;
end % if

%--------------------------------------------------------------------------
% Get matrices from the files
%--------------------------------------------------------------------------

listLHS = cell(1,length(matrix_names)); 

tic;
% Matrices of the FE problem
for ii=1:length(matrix_names)
    fid = fopen(matrix_names(ii),'rt');
    for jj=1:3
        line = fgets(fid);
        if jj==3
           values = str2num(strtrim(line));
           nrows = values(1);
           ncolums = values(2);
        end
    end
    fclose(fid);
    matrix_data = importdata(matrix_names(ii)," ",3);
    matrix_data = matrix_data.data;
    %FEmatrices.listLHS{ii} = sparse(matrix_data(:,1)+1,matrix_data(:,2)+1,matrix_data(:,3));
    listLHS{ii} = sparse([matrix_data(:,1);nrows-1]+1,[matrix_data(:,2);ncolums-1]+1,[matrix_data(:,3);0]);
end

% Nodes
Nodes = load(strcat(mesh.file,'/',"Nodes.txt"));
ndof = size(Nodes,1);
toc;
%--------------------------------------------------------------------------
% return
%--------------------------------------------------------------------------

FEmatrices.Nodes = Nodes; 
FEmatrices = build_global(FEmatrices,listLHS,param,mesh.file);
end



function FEmatrices = build_global(FEmatrices,listLHS,param,FILENAME)

% definition of the label numbering
BGL_label = 1;
BGR_label = 2;

ndof = size(FEmatrices.Nodes,1);

H = listLHS{1}; %Stiffness matrix acoustic domain
Q = listLHS{2}; %Mass matrix acoustic domain
Hpmlr = listLHS{3}; %Stiffness matrix PML domain
Hpmli = listLHS{4};
Qpmlr = listLHS{5}; %Mass matrix PML domain
Qpmli = listLHS{6};

FEmatrices.H = H;
FEmatrices.Q = Q; % Q may be useful for Mean Quadratic Pressure calculation

% label of the different region of the mesh
region_labels = load(['Matrices/',FILENAME,'/labels.txt']); 

% initialisation arrays of respective nodes
acoustic_nodes = zeros(ndof,1);
BGL_nodes = zeros(ndof,1);
BGR_nodes = zeros(ndof,1);


for ii=1:length(region_labels)
    if region_labels(ii,2) == BGL_label || region_labels(ii,2) == BGR_label
        acoustic_nodes(region_labels(ii,1)+1) = 1;
    end
    if region_labels(ii,2) == BGL_label
        BGL_nodes(region_labels(ii,1)+1) = 1;
    elseif region_labels(ii,2) == BGR_label
        BGR_nodes(region_labels(ii,1)+1) = 1;
    end
end

acoustic_nodes = find(acoustic_nodes);
BGL_nodes = find(BGL_nodes);
BGR_nodes = find(BGR_nodes);
wall_nodes = find(abs(FEmatrices.Nodes(:,3))<(1e-10));

% save arrays of Nodes, needed for partionning
FEmatrices.acoustic_nodes = acoustic_nodes;
FEmatrices.BGL_nodes = BGL_nodes;
FEmatrices.BGR_nodes = BGR_nodes;
FEmatrices.wall_nodes = wall_nodes;


% indexing of the differents subspaces for partitionning
FEmatrices.indexp = acoustic_nodes;

Hpml = Hpmlr+1i*Hpmli;
Qpml = Qpmlr+1i*Qpmli;

Kglob = H + Hpml;
Mglob = Q + Qpml;

FEmatrices.LHS = {Kglob,Mglob};
% size of the "reduced" system < 4*ndof
FEmatrices.size_system = size(FEmatrices.LHS{1},1);
FEmatrices = get_RHS(FEmatrices,param);

% Verifcation of the global matrices
if ( size(Kglob,1)-length(find(diag(Kglob))) ) ~= 0
    disp('Stiffness matrix is singular');
elseif ( size(Mglob,1)-length(find(diag(Mglob))) ) ~= 0
    disp('Mass matrix is singular');
end

end

function FEmatrices = get_RHS(FEmatrices,param)
% Cell of RHS against f and theta
RHS_BG = cell(param.nfreq,param.ntheta);

BG_nodes = FEmatrices.wall_nodes;

xbg = FEmatrices.Nodes(BG_nodes,1);
zbg = FEmatrices.Nodes(BG_nodes,3);

FEmatrices.BG_pressure = zeros(FEmatrices.size_system,length(param.freq),length(param.theta));

P0 = 1;

for ii=1:param.nfreq
    for jj=1:param.ntheta
        U_inc = zeros(FEmatrices.size_system,1);
        k = 2*pi*param.freq(ii)/param.c0;
        BG_Pressure_tmp = P0*exp(-1i*k*(xbg*cos(param.theta(jj))+zbg*sin(param.theta(jj))));%BG_Pressure_tmp = P0*cos(2*pi*param.freq(ii)/param.c0*( xbg*cos(param.theta(jj))+zbg*sin(param.theta(jj)) )); % the array should just contain 1 component
        U_inc(BG_nodes,1) = BG_Pressure_tmp;
        FEmatrices.BG_pressure(:,ii,jj) = U_inc;
        Z = FEmatrices.H - (2*pi*param.freq(ii)/param.c0)^2*FEmatrices.Q;
        RHS_BG{ii,jj} = Z*U_inc;
    end
end

FEmatrices.RHS_BG = RHS_BG;

% % RHS
% 
% [~,node_zeros] = min(FEmatrices.Nodes(:,1).*FEmatrices.Nodes(:,1) + ...
%                      FEmatrices.Nodes(:,2).*FEmatrices.Nodes(:,2) + ...
%                      FEmatrices.Nodes(:,3).*FEmatrices.Nodes(:,3));
% RHStmp = zeros(FEmatrices.size_system,1);
% RHStmp(node_zeros) = 1;
% FEmatrices.RHS = RHStmp;

end

































