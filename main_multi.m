%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Test case using WCAWE                             %
%                                                                         %
%                             March 2020                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% Init main program
%--------------------------------------------------------------------------

clear FEmatrices param flag;
warning('off', 'MATLAB:nearlySingularMatrix');

%--------------------------------------------------------------------------
% Folders
%--------------------------------------------------------------------------

% Add folders for Mumps, WCAWE, and Mesh functions/files
global Meshfold WCAWEfold Mumps DataMap Derivatives
Meshfold = 'Matrices';
WCAWEfold = 'WCAWE';
Mumps = 'Mumps';
DataMap = 'DataMap';
Derivatives = 'Derivatives';


addpath(genpath(strcat(pwd,'/',Meshfold)));
addpath(genpath(strcat(pwd,'/',WCAWEfold)));
addpath(genpath(strcat(pwd,'/',Mumps)));
addpath(genpath(strcat(pwd,'/',DataMap)));
addpath(genpath(strcat(pwd,'/',Derivatives)));




%--------------------------------------------------------------------------
% Input data for the problem
%--------------------------------------------------------------------------

% Input parameters for Matlab calculation
flag.rerun = 0; % to recalculate FreeFem++ matrices
flag.recalculated = 0; % allow WCAWE and/or FE recalculation
flag.calculateFE = 1;  % calculate FE solution
flag.calculateWCAWE = 0; % calculate WCAWE solution

flag.plotcomparison = 0; % plot comparison between FE and WCAWE
flag.comparisonMULTI = 0;

flag.converge = 0;
flag.convert2VTK = 1; % convert SOLFE.mat into a .vkt file
flag.plotMQP = 0;

flag.getmatrices = 1;

if flag.converge || flag.plotMQP || flag.convert2VTK
    flag.getmatrices = 0;
end



% Input files for mesh and FE matrices
mesh.file = 'Wallv2';
sizemesh = load('sizemesh.txt');
sizemesh = sizemesh(end);

% define timing struct to access to time calculation of each method                                                    
timing.freefem = 0;
timing.WCAWE = 0;
timing.computeFE = 0;                                                    


% Material parameters
param.rho = 1.213;
param.rhoS = 1188;
param.c0 = 342.2;


% Frequency range
param.fmin = 400; % 300
param.fmax = 400; % 6000
param.f_range = [param.fmin param.fmax];
param.freqincr = 100; % 20
param.freq = param.fmin : param.freqincr : param.fmax; % frequency range
param.nfreq = length(param.freq);

% Angle range
param.thetamin = pi/4; % 300
param.thetamax = pi/4; % 6000
param.theta_range = [param.thetamin param.thetamax];
param.thetaincr = 0.1; % 20
param.theta = param.thetamin : param.thetaincr : param.thetamax; % frequency range
param.ntheta = length(param.theta);


% those frequencies are the frequencies point for Padé expension
param.freqref = [400];
param.nfreqref = length(param.freqref);

param.thetaref = [pi/4];
param.nthetaref = length(param.freqref);

% Input data for the loop over expansion orders. Note that for each
% frequency sweep the number of vectors from the WCAWE basis will depend on
% the number of point for Padé expension. For instance, if we have 2 points
% for expansion, and nvecfreq=5 (order of expansion), we will have 15
% vectors in the basis, 5 by intervals.
param.nvecfreqmin = 25;
param.nvecfreqmax = 25;
param.incrvec = 20;
param.vecfreqrange = param.nvecfreqmin : param.incrvec : param.nvecfreqmax;

param.nvecthetamin = 25;
param.nvecthetamax = 25;
param.vecfreqrange = param.nvecthetamin : param.incrvec : param.nvecthetamax;

%Identificator
param.idData = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']['...
                num2str(180*param.theta_range(1)/pi) '_' num2str(180*param.theta_range(2)/pi) ']'];

% generation of the different folder to store data if they don't already exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genfolders(mesh,param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% Build intervals for the Parametric sweep
%--------------------------------------------------------------------------

param = build_interval(param);



%--------------------------------------------------------------------------
% Matrices calculated with Freefem++
%--------------------------------------------------------------------------
if flag.getmatrices
    matrix_names = ["H.txt","Q.txt",... % matrices defined on the incident domain
                    "Hpmlr.txt","Hpmli.txt",...
                    "Qpmlr.txt","Qpmli.txt"];

    [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param);
    Nodes = FEmatrices.Nodes;
    LHS = FEmatrices.LHS;
    nLHS = length(LHS);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeff_LHS = {@(f) 1,@(f) -(2*pi*f/param.c0)^2};
coeff_RHS = @(f) 1;

if flag.recalculated
    
%--------------------------------------------------------------------------
% Calculate reference finite element parametric sweep
%--------------------------------------------------------------------------

t_0 = cputime;
if flag.calculateFE == 1
   disp('#######################');
   disp('Recalculate FE solution');
   disp('#######################');
   SOLFE = zeros(FEmatrices.size_system,param.nfreq,param.ntheta); %size ndof,nfreq
   % Parametric loop calculation
   for ii=1:param.nfreq
       for jj=1:param.ntheta
           tic;
           disp(['[FE] Frequency : ',num2str(param.freq(ii))]);
           Aglob = sparse(size(LHS{1},1),size(LHS{1},2));
           for kk = 1:nLHS
              Aglob = Aglob + coeff_LHS{kk}(param.freq(ii))*LHS{kk};
           end
           resp_P = Aglob\FEmatrices.RHS_BG{ii,jj};
           SOLFE(:,ii,jj) = resp_P;
           toc;
       end
   end
   timing.computeFE = cputime-t_0;
end



%--------------------------------------------------------------------------
% Calculate WCAWE parametric sweep
%--------------------------------------------------------------------------
if flag.calculateWCAWE == 1
    
    for nvecfreq=param.vecfreqrange

        % Initialize array of RHSderiv, Cell array linear comb coefficients derivative functions
        deriv_deg = [param.nvecfreqmax];

        if exist('Derivatives/derivative_orders.mat','file') ~= 2
            disp('#################################');
            disp('Recalculate all cross derivatives');
            disp('#################################');
            create_cross_derivatives(LHS,coeff_LHS,...
                                     coeff_RHS,deriv_deg,'f');
        end
        load('Derivatives/derivative_orders.mat');
        if ~isempty(find(derivative_orders-deriv_deg<0))
            disp('#################################');
            disp('Recalculate all cross derivatives');
            disp('#################################');
            create_cross_derivatives(LHS,coeff_LHS,...
                                     coeff_RHS,deriv_deg,'f');
        end 



        coeff_deriv_fun = cell(nLHS,nvecfreq+1);
        RHScoeffderiv = cell(1,nvecfreq+1);
        [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices(...
                                coeff_deriv_fun,RHScoeffderiv,nvecfreq,nLHS);

        coeff_derivgen_fun = @(freq) cellfun(@(cellfunc) cellfunc(freq),coeff_deriv_fun);

        % RHS

        RHSderivmulti = cell(1,param.nfreqref);
        for ii=1:param.nfreqref
            RHSderiv = cell(1,nvecfreq+1);
            for kk=1:nvecfreq+1
                RHSderiv{kk} = RHScoeffderiv{kk}(param.freqref(ii))*RHS;    %derivatives of RHS at freq0
            end
            RHSderivmulti{ii}=RHSderiv;
        end



        % Fill Cell array of linear comb coefficients derivatives at
        % freqref(ii)
        coeff_deriv_multi = cell(1,param.nfreqref);
        for ii=1:param.nfreqref
            coeff_deriv_multi{ii}=coeff_derivgen_fun(param.freqref(ii));
            % Fix coeff_deriv by replacing all NaN values by 0
            tmpidxnan = find(isnan(coeff_deriv_multi{ii}));
            coeff_deriv_multi{ii}(tmpidxnan) = 0;
            tmpidxinf = find(isinf(coeff_deriv_multi{ii}));
            coeff_deriv_multi{ii}(tmpidxinf) = 0;
        end




           %-----------------------------------------------------------------------
           %Recalculation of WCAWE basis
           %-----------------------------------------------------------------------


               %if exist(strcat('Matrices/',mesh.file,'/','[',num2str(param.f_range),']/',...
               %                'SOLWCAWE_',mesh.file,'_nvec',num2str(param.n_sub_range*nvecfreq),...
               %                '_[',num2str(param.freqref),']','.mat'),'file')~=2
               Wtrans = [];
               for ii=1:param.nfreqref
                  [Wtranstmp,Ucoeff,timing] = WCAWE_basis(LHS,coeff_deriv_multi{ii},RHSderivmulti{ii},nvecfreq,timing);
                   Wtrans = [Wtrans Wtranstmp];
               end
               [uu,vv,ww] = svd(Wtrans,0);
               iiselect = find(diag(vv)>vv(1,1)*1e-15);
               Wtranssvd = uu(:,iiselect);
               nsvd = size(Wtranssvd,2);
               output = sprintf("[SVD:Info] Number of selcted vector %d/%d",nsvd,size(Wtrans,2));
               disp(output);
               [SOLWCAWE] = Solve_WCAWE(LHS,coeff_deriv_fun,RHS,Wtranssvd,param.freq);
               %end
    end
end
end
%--------------------------------------------------------------------------
% Saves
%--------------------------------------------------------------------------
% FE solution
if flag.recalculated
    if flag.calculateFE
        save(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_sizemesh_',num2str(sizemesh),'.mat'],'SOLFE');
    end

    if flag.calculateWCAWE
        save(['Matrices/',mesh.file,'/',param.idData,'/SOLWCAWE_',mesh.file,'_nvec_',num2str(param.n_sub_range*nvecfreq),'_sizemesh_',num2str(sizemesh),'.mat'],'SOLWCAWE');
    end
    
    % save FEmatrices which contains all the data of the simulation for each
    % mesh
    save(['Matrices/',mesh.file,'/',param.idData,'/','DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param');
end


%--------------------------------------------------------------------------
% Post processing
%--------------------------------------------------------------------------

if flag.converge
    clear FEmatrices SOLFE SOLWCAWE;
    sizemesh_ARRAY = load('sizemesh.txt');
    fid1 = fopen('converge1.txt','wt');
    fid2 = fopen('converge2.txt','wt');
    if flag.calculateFE == 1
        for ii=1:length(sizemesh_ARRAY)
            sizemeshtmp = sizemesh_ARRAY(ii);    
            SOL = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_sizemesh_',num2str(sizemesh),'.mat']));
            SOL = SOL{1};
        end
    elseif flag.calculateWCAWE
        for ii=1:length(sizemesh_ARRAY)
            sizemeshtmp = sizemesh_ARRAY(ii);
            SOL = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLWCAWE_',mesh.file,'_nvec_',num2str(param.n_sub_range*nvecfreq),'_sizemesh_',num2str(sizemesh),'.mat']));
            SOL = SOL{1};
        end
    end
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    ndof = size(FEmatrices.Nodes,1);
    ndof_acoustic = length(FEmatrices.indexp);
    %ndof_acoustic = length(FEmatrices.acoustic_nodes);
    acoustic_volume = load('acoustic_volume.txt');
    Pressure = zeros(ndof,param.nfreq);
    Pressure(FEmatrices.acoustic_nodes,:) = real(SOL(FEmatrices.indexp,:));
    MQP1 = Pressure(FEmatrices.acoustic_nodes,:)'*FEmatrices.Q*Pressure(FEmatrices.acoustic_nodes,:)/((4e-10)*acoustic_volume);
    MQP1 = 10*log10(MQP1);
    fprintf(fid1,strcat(num2str(ndof),'\t',num2str(real(diag(MQP1)')),'\n'));
    MQP2 = Pressure'*Pressure/4/(10e-10)/2/ndof_acoustic;
    MQP2 = 10*log10(MQP2);
    fprintf(fid2,strcat(num2str(ndof),'\t',num2str(real(diag(MQP2)')),'\n'));
    fclose(fid1);
    fclose(fid2);
end

if flag.plotMQP
    fid = fopen('converge.txt','rt');
    while true
        line = fgets(fid);
        if line == -1
            break;
        else
            line = str2num(strtrim(line));
            plot(param.freq,line(2:end),'DisplayName',['ndof = ' num2str(line(1))]);
            hold on
        end
    end
    
    xlabel('Frequency (Hz)');
    ylabel('Mean quadratic pressure (dB)');
    legend();
    hold off
    
    fclose(fid);
end


%--------------------------------------------------------------------------
% convert
%--------------------------------------------------------------------------
if flag.convert2VTK
    clear FEmatrices;
    sizemesh_ARRAY = load('sizemesh.txt');
    sizemesh = sizemesh_ARRAY(end);
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    if flag.calculateWCAWE
        SOLWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLWCAWE_',mesh.file,'_nvec_',num2str(param.n_sub_range*nvecfreq),'_sizemesh_',num2str(sizemesh),'.mat']));
        SOL = SOLWCAWE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLWCAWE,param,{1:1:param.nfreq, 1:1:param.ntheta})
    elseif flag.calculateFE
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',param.idData,'/SOLFE_',mesh.file,'_sizemesh_',num2str(sizemesh),'.mat']));
        SOLFE = SOLFE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLFE,param,{1:1:param.nfreq, 1:1:param.ntheta})
    end
end














