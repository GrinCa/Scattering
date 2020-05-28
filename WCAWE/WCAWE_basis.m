function [Wtrans,Ucoeff,timing] = WCAWE_basis(LHS,coeff_deriv,RHSderiv,nvecfreq,timing)

disp('#######################');
disp('Calculating WCAWE basis');
disp('#######################');


tic;

ndof = size(LHS{1},1);
nmatglob = length(LHS);

%warning off all

% ndof = 100;
% nvecfreq = 4;

Vtilde = zeros(ndof,nvecfreq);
V = zeros(ndof,nvecfreq);
U = [];



%--------------------------------------------------------------------------
%Initialisation of random cell matrices
%--------------------------------------------------------------------------


% RHSderiv = cell(1,4);
% Kglob = rand(ndof);
% Mglob = rand(ndof);
% listLHS = {Kglob,Mglob};
% Aglob = rand(ndof,ndof);
% nmatglob = length(listLHS);
% coeff_deriv = rand(2,nvecfreq);

% for ii=1:nvecfreq
%     Vtilde(:,ii) = zeros(ndof,1);
%     V(:,ii) = zeros(ndof,1);
%     %RHSderiv{ii} = rand(ndof,1);
% end

Aglob = sparse(ndof,ndof);
for kk=1:nmatglob
    Aglob = Aglob + coeff_deriv(kk,1)*LHS{kk};
end

  

% RHSderiv{1,1} = ones(ndof,1);


%--------------------------------------------------------------------------
% MDWCAWE Algorithm
%--------------------------------------------------------------------------


Vtilde(:,1) = Aglob\RHSderiv{1}; %Aglob\RHS_0  RHS_0 : 0 order derivative of RHS
U(1,1) = norm(Vtilde(:,1));
V(:,1) = Vtilde(:,1)/U(1,1);    % normalization

for nn=2:nvecfreq
    term1 = 0;
    term2 = 0;
    for mm=1:nn-1
        PU1 = Pu1(nn,mm,U);
        term1 = term1 + RHSderiv{mm+1}*PU1(1,nn-mm);
    end% mm
    for mm=2:nn-1
        PU2 = Pu2(nn,mm,U);
        A_m = sparse(ndof,ndof);
        for ll=1:nmatglob
            A_m = A_m + coeff_deriv(ll,mm+1)*LHS{ll};
        end
        term2 = term2 + A_m*V(:,1:nn-mm)*PU2(:,nn-mm);
    end% mm
    A_1 = sparse(ndof,ndof);
    for ll=1:nmatglob
        A_1 = A_1 + coeff_deriv(ll,2)*LHS{ll};
    end
    Vtilde(:,nn) = Aglob\(term1-term2-A_1*V(:,nn-1));
    for alpha=1:nn-1
        U(alpha,nn) = V(:,alpha)'*Vtilde(:,nn);
        Vtilde(:,nn) = Vtilde(:,nn) - U(alpha,nn)*V(:,alpha);
    end %alpha
    U(nn,nn) = norm(Vtilde(:,nn));
    V(:,nn) = Vtilde(:,nn)/U(nn,nn);
end %nn

Wtrans = V;
Ucoeff = U;

timing.WCAWE = toc;
outputdisplay = sprintf('[WCAWE] CPUtime for building of WCAWE basis (%d vectors): %.4f s',size(Wtrans,2),timing.WCAWE);
disp(outputdisplay);


function Pu = Pu1(nn,mm,U)
    Pu = eye(nn-mm);
    for t=1:mm
        Pu = Pu/U(t:nn-mm+t-1,t:nn-mm+t-1);
    end
end

function Pu = Pu2(nn,mm,U)
    Pu = eye(nn-mm);
    for t=2:mm
        Pu = Pu/U(t:nn-mm+t-1,t:nn-mm+t-1);
    end
end


end





