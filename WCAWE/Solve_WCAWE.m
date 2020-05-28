function [Pout_WCAWE] = Solve_WCAWE(LHS,coeff_deriv_fun,RHS,Wtrans,freq)


nmat_glob = length(LHS);
ndof = size(LHS{1},1);

%--------------------------------------------------------------------------
% Project matrices on WCAWE basis and setup RHS
%--------------------------------------------------------------------------
for ii = 1:nmat_glob
   LHS{ii} = sparse(Wtrans'*LHS{ii}*Wtrans);
end %ii

%--------------------------------------------------------------------------
% Setup projected RHS
%--------------------------------------------------------------------------
RHS = sparse(Wtrans'*RHS);

%--------------------------------------------------------------------------
% Frequency/Fi loops calculation
%--------------------------------------------------------------------------

nfreq = length(freq);

Pout_WCAWE = zeros(ndof,nfreq);


for ii=1:nfreq
      Aglob_red = sparse(size(LHS{1},1),size(LHS{1},2));
      for kk = 1:nmat_glob
         Aglob_red = Aglob_red + coeff_deriv_fun{kk,1}(freq(ii))*LHS{kk};
      end %kk
      solp = Aglob_red\RHS;
      resp_P = Wtrans*solp;
      Pout_WCAWE(:,ii) = resp_P;
end % ii


end

