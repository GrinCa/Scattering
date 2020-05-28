function create_cross_derivatives(listLHS,coeff_matrices,coeff_RHS,derivative_orders,variables)

% this function creates a file called 'get_coeff_deriv_matrices.m' which
% enables us to fill the arrays of cross derivatives coeff.

% save array of derivative order for each variable
save('Derivatives/derivative_orders.mat','derivative_orders');

% number of matrices
nmat = length(listLHS);

%number of variables of the given problem
nvariables = length(variables);

% char_variable = 'name_var1,name_var2,name_var3,...' 
char_variable = ''; 

%--------------------------------------------------------------------------
sym_cell = cell(1,nvariables); % cell of symbols created with variables strings
for ii=1:nvariables
    sym_cell{ii}=sym(variables(ii));
    char_variable = strcat(char_variable,variables(ii),',');
    char_out = char_variable(1:length(char_variable)-1);
end
%--------------------------------------------------------------------------

%definition of global variable needed for creating array of index, for
%derivative

% matrix_tuples is a matrix made of all the different tupples created with
% a given set (here the set is derivative_order). For instance, if the set
% is [2,3], matrix_tupple will be :
%[0 0;
% 0 1
% ...
% 2 2
% 2 3]
%Note that it varies from 0 to n, because we want the 0 order derivative 

global matrix_tuples
matrix_tuples = [];

%get_array_index build matrix_tuples. It is a recursive function
for ii=0:derivative_orders(1)
    get_array_index(ii,variables,derivative_orders);
end

%--------------------------------------------------------------------------
% Writing file get_coeff_deriv_matrices.m
%--------------------------------------------------------------------------

filID = fopen('Derivatives/get_coeff_deriv_matrices.m','w');

fprintf(filID,"function [coeff_deriv_fun,RHScoeffderiv] = get_coeff_deriv_matrices(coeff_deriv_fun,RHScoeffderiv,derivative_orders,nmat)");

chartmp = '';
for ii=1:length(derivative_orders)
    chartmp = strcat(chartmp,num2str(derivative_orders(ii)),',');
end

chartmp = strcat('\n\nderiv_orders_pre_calc = [',chartmp(1:length(chartmp)-1),'];\n');

fprintf(filID,chartmp);
fprintf(filID,strcat('nmat_pre_calc=',num2str(nmat),';\n\n'));

fprintf(filID,strcat('if nmat > nmat_pre_calc\n\t',...
                    "disp('[Error] : data mismatch')\n\t",...
                    "return;\n",...
                    "end\n\n"));

fprintf(filID,strcat("for ii=1:length(derivative_orders)\n\t",...
                        "if derivative_orders(ii)>deriv_orders_pre_calc(ii)\n\t\t",...
                            "disp('[Error] : data mismatch')\n\t\t",...
                            "return;\n\t",...
                        "end\n",...
                    "end\n\n"));
%calculate all the cross derivatives of the left hand side
for ii=1:nmat
    for jj=1:size(matrix_tuples,1)
        cross_deriv = coeff_matrices{ii};
        index = num2str(ii);
        for kk=1:size(matrix_tuples,2)
            cross_deriv = simplify(diff(cross_deriv,sym_cell{kk},matrix_tuples(jj,kk)));
            index = strcat(index,',',num2str(matrix_tuples(jj,kk)+1));
        end
        fprintf(filID,strcat("coeff_deriv_fun{",index,"} = @(",char_out,") ",char(cross_deriv),";\n\n"));
        
    end
end

%calculate all the cross derivatives of the right hand side
for jj=1:size(matrix_tuples,1)
    cross_deriv = coeff_RHS;
    index = '1';
    for kk=1:size(matrix_tuples,2)
        cross_deriv = simplify(diff(cross_deriv,sym_cell{kk},matrix_tuples(jj,kk)));
        index = strcat(index,',',num2str(matrix_tuples(jj,kk)+1));
    end
    fprintf(filID,strcat("RHScoeffderiv{",index,"} = @(",char_out,") ",char(cross_deriv),";\n\n"));
end
fprintf(filID,"end");
fclose(filID);
end

%--------------------------------------------------------------------------
% End writing
%--------------------------------------------------------------------------


% here is the recursive function that creates all the tuples given a limit
% range. For a range [l1 l2...ln], we have l1*l2*...*ln tuples. This will
% be useful for cross derivatives.
function get_array_index(array,variables,derivative_orders)
    global matrix_tuples;
    if length(array)>=length(variables)
        matrix_tuples = [matrix_tuples;array];
    else
        for ii=0:derivative_orders(length(array)+1)
            get_array_index([array ii],variables,derivative_orders);
        end
    end
end




