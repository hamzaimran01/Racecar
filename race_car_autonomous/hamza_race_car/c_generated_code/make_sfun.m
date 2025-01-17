%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;
%

SOURCES = [ ...
            'Spatialbycicle_model_model/Spatialbycicle_model_expl_ode_fun.c ', ...
            'Spatialbycicle_model_model/Spatialbycicle_model_expl_vde_forw.c ',...
            'Spatialbycicle_model_constraints/Spatialbycicle_model_constr_h_fun.c ', ...
            'Spatialbycicle_model_constraints/Spatialbycicle_model_constr_h_fun_jac_uxt_hess.c ', ...
            'Spatialbycicle_model_constraints/Spatialbycicle_model_constr_h_fun_jac_uxt_zt.c ', ...
            'acados_solver_sfunction_Spatialbycicle_model.c ', ...
            'acados_solver_Spatialbycicle_model.c '
          ];

INC_PATH = '/home/hamza/Desktop/THESIS/venv/env/acados/include';

INCS = [ ' -I', fullfile(INC_PATH, 'blasfeo', 'include'), ...
         ' -I', fullfile(INC_PATH, 'hpipm', 'include'), ...
        ' -I', INC_PATH, ' -I', fullfile(INC_PATH, 'acados'), ' '];



CFLAGS  = ' -O';



LIB_PATH = '/home/hamza/Desktop/THESIS/venv/env/acados/lib';

LIBS = '-lacados -lhpipm -lblasfeo';



eval( [ 'mex -v -output  acados_solver_sfunction_Spatialbycicle_model ', ...
    CFLAGS, INCS, ' ', SOURCES, ' -L', LIB_PATH, ' ', LIBS ]);

fprintf( [ '\n\nSuccessfully created sfunction:\nacados_solver_sfunction_Spatialbycicle_model', '.', ...
    eval('mexext')] );


%% print note on usage of s-function
fprintf('\n\nNote: Usage of Sfunction is as follows:\n')
input_note = 'Inputs are:\n1) x0, initial state, size [6]\n ';
i_in = 2;
input_note = strcat(input_note, num2str(i_in), ') y_ref - concatenated for intermediate stages,',...
                    ' size [400]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') y_ref_e, size [6]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') lbx, size [1]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') ubx, size [1]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') lbu, size [2]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') ubu, size [2]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') lh, size [5]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') uh, size [5]\n ');
i_in = i_in + 1;

fprintf(input_note)

disp(' ')

output_note = strcat('Outputs are:\n', ...
                '1) u0 - optimal input, size [2]\n',...
                '2) acados solver status (0 = SUCCESS)\n3) KKT residual\n',...
                '4) first state \n5) CPU time\n6) sqp iter\n');

fprintf(output_note)
