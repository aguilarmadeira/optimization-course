% EX03_CAPITULO  Corre, por ordem, os 9 exemplos MATLAB do capítulo 3 (Otimização multidimensional sem restrições).
%
%   É o ficheiro que o link «Open in MATLAB Online» abre. Carregue em Run
%   (e, se o MATLAB perguntar, escolha «Change Folder»): cada exemplo imprime
%   os resultados e termina com «confere com os slides: sim».
%
%   Toolboxes: ex03_3_comparison.m usa o fminunc (Optimization Toolbox) só
%   numa parte informativa, protegida por try/catch: sem a toolbox deixa um
%   aviso e o exemplo corre até ao fim.
%
%   Otimização — capítulo 3.  J. F. A. Madeira — Licença MIT.

aqui = fileparts(mfilename('fullpath'));
exemplos = {'03_2_1_random_search/ex03_2_1_random_search.m', ...
            '03_2_2_grid_search/ex03_2_2_grid_search.m', ...
            '03_2_3_nelder_mead/ex03_2_3_nelder_mead.m', ...
            '03_2_4_box/ex03_2_4_box.m', ...
            '03_2_5_hooke_jeeves/ex03_2_5_hooke_jeeves.m', ...
            '03_3_1_steepest_descent/ex03_3_1_steepest_descent.m', ...
            '03_3_2_newton/ex03_3_2_newton.m', ...
            '03_3_3_conjugate_gradients/ex03_3_3_conjugate_gradients.m', ...
            '03_3_comparison/ex03_3_comparison.m'};
for iex = 1:numel(exemplos)
  fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 72), exemplos{iex}, repmat('=', 1, 72));
  try
    run(fullfile(aqui, exemplos{iex}));
  catch erro_ex
    fprintf('\n*** %s não correu: %s\n', exemplos{iex}, erro_ex.message);
  end
end
fprintf('\n%s\nFim do capítulo 3: %d exemplos.\n', repmat('=', 1, 72), numel(exemplos));
