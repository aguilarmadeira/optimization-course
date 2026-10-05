% EX04_CAPITULO  Corre, por ordem, os 5 exemplos MATLAB do capítulo 4 (Otimização com restrições).
%
%   É o ficheiro que o link «Open in MATLAB Online» abre. Carregue em Run
%   (e, se o MATLAB perguntar, escolha «Change Folder»): cada exemplo imprime
%   os resultados e termina com «confere com os slides: sim».
%
%   Toolboxes: ex04_1_lagrange.m e ex04_2_kkt.m usam o fmincon (e o 4.2
%   também o linprog), da Optimization Toolbox; no GNU Octave usam o sqp e o
%   glpk.
%   No MATLAB Online, esses exemplos só correm se a licença incluir a
%   toolbox; se um exemplo falhar, o script diz qual e continua.
%
%   Otimização — capítulo 4.  J. F. A. Madeira — Licença MIT.

aqui = fileparts(mfilename('fullpath'));
exemplos = {'04_1_lagrange/ex04_1_lagrange.m', ...
            '04_2_kkt/ex04_2_kkt.m', ...
            '04_3_1_exterior_penalty/ex04_3_1_exterior_penalty.m', ...
            '04_3_2_barrier/ex04_3_2_barrier.m', ...
            '04_3_3_constraint_handling/ex04_3_3_constraint_handling.m'};
for iex = 1:numel(exemplos)
  fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 72), exemplos{iex}, repmat('=', 1, 72));
  try
    run(fullfile(aqui, exemplos{iex}));
  catch erro_ex
    fprintf('\n*** %s não correu: %s\n', exemplos{iex}, erro_ex.message);
  end
end
fprintf('\n%s\nFim do capítulo 4: %d exemplos.\n', repmat('=', 1, 72), numel(exemplos));
