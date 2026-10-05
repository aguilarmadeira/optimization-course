% EX06_CAPITULO  Corre, por ordem, os 4 exemplos MATLAB do capítulo 6 (Otimização multiobjetivo).
%
%   É o ficheiro que o link «Open in MATLAB Online» abre. Carregue em Run
%   (e, se o MATLAB perguntar, escolha «Change Folder»): cada exemplo imprime
%   os resultados e termina com «confere com os slides: sim» (exceto
%   ex06_6_gamultiobj.m, uma variante com o gamultiobj, sem verificação).
%
%   Toolboxes: ex06_3_aggregation.m e ex06_4_epsilon_constraint.m usam o
%   fmincon (Optimization Toolbox; no GNU Octave, o sqp);
%   ex06_6_gamultiobj.m usa o gamultiobj (Global Optimization Toolbox) e,
%   sem ela, só avisa.
%   No MATLAB Online, esses exemplos só correm se a licença incluir a
%   toolbox; se um exemplo falhar, o script diz qual e continua.
%
%   Otimização — capítulo 6.  J. F. A. Madeira — Licença MIT.

aqui = fileparts(mfilename('fullpath'));
exemplos = {'06_3_aggregation/ex06_3_aggregation.m', ...
            '06_4_epsilon_constraint/ex06_4_epsilon_constraint.m', ...
            '06_6_nsga2/ex06_6_nsga2.m', ...
            '06_6_nsga2/ex06_6_gamultiobj.m'};
for iex = 1:numel(exemplos)
  fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 72), exemplos{iex}, repmat('=', 1, 72));
  try
    run(fullfile(aqui, exemplos{iex}));
  catch erro_ex
    fprintf('\n*** %s não correu: %s\n', exemplos{iex}, erro_ex.message);
  end
end
fprintf('\n%s\nFim do capítulo 6: %d exemplos.\n', repmat('=', 1, 72), numel(exemplos));
