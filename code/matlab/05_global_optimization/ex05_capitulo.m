% EX05_CAPITULO  Corre, por ordem, os 2 exemplos MATLAB do capítulo 5 (Otimização global).
%
%   É o ficheiro que o link «Open in MATLAB Online» abre. Carregue em Run
%   (e, se o MATLAB perguntar, escolha «Change Folder»): cada exemplo imprime
%   os resultados e termina com «confere com os slides: sim».
%
%   Não usa toolboxes: corre no MATLAB, no MATLAB Online e no GNU Octave.
%
%   Otimização — capítulo 5.  J. F. A. Madeira — Licença MIT.

aqui = fileparts(mfilename('fullpath'));
exemplos = {'05_2_simulated_annealing/ex05_2_simulated_annealing.m', ...
            '05_3_genetic_algorithms/ex05_3_genetic_algorithms.m'};
for iex = 1:numel(exemplos)
  fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 72), exemplos{iex}, repmat('=', 1, 72));
  try
    run(fullfile(aqui, exemplos{iex}));
  catch erro_ex
    fprintf('\n*** %s não correu: %s\n', exemplos{iex}, erro_ex.message);
  end
end
fprintf('\n%s\nFim do capítulo 5: %d exemplos.\n', repmat('=', 1, 72), numel(exemplos));
