% EX02_CAPITULO  Corre, por ordem, os 6 exemplos MATLAB do capítulo 2 (Otimização unidimensional).
%
%   É o ficheiro que o link «Open in MATLAB Online» abre. Carregue em Run
%   (e, se o MATLAB perguntar, escolha «Change Folder»): cada exemplo imprime
%   os resultados e termina com «confere com os slides: sim».
%
%   Não usa toolboxes: corre no MATLAB, no MATLAB Online e no GNU Octave.
%
%   Otimização — capítulo 2.  J. F. A. Madeira — Licença MIT.

aqui = fileparts(mfilename('fullpath'));
exemplos = {'02_2_golden_section/ex02_2_golden_section.m', ...
            '02_3_powell/ex02_3_powell.m', ...
            '02_4_bisection/ex02_4_bisection.m', ...
            '02_5_newton/ex02_5_newton.m', ...
            '02_5_newton/ex02_5_newton_safeguarded.m', ...
            '02_6_finite_differences/ex02_6_finite_differences.m'};
for iex = 1:numel(exemplos)
  fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 72), exemplos{iex}, repmat('=', 1, 72));
  try
    run(fullfile(aqui, exemplos{iex}));
  catch erro_ex
    fprintf('\n*** %s não correu: %s\n', exemplos{iex}, erro_ex.message);
  end
end
fprintf('\n%s\nFim do capítulo 2: %d exemplos.\n', repmat('=', 1, 72), numel(exemplos));
