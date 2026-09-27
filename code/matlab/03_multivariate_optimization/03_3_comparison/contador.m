function varargout = contador(acao, x)
%CONTADOR  Conta as avaliações de f, g e H da função fg (Rosenbrock) numa corrida.
%
%   contador('inicio', ftol) põe os contadores a zero; ftol = alvo
%   v = contador('f', x)     f(x)       (n_f = n_f + 1)
%   g = contador('g', x)     gradiente  (n_g = n_g + 1)
%   H = contador('H', x)     Hessiana   (n_H = n_H + 1)
%   C = contador('estado')   estrutura com .nf, .ng, .nH (totais), .fmin e
%                            .hit = [n_f n_g n_H] na primeira avaliação com
%                            f < ftol ([] se nunca)
%
%   Usa-se passando aos métodos F = @(x) contador('f', x),
%   G = @(x) contador('g', x) e H = @(x) contador('H', x) (ver
%   ex03_3_comparison).
%   O estado fica numa variável global (CONTADOR), porque uma função anónima
%   do MATLAB/Octave não pode alterar variáveis.
%
%   Otimização — comparação do cap. 3.  J. F. A. Madeira — Licença MIT.

global CONTADOR
switch acao
  case 'inicio'
    CONTADOR = struct('nf', 0, 'ng', 0, 'nH', 0, 'ftol', x, 'fmin', Inf, 'hit', []);
  case 'f'
    v = fg(x);                                  % só f (nargout = 1)
    CONTADOR.nf = CONTADOR.nf + 1;
    CONTADOR.fmin = min(CONTADOR.fmin, v);
    if v < CONTADOR.ftol && isempty(CONTADOR.hit)   % primeira avaliação abaixo do alvo
      CONTADOR.hit = [CONTADOR.nf, CONTADOR.ng, CONTADOR.nH];
    end
    varargout{1} = v;
  case 'g'
    [~, g] = fg(x);                             % f e gradiente (nargout = 2)
    CONTADOR.ng = CONTADOR.ng + 1;
    varargout{1} = g;
  case 'H'
    [~, ~, H] = fg(x);                          % f, gradiente e Hessiana
    CONTADOR.nH = CONTADOR.nH + 1;
    varargout{1} = H;
  case 'estado'
    varargout{1} = CONTADOR;
  otherwise
    error('contador: ação desconhecida ''%s''', acao);
end
end
