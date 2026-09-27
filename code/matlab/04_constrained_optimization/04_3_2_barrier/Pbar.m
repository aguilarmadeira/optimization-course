function v = Pbar(f, g, x, R)
%PBAR  Função de barreira logarítmica: P(x,R) = f(x) - R sum_j ln g_j(x).
%
%   v = Pbar(f, g, x, R)   Inf se algum g_j(x) <= 0 (sem avaliar f);
%                          senão f(x) - R*sum(log(g(x)))
%   Pbar('reset')          põe a zero o contador de chamadas fora de X
%   n = Pbar('fora')       número de chamadas com x fora de X desde o
%                          último 'reset' (essas chamadas avaliam g, não f)
%
%   O contador (variável persistente) serve para BarrierLog separar as
%   chamadas a P que avaliam f das que caem fora de X: funcCount (ou
%   info.nfev do minimizador interno) conta todas as chamadas a P.
%
%   Otimização — deck 4.3.2 (Método da função de penalização interior).
%   Implementação didática.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

persistent nfora
if isempty(nfora), nfora = 0; end

if nargin == 1                              % comandos do contador
  switch f
    case 'reset', nfora = 0;  v = 0;
    case 'fora',  v = nfora;
    otherwise, error('Pbar: comando desconhecido.');
  end
  return
end

gx = g(x);
if any(gx <= 0)                             % nunca avaliar log de um número <= 0
  nfora = nfora + 1;
  v = Inf;
else
  v = f(x) - R*sum(log(gx));
end
end
