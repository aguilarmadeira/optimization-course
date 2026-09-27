function [ok, res] = KKTCheck(df, g, Jg, u, h, Jh, lam, tol)
%KKTCHECK  Verifica as condições KKT num ponto, na convenção da UC.
%
%   [ok, res] = KKTCheck(df, g, Jg, u)                 só desigualdades
%   [ok, res] = KKTCheck(df, g, Jg, u, h, Jh, lam)     com igualdades
%   [ok, res] = KKTCheck(df, g, Jg, u, h, Jh, lam, tol)
%
%   Problema: min f(x) s.a. g_j(x) >= 0 (j = 1..J), h_l(x) = 0 (l = 1..K).
%   Lagrangiano da UC: L = f - u'g - lambda'h.  Tudo avaliado no ponto x:
%     df    gradiente de f (n x 1)
%     g     valores g_j(x) (J x 1);  Jg  jacobiana (J x n): linha j = grad g_j'
%     u     multiplicadores das desigualdades (J x 1)
%     h, Jh, lam   o mesmo para as igualdades (K x 1, K x n, K x 1);
%                  [] se não houver igualdades (ou se g = [] e só houver h)
%     tol   tolerância (opcional, 1e-8)
%
%   Saídas
%     ok    true se as quatro condições se cumprem (todas <= tol)
%     res   estrutura com os resíduos
%       .estac   estacionariedade  || df - Jg'u - Jh'lambda ||_inf
%       .admis   admissibilidade   max( max(-g_j, 0), |h_l| )
%       .compl   complementaridade max |u_j g_j|
%       .sinal   sinal             max( -u_j, 0 )   (u_j >= 0; lambda livre)
%       .ativas  índices j com |g_j| <= tol (conjunto ativo I(x))
%       .ok      o mesmo que ok
%
%   Otimização — deck 4.2 (Condições KKT). Usada em ex04_2_kkt.m.
%   Implementação didática: só verifica as condições de 1.ª ordem; não
%   verifica a qualificação das restrições (LICQ) nem classifica o ponto.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 5, h = []; Jh = []; lam = []; end
if nargin < 8 || isempty(tol), tol = 1e-8; end
df = df(:);  n = numel(df);
g = g(:);  u = u(:);  h = h(:);  lam = lam(:);
if isempty(g), Jg = zeros(0, n); u = zeros(0, 1); end
if isempty(h), Jh = zeros(0, n); lam = zeros(0, 1); end

r = df - Jg'*u - Jh'*lam;
res.estac = norm(r, inf);
res.admis = max([0; -g; abs(h)]);
res.compl = max([0; abs(u.*g)]);
res.sinal = max([0; -u]);
res.ativas = find(abs(g) <= tol)';
ok = res.estac <= tol && res.admis <= tol && res.compl <= tol && res.sinal <= tol;
res.ok = ok;
end
