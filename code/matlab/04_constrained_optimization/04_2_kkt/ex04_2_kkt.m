% EX04_2_KKT  Reproduz os exemplos do deck 4.2 (condições KKT).
%
%   (a) Resolve os exemplos dos slides pelas condições KKT, com as contas do
%       slide, e verifica-as numericamente com KKTCheck (estacionariedade,
%       admissibilidade, complementaridade u_j g_j = 0 e sinal u_j >= 0).
%       Convenção da UC: g_j(x) >= 0, h_l(x) = 0, L = f - u'g - lambda'h.
%         1. De 4.1 para 4.2 (círculo largo / apertado) e Geometria (2) (u < 0).
%         2. Exemplo 1: os 2^J = 4 casos da receita; sobrevive o caso 3.
%         3. «Cuidado: ativa não implica u > 0»: min x^2 s.a. x >= 0.
%         4. Exemplo 2: verificar o candidato (2,1): u1 = 1/3, u2 = 2/3.
%         5. «KKT dá candidatos»: min -x^2 s.a. x + 1 >= 0, 2 - x >= 0.
%         6. O que medem os multiplicadores: exemplo 1 com b = 3.1 e a produção
%            do 1.1 (preços-sombra 10, 20, 0).
%         7. Quando KKT falha: LICQ.
%         8. Exercícios 1, 2 e 4 (as respostas a cinzento).
%   (b) Resolve o exemplo 1, o exemplo 2, os exercícios 2 e 4 com fmincon
%       (MATLAB, Optimization Toolbox) se existir; senão, em GNU Octave, com
%       sqp (Octave base). A produção com linprog (MATLAB) ou glpk (Octave).
%       Sinais dos multiplicadores (deck 4.2, «No computador»):
%         fmincon: c(x) <= 0 com c = -g  =>  u_UC = lambda.ineqnonlin;
%                  L = f + lambda.eqnonlin'*ceq, ceq = h  =>  lambda_UC = -lambda.eqnonlin.
%         sqp (Octave): h(x) = 0, g(x) >= 0 e L = f - lambda'[h; g], a convenção
%                  da UC: o vetor lambda é [lambda_UC; u_UC], sem troca de sinal.
%         linprog: lambda.ineqlin = u (preços-sombra, >= 0).
%         glpk (maximizar, sense = -1): extra.lambda = u (preços-sombra).
%
%   Usa KKTCheck.m (nesta pasta). Os slides usam vírgula decimal; aqui o ponto.
%   Não há função da UC neste deck: é um exemplo. Implementação didática.
%
%   Otimização — deck 4.2.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
% v está a menos de tol do valor exato s?
perto = @(v, s, tol) all(abs(v(:) - s(:)) <= tol);
simnao = {'não', 'sim'};
vec = @(v) sprintf(' %.4f', v);
ok = true;
tolA = 1e-10;       % contas exatas (a)
tolB = 1e-4;        % solver (b)

fprintf('Otimização — deck 4.2: condições de Karush-Kuhn-Tucker\n');
fprintf('Convenção da UC: g_j >= 0, h_l = 0, L = f - u''g - lambda''h\n');
fprintf('\n=============== (a) pelas condições KKT ===============\n');

% ------------------------------------------------ 1. De 4.1 para 4.2
f  = @(x) (x(1) - 1)^2 + (x(2) - 1)^2;
df = @(x) [2*(x(1) - 1); 2*(x(2) - 1)];
dg = @(x) [-2*x(1), -2*x(2)];                  % grad g' para g = R2 - x1^2 - x2^2
fprintf('\n1. min (x1-1)^2 + (x2-1)^2  s.a.  g = R^2 - x1^2 - x2^2 >= 0\n');
% círculo largo, R^2 = 4: x* = (1,1), inativa, u = 0
x = [1; 1];  g = 4 - x'*x;  u = 0;
[c1, r] = KKTCheck(df(x), g, dg(x), u, [], [], [], tolA);
fprintf('   largo (R^2 = 4):    x* = (1, 1), g = %g > 0 (inativa), u* = 0, grad f = 0;  KKT: %s\n', ...
        g, simnao{c1 + 1});
% círculo apertado, R^2 = 1: x* = (1,1)/sqrt(2), ativa, grad f = u grad g
x = [1; 1]/sqrt(2);  g = 1 - x'*x;
u = dg(x)' \ df(x);  uap = u;
[c2, r] = KKTCheck(df(x), g, dg(x), u, [], [], [], tolA);
fprintf('   apertado (R^2 = 1): x* = (%.4f, %.4f), g = %.1e (ativa), u* = %.4f = sqrt(2) - 1;  KKT: %s\n', ...
        x, g, u, simnao{c2 + 1});
% Geometria (2): (sqrt(2), sqrt(2)) no círculo largo tem u = 1/sqrt(2) - 1 < 0
x = [1; 1]*sqrt(2);  g = 4 - x'*x;
u = dg(x)' \ df(x);
[c3, r] = KKTCheck(df(x), g, dg(x), u, [], [], [], tolA);
fprintf('   Geometria (2): x = (sqrt 2, sqrt 2) no largo: u = %.4f = 1/sqrt(2) - 1 < 0;  KKT: %s (falha o sinal: %.4f)\n', ...
        u, simnao{c3 + 1}, r.sinal);
c = [c1, c2, perto(uap, sqrt(2) - 1, tolA), perto(u, 1/sqrt(2) - 1, tolA), ~c3];
fprintf('   largo KKT: %s | apertado KKT: %s | u* = sqrt2 - 1: %s | u = 1/sqrt2 - 1: %s | (sqrt2, sqrt2) não é KKT: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 2. Exemplo 1: 2^J casos
f1  = @(x) (x(1) - 2)^2 + 2*(x(2) - 1)^2;
df1 = @(x) [2*(x(1) - 2); 4*(x(2) - 1)];
g1  = @(x) [3 - x(1) - 4*x(2); x(1) - x(2)];
Jg1 = [-1 -4; 1 -1];
fprintf('\n2. Exemplo 1: min (x1-2)^2 + 2(x2-1)^2  s.a.  g1 = 3 - x1 - 4x2 >= 0,  g2 = x1 - x2 >= 0\n');
fprintf('   2(x1-2) + u1 - u2 = 0,  4(x2-1) + 4u1 + u2 = 0,  u1 g1 = 0,  u2 g2 = 0,  u >= 0\n');
casos = [0 0; 0 1; 1 0; 1 1];                  % 1 = restrição ativa (g_j = 0); 0 = u_j = 0
hip = {'u1 = u2 = 0     ', 'u1 = 0, g2 = 0  ', 'g1 = 0, u2 = 0  ', 'g1 = g2 = 0     '};
Xs = [2 1; 4/3 4/3; 5/3 1/3; 3/5 3/5];         % slide
Us = [0 0; 0 -4/3; 2/3 0; 22/25 -48/25];
fprintf('   %-5s %-17s %-19s %-19s %-16s %s\n', 'caso', 'hipótese', 'x', 'u', 'g', 'verificação');
c = [];
for k = 1:4
  M = [2 0 1 -1; 0 4 4 1];  rhs = [4; 4];      % estacionariedade em (x1, x2, u1, u2)
  if casos(k, 1), M = [M; 1 4 0 0]; rhs = [rhs; 3]; else, M = [M; 0 0 1 0]; rhs = [rhs; 0]; end
  if casos(k, 2), M = [M; 1 -1 0 0]; rhs = [rhs; 0]; else, M = [M; 0 0 0 1]; rhs = [rhs; 0]; end
  z = M \ rhs;  x = z(1:2);  u = z(3:4);  g = g1(x);
  ver = '';                                    % todas as falhas do caso
  for j = find(g' < -tolA), ver = [ver, sprintf('g%d < 0 (inadmissível); ', j)]; end
  for j = find(u' < -tolA), ver = [ver, sprintf('u%d < 0; ', j)]; end
  if isempty(ver), ver = 'ponto KKT  V'; else, ver = [ver(1:end-2), '  x']; end
  fprintf('   %-5d %s (%7.4f,%7.4f) (%7.4f,%7.4f) (%6.3f,%6.3f) %s\n', k, hip{k}, x, u, g, ver);
  c = [c, perto(x, Xs(k, :)', tolA), perto(u, Us(k, :)', tolA)];
  if k == 3, x3 = x; u3 = u; end
end
[cK, r] = KKTCheck(df1(x3), g1(x3), Jg1, u3, [], [], [], tolA);
gx3 = g1(x3);
fprintf('   Sobrevive o caso 3: x* = (5/3, 1/3), f* = %.4f; ativas I = {%s}; u1 g1 = %.1e, u2 g2 = %.1e\n', ...
        f1(x3), num2str(r.ativas), u3.*gx3);
fprintf('   Geometria (2): grad f(x*) = (%.4f, %.4f) = u1 grad g1 = (2/3)(-1,-4)\n', df1(x3));
c = [all(c), cK, perto(f1(x3), 1, tolA), perto(gx3(2), 4/3, tolA), perto(df1(x3), [-2/3; -8/3], tolA)];
fprintf('   tabela dos 4 casos: %s | caso 3 é KKT: %s | f* = 1: %s | g2 = 4/3: %s | grad f = (2/3)(-1,-4): %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 3. Ativa com u = 0
fprintf('\n3. min x^2  s.a.  g = x >= 0:  x* = 0, g(x*) = 0 (ativa), 2x - u = 0 => u* = 0\n');
[cK, r] = KKTCheck(0, 0, 1, 0, [], [], [], tolA);
fprintf('   KKT: %s; ativas I = {%s} com u = 0 (ativa não implica u > 0)\n', simnao{cK + 1}, num2str(r.ativas));
ok = ok && cK && isequal(r.ativas, 1);

% ------------------------------------------------ 4. Exemplo 2
f2  = @(x) (x(1) - 3)^2 + (x(2) - 2)^2;
df2 = @(x) [2*(x(1) - 3); 2*(x(2) - 2)];
g2  = @(x) [5 - x(1)^2 - x(2)^2; 4 - x(1) - 2*x(2); x(1); x(2)];
Jg2 = @(x) [-2*x(1) -2*x(2); -1 -2; 1 0; 0 1];
fprintf('\n4. Exemplo 2: min (x1-3)^2 + (x2-2)^2  s.a.  5 - x1^2 - x2^2 >= 0, 4 - x1 - 2x2 >= 0, x1 >= 0, x2 >= 0\n');
x = [2; 1];  g = g2(x);  J = Jg2(x);
I = find(abs(g) <= tolA);                       % ativas
u = zeros(4, 1);
u(I) = J(I, :)' \ df2(x);                       % grad f = sum_{j em I} u_j grad g_j
[cK, r] = KKTCheck(df2(x), g, J, u, [], [], [], tolA);
fprintf('   x = (2, 1): g = (%s), ativas I = {%s}; grad f = (%g, %g), grad g1 = (%g, %g), grad g2 = (%g, %g)\n', ...
        num2str(g', '%g '), num2str(I'), df2(x), J(1, :), J(2, :));
fprintf('   u = (%s)  (slide: 1/3, 2/3, 0, 0);  LICQ: característica de [grad g1 grad g2] = %d;  f = %g;  KKT: %s\n', ...
        vec(u), rank(J(I, :)), f2(x), simnao{cK + 1});
c = [perto(u, [1/3; 2/3; 0; 0], tolA), rank(J(I, :)) == 2, perto(f2(x), 2, tolA), cK];
fprintf('   u1 = 1/3, u2 = 2/3: %s | LICQ: %s | f = 2: %s | ponto KKT: %s\n', simnao{c + 1});
ok = ok && all(c);
x2s = x;  u2s = u;

% ------------------------------------------------ 5. KKT dá candidatos
fprintf('\n5. min -x^2  s.a.  g1 = x + 1 >= 0,  g2 = 2 - x >= 0  (os 4 casos; «ambas ativas» é impossível)\n');
fm = @(x) -x^2;  dfm = @(x) -2*x;  gm = @(x) [x + 1; 2 - x];  Jgm = [1; -1];
cand = [0, -1, 2];  ucand = [0 0; 2 0; 0 4];  nome = {'máximo', 'mínimo local', 'mínimo global'};
c = [];
for i = 1:3
  x = cand(i);  u = ucand(i, :)';
  [cK, r] = KKTCheck(dfm(x), gm(x), Jgm, u, [], [], [], tolA);
  fprintf('   x = %2g: u = (%g, %g), f = %2g, f'''' = -2;  KKT: %s -> %s\n', x, u, fm(x) + 0, simnao{cK + 1}, nome{i});
  c = [c, cK];
end
[~, imin] = min(arrayfun(fm, cand));
c = [c, cand(imin) == 2];
fprintf('   os três são pontos KKT: %s | o menor f (global) é x = 2: %s\n', simnao{all(c(1:3)) + 1}, simnao{c(4) + 1});
ok = ok && all(c);

% ------------------------------------------------ 6. O que medem os multiplicadores
fprintf('\n6. Sensibilidade no exemplo 1: g1 = b - x1 - 4x2 >= 0, b = 3 -> 3.1 (caso 3: g1 ativa, u2 = 0)\n');
sel = @(v, i) v(i);                            % (MATLAB não indexa o resultado de uma chamada)
fb = @(b) f1(sel([2 0 1 -1; 0 4 4 1; 1 4 0 0; 0 0 0 1] \ [4; 4; b; 0], 1:2));
f31 = fb(3.1);
fprintf('   f*(3) = %.4f, f*(3.1) = %.5f (exato: (6 - b)^2/9); aproximação 1 - u1*0.1 = %.4f\n', ...
        fb(3), f31, 1 - 2/3*0.1);
c = [perto(fb(3), 1, tolA), confere(f31, 0.934, 3), confere(1 - 2/3*0.1, 0.933, 3), ...
     perto(f31, (6 - 3.1)^2/9, tolA)];
fprintf('   f*(3) = 1: %s | f*(3.1) = 0.934: %s | aproximação 0.933: %s | (6-b)^2/9: %s\n', simnao{c + 1});
ok = ok && all(c);

% produção do 1.1: max 40x1 + 30x2  <=>  min -(40x1 + 30x2)
fprintf('\n   Produção do 1.1: max 40x1 + 30x2 s.a. 2x1 + x2 <= 100 (máquina), x1 + x2 <= 80 (mão de obra),\n');
fprintf('   x1 <= 40 (procura), x >= 0.  Em g_j >= 0: g = b - A x >= 0 e x >= 0;  f = -(40x1 + 30x2)\n');
Ap = [2 1; 1 1; 1 0];  bp = [100; 80; 40];  cp = [40; 30];
x = [20; 60];
gp = [bp - Ap*x; x];  Jgp = [-Ap; eye(2)];
I = find(abs(gp) <= tolA);
u = zeros(5, 1);
u(I) = Jgp(I, :)' \ (-cp);                      % grad f = -c = sum u_j grad g_j
[cK, r] = KKTCheck(-cp, gp, Jgp, u, [], [], [], tolA);
fprintf('   x* = (20, 60), L* = %g;  ativas I = {%s} (máquina, mão de obra): 2u1 + u2 = 40, u1 + u2 = 30\n', ...
        cp'*x, num2str(I'));
fprintf('   u = (%s)  -> preços-sombra 10 EUR/h, 20 EUR/h, 0;  KKT: %s\n', vec(u), simnao{cK + 1});
c = [perto(cp'*x, 2600, tolA), perto(u, [10; 20; 0; 0; 0], tolA), cK];
fprintf('   L* = 2600: %s | u = (10, 20, 0): %s | ponto KKT: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 7. Quando KKT falha: LICQ
fprintf('\n7. min -x1  s.a.  g1 = (1 - x1)^3 - x2 >= 0,  g2 = x2 >= 0;  x* = (1, 0)\n');
x = [1; 0];  gf = [-1; 0];
J = [-3*(1 - x(1))^2 + 0, -1; 0, 1];         % (+ 0: evita imprimir -0)
u = pinv(J')*gf;                                % melhor u (mínimos quadrados: J' é singular)
resmin = norm(gf - J'*u);
fprintf('   grad g1 = (%g, %g), grad g2 = (%g, %g): característica %d (dependentes)\n', J(1, :), J(2, :), rank(J));
fprintf('   min_u |grad f - J''u| = %.4f: a componente x1 de grad f = (-1, 0) não se obtém -> não há u\n', resmin);
c = [rank(J) == 1, perto(resmin, 1, tolA)];
fprintf('   gradientes dependentes: %s | sem multiplicadores (resíduo 1): %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 8. Exercícios
fprintf('\n8. Exercícios\n');
% Ex. 1: min -ln(x1+1) - x2 s.a. 3 - 2x1 - x2 >= 0, x1 >= 0, x2 >= 0
fe1 = @(x) -log(x(1) + 1) - x(2);  dfe1 = @(x) [-1/(x(1) + 1); -1];
ge1 = @(x) [3 - 2*x(1) - x(2); x(1); x(2)];  Jge1 = [-2 -1; 1 0; 0 1];
x = [0; 3];  g = ge1(x);  I = find(abs(g) <= tolA);
u = zeros(3, 1);  u(I) = Jge1(I, :)' \ dfe1(x);
[cK, r] = KKTCheck(dfe1(x), g, Jge1, u, [], [], [], tolA);
fprintf('   Ex. 1: x* = (0, 3), ativas I = {%s}, u = (%s), f* = %.4f;  KKT: %s\n', num2str(I'), vec(u), fe1(x), simnao{cK + 1});
c = [isequal(I', [1 2]), cK];
% Ex. 2: min x1^2 + x2^2 s.a. 5 - x1^2 - x2^2 >= 0, x1 + 2x2 = 4, x >= 0
fe2 = @(x) x(1)^2 + x(2)^2;  dfe2 = @(x) 2*x;
ge2 = @(x) [5 - x(1)^2 - x(2)^2; x(1); x(2)];  Jge2 = @(x) [-2*x(1) -2*x(2); 1 0; 0 1];
he2 = @(x) x(1) + 2*x(2) - 4;  Jhe2 = [1 2];
x = [4/5; 8/5];  g = ge2(x);
lam = Jhe2' \ dfe2(x);                          % só a igualdade é ativa: u = 0
[cK, r] = KKTCheck(dfe2(x), g, Jge2(x), zeros(3, 1), he2(x), Jhe2, lam, tolA);
fprintf('   Ex. 2: x* = (4/5, 8/5), g = (%s) > 0 (inativas, u = 0), lambda = %.4f (slide: 8/5);  KKT: %s\n', ...
        vec(g), lam, simnao{cK + 1});
c = [c, isempty(r.ativas), perto(lam, 8/5, tolA), cK];
% Ex. 4: a lata com h <= 6 (g = 6 - h >= 0), z = (r, h)
V = 330;
A  = @(z) 2*pi*z(1)^2 + 2*pi*z(1)*z(2);  dA = @(z) [4*pi*z(1) + 2*pi*z(2); 2*pi*z(1)];
hv = @(z) pi*z(1)^2*z(2) - V;            Jhv = @(z) [2*pi*z(1)*z(2), pi*z(1)^2];
gh = @(z) 6 - z(2);                      Jgh = [0 -1];
hh = 6;  rr = sqrt(V/(6*pi));  z = [rr; hh];
lam = (2*rr + hh)/(rr*hh);                      % de dL/dr = 0
u = lam*pi*rr^2 - 2*pi*rr;                      % de dL/dh = 0: 2 pi r - lambda pi r^2 + u = 0
[cK, r] = KKTCheck(dA(z), gh(z), Jgh, u, hv(z), Jhv(z), lam, 1e-9);
fprintf('   Ex. 4: h* = 6, r* = %.4f cm, A* = %.4f cm^2, lambda = %.4f, u* = %.4f cm^2/cm;  KKT: %s\n', ...
        rr, A(z), lam, u, simnao{cK + 1});
c = [c, confere(rr, 4.18, 2), confere(A(z), 267.7, 1), confere(u, 5.19, 2), cK];
fprintf('   Ex. 1 ativas {1,2}: %s | KKT: %s | Ex. 2 inativas: %s | lambda = 8/5: %s | KKT: %s | Ex. 4 r* = 4.18: %s | A* = 267.7: %s | u* = 5.19: %s | KKT: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% =================================================== (b) com um solver
fprintf('\n=============== (b) com um solver ===============\n');
temFmincon = exist('fmincon', 'file') == 2;
temSqp = exist('sqp', 'file') == 2;
if temFmincon
  fprintf('Solver: fmincon. c = -g <= 0, ceq = h:  u_UC = lambda.ineqnonlin,  lambda_UC = -lambda.eqnonlin.\n');
  opts = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12, 'TolCon', 1e-12);
elseif temSqp
  fprintf('Solver: sqp (GNU Octave; fmincon não existe). O sqp usa g(x) >= 0, h(x) = 0 e L = f - lambda''[h; g]:\n');
  fprintf('o vetor lambda é [lambda_UC; u_UC], sem troca de sinal.\n');
  fprintf('Com o fmincon (MATLAB): u_UC = lambda.ineqnonlin, lambda_UC = -lambda.eqnonlin (deck 4.2, «No computador»).\n');
else
  fprintf('Nem fmincon nem sqp disponíveis: problemas não lineares não corridos.\n');
end
if temFmincon || temSqp
  % f, g (>= 0), h (= 0, ou []), x0, x* exato, u* exato, lambda* exato, nome
  probs = {f1, g1, [], [0; 0], [5/3; 1/3], [2/3; 0], [], 'Exemplo 1';
           f2, g2, [], [0.5; 0.5], [2; 1], [1/3; 2/3; 0; 0], [], 'Exemplo 2';
           fe2, ge2, he2, [1; 1], [4/5; 8/5], [0; 0; 0], 8/5, 'Exercício 2';
           A, gh, hv, [3; 6], [], [], [], 'Exercício 4 (lata, h <= 6)'};
  for i = 1:size(probs, 1)
    [fo, gc, hc, x0, xs, us, ls, nome] = probs{i, :};
    if temFmincon
      if isempty(hc), hcf = @(x) []; else, hcf = hc; end
      [x, fx, flag, out, lambda] = fmincon(fo, x0, [], [], [], [], [], [], ...
                                           @(x) deal(-gc(x), hcf(x)), opts);
      u = lambda.ineqnonlin;  lam = -lambda.eqnonlin;  nf = out.funcCount;  nit = out.iterations;
    else
      [x, fx, flag, nit, nf, lambda] = sqp(x0, fo, hc, gc);
      K = 0;  if ~isempty(hc), K = numel(hc(x0)); end
      lam = lambda(1:K);  u = lambda(K + 1:end);
    end
    fprintf('   %s: x = (%s), f = %.4f, u = (%s)', nome, vec(x), fx, vec(u));
    if ~isempty(lam), fprintf(', lambda = (%s)', vec(lam)); end
    fprintf('  [flag %d, %d it., %d aval.]\n', flag, nit, nf);
    if isempty(xs)          % exercício 4: comparar com os valores do slide
      c = [confere(x(2), 6, 2), confere(x(1), 4.18, 2), confere(fx, 267.7, 1), confere(u, 5.19, 2)];
      fprintf('      h* = 6: %s | r* = 4.18: %s | A* = 267.7: %s | u* = 5.19: %s\n', simnao{c + 1});
    else
      c = [perto(x, xs, tolB), perto(u, us, tolB)];
      if ~isempty(ls), c = [c, perto(lam, ls, tolB)]; end
      fprintf('      x*: %s | u*: %s', simnao{c(1:2) + 1});
      if ~isempty(ls), fprintf(' | lambda*: %s', simnao{c(3) + 1}); end
      fprintf('\n');
    end
    ok = ok && all(c);
  end
end

% produção (PL): linprog (MATLAB) ou glpk (Octave)
if exist('linprog', 'file') == 2
  [x, fv, flag, out, lambda] = linprog(-cp, Ap, bp, [], [], [0; 0]);
  u = lambda.ineqlin;  L = -fv;
  fprintf('   Produção com linprog (min -L): x = (%s), L = %.4f, lambda.ineqlin = (%s)\n', vec(x), L, vec(u));
elseif exist('glpk', 'file') == 2
  [x, L, err, extra] = glpk(cp, Ap, bp, [0; 0], [], 'UUU', 'CC', -1);
  u = extra.lambda;
  fprintf('   Produção com glpk (Octave; linprog não existe; sense = -1 maximiza L):\n');
  fprintf('      x = (%s), L = %.4f, extra.lambda = (%s) = preços-sombra\n', vec(x), L, vec(u));
else
  u = [];
  fprintf('   Produção: nem linprog nem glpk disponíveis.\n');
end
if ~isempty(u)
  c = [perto(x, [20; 60], tolB), perto(L, 2600, tolB), perto(u, [10; 20; 0], tolB)];
  fprintf('      x* = (20, 60): %s | L* = 2600: %s | preços-sombra (10, 20, 0): %s\n', simnao{c + 1});
  ok = ok && all(c);
end

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
