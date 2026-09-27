% EX04_1_LAGRANGE  Reproduz os exemplos do deck 4.1 (multiplicadores de Lagrange).
%
%   (a) Resolve os exemplos dos slides pelas condições de Lagrange, com as
%       contas do slide, e verifica numericamente as condições:
%       grad f(x*) - sum_l lambda_l grad h_l(x*) = 0  e  h_l(x*) = 0.
%       Convenção da UC:  L(x, lambda) = f(x) - sum_l lambda_l h_l(x).
%         1. O exemplo resolvido: min 2x1^2 + x2^2 s.a. x1 + x2 = 1
%            (e «Geometria (1)» e «O que significa lambda?»).
%         2. Candidatos: min x1 + x2 s.a. x1^2 + x2^2 = 2.
%         3. Exemplo 2, a lata: min 2 pi r^2 + 2 pi r h s.a. pi r^2 h = V.
%         4. Várias restrições: min x1^2+x2^2+x3^2 s.a. x1+x2+x3 = 3, x1-x2 = 1.
%         5. Quando a condição de Lagrange falha (gradientes dependentes).
%   (b) Resolve 1, 3 e 4 com um solver: fmincon (MATLAB, Optimization
%       Toolbox) se existir; senão, em GNU Octave, sqp (Octave base).
%       Sinais dos multiplicadores:
%         fmincon: L = f + lambda.eqnonlin' * ceq, com ceq = h
%                  =>  lambda_UC = -lambda.eqnonlin   (deck 4.2, «No computador»);
%         sqp (Octave): usa a convenção da UC (h(x) = 0, L = f - lambda' h)
%                  =>  lambda_UC = lambda   (sem troca de sinal).
%
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%   Não há função da UC neste deck: é um exemplo. Implementação didática.
%
%   Otimização — deck 4.1.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
% v está a menos de tol do valor exato s?
perto = @(v, s, tol) all(abs(v(:) - s(:)) <= tol);
simnao = {'não', 'sim'};
ok = true;
tolA = 1e-10;       % contas exatas (a)
tolB = 1e-5;        % solver (b)

% Problema quadrático com igualdades lineares: min 1/2 x'Qx + c'x s.a. Ax = b.
% grad_x L = Qx + c - A'lambda = 0, Ax = b  -> sistema linear em (x, lambda).
lagrange_lin = @(Q, c, A, b) [Q, -A'; A, zeros(size(A, 1))] \ [-c; b];

fprintf('Otimização — deck 4.1: multiplicadores de Lagrange\n');
fprintf('Convenção da UC: L = f - sum lambda_l h_l\n');
fprintf('\n=============== (a) pelas condições de Lagrange ===============\n');

% ------------------------------------------------ 1. O exemplo resolvido
f1  = @(x) 2*x(1)^2 + x(2)^2;
df1 = @(x) [4*x(1); 2*x(2)];
h1  = @(x) x(1) + x(2) - 1;
dh1 = [1; 1];
fprintf('\n1. min 2x1^2 + x2^2  s.a.  h = x1 + x2 - 1 = 0\n');
% Geometria (1): em (0.8, 0.2) a curva de nível atravessa a reta
xg = [0.8; 0.2];  d = [1; -1];
fprintf('   Geometria (1): em x = (0.8, 0.2), grad f = (%.1f, %.1f), grad f''d = %.1f (d = (1,-1))\n', ...
        df1(xg), df1(xg)'*d);
c = [confere(df1(xg), [3.2; 0.4], 1), confere(df1(xg)'*d, 2.8, 1)];
% 4x1 - lambda = 0, 2x2 - lambda = 0, x1 + x2 = 1
z = lagrange_lin(diag([4 2]), [0; 0], [1 1], 1);
x = z(1:2);  lam = z(3);
res = norm(df1(x) - lam*dh1, inf);
fprintf('   Sistema: 4x1 - lambda = 0, 2x2 - lambda = 0, x1 + x2 = 1\n');
fprintf('   x* = (%.4f, %.4f), lambda* = %.4f, f* = %.4f  (slide: (1/3, 2/3), 4/3, 2/3)\n', ...
        x, lam, f1(x));
fprintf('   grad f(x*) = (%.4f, %.4f) = lambda*(1,1);  |grad_x L| = %.1e, |h| = %.1e\n', ...
        df1(x), res, abs(h1(x)));
% O que significa lambda?  x1 + x2 = 1 + delta  ->  df*/d(delta) = lambda*
sel = @(v, i) v(i);                       % (MATLAB não indexa o resultado de uma chamada)
fd = @(delta) f1(sel(lagrange_lin(diag([4 2]), [0; 0], [1 1], 1 + delta), 1:2));
dl = 1e-4;
decl = (fd(dl) - fd(-dl))/(2*dl);
fprintf('   Sensibilidade: f*(delta) = (2/3)(1+delta)^2, declive em 0 (dif. centrais) = %.6f = lambda*\n', decl);
c = [c, perto(x, [1/3; 2/3], tolA), perto(lam, 4/3, tolA), perto(f1(x), 2/3, tolA), ...
     res <= tolA, abs(h1(x)) <= tolA, perto(decl, 4/3, 1e-8)];
fprintf('   Geometria (1): %s %s | x*: %s | lambda*: %s | f*: %s | grad_x L = 0: %s | h = 0: %s | declive = lambda*: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 2. Candidatos
f2  = @(x) x(1) + x(2);
df2 = @(x) [1; 1];
h2  = @(x) x(1)^2 + x(2)^2 - 2;
dh2 = @(x) [2*x(1); 2*x(2)];
fprintf('\n2. min x1 + x2  s.a.  h = x1^2 + x2^2 - 2 = 0\n');
fprintf('   (1,1) = lambda (2x1, 2x2)  =>  x1 = x2 = 1/(2 lambda);  na restrição, x1 = x2 = +-1\n');
dT = [1; -1]/sqrt(2);                       % direção tangente (a mesma nos dois)
c = [];
cand = [-1 -1; 1 1];  slam = [-1/2, 1/2];  sf = [-2, 2];  nome = {'mínimo', 'máximo'};
for i = 1:2
  x = cand(i, :)';
  lam = 1/(2*x(1));
  res = norm(df2(x) - lam*dh2(x), inf);
  H = -2*lam*eye(2);                        % Hessiana de L em x
  q = dT'*H*dT;                             % 2.ª ordem em T
  fprintf('   x = (%2d, %2d): lambda = %5.2f, f = %2d, |grad_x L| = %.1e, |h| = %.1e, d''H_L d = %5.2f -> %s\n', ...
          x, lam, f2(x), res, abs(h2(x)), q, nome{i});
  c = [c, perto(lam, slam(i), tolA), perto(f2(x), sf(i), tolA), res <= tolA, ...
       abs(h2(x)) <= tolA, sign(q) == -sign(lam)];
end
fprintf('   lambda, f, grad_x L = 0, h = 0 e classificação nos dois candidatos: %s\n', simnao{all(c) + 1});
ok = ok && all(c);

% ------------------------------------------------ 3. A lata
V = 330;
A3  = @(z) 2*pi*z(1)^2 + 2*pi*z(1)*z(2);            % z = (r, h), h = altura
dA3 = @(z) [4*pi*z(1) + 2*pi*z(2); 2*pi*z(1)];
h3  = @(z) pi*z(1)^2*z(2) - V;
dh3 = @(z) [2*pi*z(1)*z(2); pi*z(1)^2];
fprintf('\n3. Lata: min A = 2 pi r^2 + 2 pi r h  s.a.  pi r^2 h = V = %g cm^3\n', V);
fprintf('   dL/dh = 2 pi r - lambda pi r^2 = 0 => lambda = 2/r;  dL/dr = 0 => h = 2r;  pi r^2 (2r) = V\n');
r = (V/(2*pi))^(1/3);  h = 2*r;  lam = 2/r;  z = [r; h];
res = norm(dA3(z) - lam*dh3(z), inf);
HL = pi*[4 - 2*lam*h, 2 - 2*lam*r; 2 - 2*lam*r, 0];  % Hessiana de L = pi [-4 -2; -2 0]
dT = [1; -4];                                        % dh3(z)'*dT = 0
fprintf('   r* = %.4f cm, h* = %.4f cm, A* = %.4f cm^2, lambda* = 2/r* = %.4f cm^2/cm^3\n', ...
        r, h, A3(z), lam);
fprintf('   |grad_x L| = %.1e, |h| = %.1e;  2.ª ordem: d = (1,-4) em T, d''H_L d = %.4f = 12 pi > 0 (mínimo)\n', ...
        res, abs(h3(z)), dT'*HL*dT);
c = [confere(r, 3.745, 3), confere(r, 3.74, 2), confere(h, 7.49, 2), confere(A3(z), 264.4, 1), ...
     confere(lam, 0.534, 3), res <= 1e-9*norm(dA3(z)), abs(h3(z)) <= 1e-9*V, ...
     perto(dT'*HL*dT, 12*pi, 1e-9)];
fprintf('   r* = 3.745 (3.74): %s %s | h* = 7.49: %s | A* = 264.4: %s | lambda* = 0.534: %s | grad_x L = 0: %s | h = 0: %s | 12 pi: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 4. Várias restrições
f4  = @(x) x(1)^2 + x(2)^2 + x(3)^2;
df4 = @(x) 2*x;
A4  = [1 1 1; 1 -1 0];  b4 = [3; 1];            % h = A4 x - b4
fprintf('\n4. min x1^2 + x2^2 + x3^2  s.a.  x1 + x2 + x3 = 3,  x1 - x2 = 1\n');
z = lagrange_lin(2*eye(3), zeros(3, 1), A4, b4);
x = z(1:3);  lam = z(4:5);
res = norm(df4(x) - A4'*lam, inf);
d = [1; 1; -2];
fprintf('   x* = (%.4f, %.4f, %.4f), lambda* = (%.4f, %.4f), f* = %.4f  (slide: (3/2, 1/2, 1), (2, 1), 7/2)\n', ...
        x, lam, f4(x));
fprintf('   grad f(x*) = (%g, %g, %g) = 2 (1,1,1) + 1 (1,-1,0);  grad f''d = %g (d = (1,1,-2));  |grad_x L| = %.1e, |h| = %.1e\n', ...
        df4(x), df4(x)'*d, res, norm(A4*x - b4, inf));
c = [perto(x, [3/2; 1/2; 1], tolA), perto(lam, [2; 1], tolA), perto(f4(x), 7/2, tolA), ...
     perto(df4(x)'*d, 0, tolA), res <= tolA, norm(A4*x - b4, inf) <= tolA];
fprintf('   x*: %s | lambda*: %s | f*: %s | grad f''d = 0: %s | grad_x L = 0: %s | h = 0: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 5. Quando Lagrange falha
fprintf('\n5. min x1 + x2 + x3^2  s.a.  x1 - 1 = 0,  x1^2 + x2^2 - 1 = 0;  x* = (1, 0, 0)\n');
x = [1; 0; 0];
gf = [1; 1; 2*x(3)];
J = [1 0 0; 2*x(1) 2*x(2) 0];                   % linhas: grad h1', grad h2'
lam = J' \ gf;                                  % melhor lambda no sentido dos mínimos quadrados
resmin = norm(gf - J'*lam);
fprintf('   grad h1 = (%g,%g,%g), grad h2 = (%g,%g,%g): característica %d (dependentes)\n', ...
        J(1, :), J(2, :), rank(J));
fprintf('   min_lambda |grad f - J''lambda| = %.4f: dL/dx2 = 1 para todos os lambda -> sistema sem solução\n', resmin);
c = [rank(J) == 1, perto(resmin, 1, tolA)];
fprintf('   gradientes dependentes: %s | sem multiplicadores (resíduo 1): %s\n', simnao{c + 1});
ok = ok && all(c);

% =================================================== (b) com um solver
fprintf('\n=============== (b) com um solver ===============\n');
temFmincon = exist('fmincon', 'file') == 2;
temSqp = exist('sqp', 'file') == 2;
if temFmincon
  fprintf('Solver: fmincon (Optimization Toolbox). lambda_UC = -lambda.eqnonlin (L_fmincon = f + lambda''ceq).\n');
  opts = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12, 'TolCon', 1e-12);
elseif temSqp
  fprintf('Solver: sqp (GNU Octave; fmincon não existe). O sqp usa h(x) = 0 e L = f - lambda''h,\n');
  fprintf('a convenção da UC: lambda_UC = lambda (sem troca de sinal).\n');
  fprintf('Com o fmincon (MATLAB) seria lambda_UC = -lambda.eqnonlin (deck 4.2, «No computador»).\n');
else
  fprintf('Nem fmincon nem sqp disponíveis: parte (b) não corrida.\n');
end
if temFmincon || temSqp
  probs = {f1, @(x) h1(x),     [0; 0],    [1/3; 2/3],    4/3,    '1. exemplo resolvido';
           A3, @(z) h3(z),     [3; 6],    [],            [],     '3. lata (x0 = (3, 6))';
           f4, @(x) A4*x - b4, [0; 0; 0], [3/2; 1/2; 1], [2; 1], '4. várias restrições'};
  for i = 1:size(probs, 1)
    [fo, hc, x0, xs, ls, nome] = probs{i, :};
    if temFmincon
      [x, fx, flag, out, lambda] = fmincon(fo, x0, [], [], [], [], [], [], ...
                                           @(x) deal([], hc(x)), opts);
      lam = -lambda.eqnonlin;  nf = out.funcCount;  nit = out.iterations;
    else
      [x, fx, flag, nit, nf, lambda] = sqp(x0, fo, hc, []);
      lam = lambda;
    end
    fprintf('   %s: x = (%s), f = %.4f, lambda_UC = (%s)  [flag %d, %d it., %d aval.]\n', ...
            nome, sprintf(' %.4f', x), fx, sprintf(' %.4f', lam), flag, nit, nf);
    if isempty(xs)          % a lata: comparar com os valores do slide
      c = [confere(x(1), 3.745, 3), confere(x(2), 7.49, 2), confere(fx, 264.4, 1), ...
           confere(lam, 0.534, 3)];
      fprintf('      r*: %s | h*: %s | A*: %s | lambda*: %s\n', simnao{c + 1});
    else
      c = [perto(x, xs, tolB), perto(lam, ls, tolB)];
      fprintf('      x*: %s | lambda*: %s\n', simnao{c + 1});
    end
    ok = ok && all(c);
  end
end

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
