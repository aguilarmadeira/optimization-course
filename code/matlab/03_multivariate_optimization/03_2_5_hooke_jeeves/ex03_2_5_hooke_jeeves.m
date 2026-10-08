% EX03_2_5_HOOKE_JEEVES  Reproduz os exemplos do deck 3.2.5 (Hooke–Jeeves,
%   algoritmo das aulas: Deb, 2012).
%
%   f(x) = 3 x1^2 + x2^2 - 12 x1 - 8 x2  (mínimo em (2,4), f = -28):
%     «Uma exploração à mão»: de x0 = (1,1), f = -16, com Delta = (0.5, 0.5):
%       x1: f+ = f(1.5,1) = -18.25, f- = f(0.5,1) = -12.25: fica (1.5,1);
%       x2: f+ = f(1.5,1.5) = -21, f- = f(1.5,0.5) = -15: fica (1.5,1.5); 4 avaliações;
%     «Exemplo completo»: x0 = (1,1), Delta = (0.5,0.5), R = 2, eps = 0.2:
%       padrão (2,2) -> (2,2.5), sucesso; (2.5,3.5) -> (2,4), sucesso;
%       (2,5.5) -> (2,5), insucesso: passo 3, Delta = (0.25,0.25); exploração em
%       (2,4) falha; Delta = (0.125,0.125); falha; ||Delta|| = 0.177 < eps: pára em
%       (2,4), f = -28, n_f = 28.
%   Himmelblau a partir de (0,0), Delta = (0.5,0.5), R = 2, eps = 0.2 (Deb, 2012,
%     exercício 3.3.3): X1 = (0.5,0.5), X2 = (1.5,1.5), X3 = (3,2); o padrão
%     (4.5,2.5) -> (4,2), f = 50, falha; pára em (3,2), f = 0.
%   Rosenbrock a partir de (-1.5, 2), Delta = 0.5, R = 2, eps = 1e-6: n_f = 627 até
%     f < 1e-4 (com a versão de 1961 do código, variante '1961': 430).
%   Custo: 2n por exploração (avaliam-se +Delta_i e -Delta_i); 1 + 2n por
%     movimento em padrão.
%
%   Contagens: f(x0) conta; f do ponto base nunca é reavaliado (fica guardado).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 3.2.5.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

quad  = @(x) 3*x(1)^2 + x(2)^2 - 12*x(1) - 8*x(2);
rosen = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;

fprintf('Otimização — deck 3.2.5: método de Hooke–Jeeves\n');

% --------------------------------------------- «Uma exploração à mão»
xb = [1 1];  fb = quad(xb);  D = [0.5 0.5];
fprintf('\nExploração à mão: x0 = (%g,%g), f(x0) = %g, Delta = (%g,%g)\n', xb, fb, D);
[xe, fe, ie] = Exploratory(quad, xb, fb, D);
sg = '+-';  res = {'não é a melhor', 'a melhor: guarda-se'};
for r = 1:size(ie.history, 1)
  h = ie.history(r, :);  s = sg((3 - h(2))/2);
  fprintf('  x%d %s Delta%d: (%g, %g) -> f%s = %g, %s\n', h(1), s, h(1), h(3:4), s, h(5), ...
          res{h(6) + 1});
end
fprintf('x1 = (%g, %g), f = %g, com %d avaliações\n', xe, fe, ie.nfev);
c = [confere(fb, -16, 0), ...
     confere(ie.history(:, 3:5), [1.5 1 -18.25; 0.5 1 -12.25; 1.5 1.5 -21; 1.5 0.5 -15], 2), ...
     confere(xe, [1.5 1.5], 1), ie.nfev == 4];
fprintf('  f(1,1) = -16: %s | -18.25, -12.25, -21, -15: %s | x1 = (1.5,1.5): %s | 4 avaliações: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% --------------------------------------------- «Exemplo completo»
x0 = [1 1];  D0 = [0.5 0.5];  R = 2;  epsD = 0.2;
fprintf('\nExemplo completo: x0 = (%g,%g), Delta = (%g,%g), R = %g, eps = %g\n', x0, D0, R, epsD);
[x, fx, info] = HookeJeeves(quad, x0, 2, D0, epsD, [], true, [], 'deb', R);
fprintf('%s\n', info.message);
fprintf('solução x = (%g, %g), f = %g;  n_f = %d\n', x, fx, info.nfev);
H = info.history;
% movimentos 2 a 4 do slide (padrão): x^P(1) x^P(2) f(x^P)  x(1) x(2) f(x)  sucesso
Tpad = [2.0 2.0 -24.00  2.0 2.5 -25.75  1
        2.5 3.5 -27.00  2.0 4.0 -28.00  1
        2.0 5.5 -25.75  2.0 5.0 -27.00  0];
c = [info.nit == 6 && isequal(H(:, 2).', [0 1 1 1 0 0]), ...   % explor., 3 padrões, 2 explor.
     confere(H(1, 8:10), [1.5 1.5 -21], 2), ...
     confere(H(2:4, [5 6 7 8 9 10 14]), Tpad, 2), ...
     confere(H(5:6, 12), [0.25; 0.125], 4) && all(H(5:6, 14) == 0), ...
     confere(norm(info.P), 0.177, 3), ...
     confere(x, [2 4], 4) && confere(fx, -28, 4), info.nfev == 28];
fprintf(['  sequência de movimentos: %s | 1.ª exploração -> (1.5,1.5): %s | padrões (2,2), (2.5,3.5), (2,5.5): %s |\n' ...
         '  falha com Delta = 0.25, 0.125: %s | ||Delta|| = 0.177 < eps: %s | (2,4), f = -28: %s | n_f = 28: %s\n'], ...
        simnao{c + 1});
ok = ok && all(c);

% --------------------------------------------- Himmelblau (Deb, 2012, exercício 3.3.3)
himmel = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
fprintf('\nHimmelblau a partir de (0,0), Delta = (0.5,0.5), R = 2, eps = 0.2\n');
[x, fx, info] = HookeJeeves(himmel, [0 0], 2, [0.5 0.5], 0.2, [], true);
fprintf('%s\n', info.message);
H = info.history;
c = [confere(H(1, 8:10), [0.5 0.5 144.125], 3), confere(H(2, 8:10), [1.5 1.5 63.125], 3), ...
     confere(H(3, 8:10), [3 2 0], 3), ...
     confere(H(4, 5:10), [4.5 2.5 152.125 4 2 50], 3) && H(4, 14) == 0, ...
     confere(x, [3 2], 4) && fx == 0];
fprintf(['  X1 = (0.5,0.5): %s | X2 = (1.5,1.5): %s | X3 = (3,2): %s | ' ...
         'padrão (4.5,2.5) -> (4,2), f = 50, falha: %s | (3,2), f = 0: %s\n'], simnao{c + 1});
ok = ok && all(c);

% --------------------------------------------- Rosenbrock
x0 = [-1.5 2];  D0 = [0.5 0.5];  epsD = 1e-6;
fprintf('\nRosenbrock a partir de (%g,%g), Delta = (%g,%g), R = 2, eps = %g\n', x0, D0, epsD);
[x, fx, info] = HookeJeeves(rosen, x0, 2, D0, epsD);
fprintf('%s\n', info.message);
hit = find(info.ftrace < 1e-4, 1);
fprintf('primeira avaliação com f < 1e-4: n_f = %d\n', hit);
fprintf('no fim: x = (%.6f, %.6f), f = %.2e;  %d movimentos, n_f total = %d\n', ...
        x, fx, info.nit, info.nfev);
% custo de cada movimento: 2n (exploração), 1 + 2n (padrão), n = 2
H = info.history;  n = 2;
dn = diff([1; H(:, end)]);
ce = dn(H(:, 2) == 0);  cp = dn(H(:, 2) == 1);
fprintf('custo por exploração: %d a %d; por movimento em padrão: %d a %d\n', ...
        min(ce), max(ce), min(cp), max(cp));
[~, ~, info61] = HookeJeeves(rosen, x0, 2, D0, [1e-6 1e-6], [], false, [], '1961');
hit61 = find(info61.ftrace < 1e-4, 1);
fprintf('(versão de 1961, variante ''1961'': n_f = %d até f < 1e-4)\n', hit61);
c = [hit == 627, hit61 == 430, all(ce == 2*n), all(cp == 1 + 2*n)];
fprintf('  n_f = 627 até f < 1e-4: %s | 1961: 430: %s | n_f = 2n por exploração: %s | 1 + 2n por padrão: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
