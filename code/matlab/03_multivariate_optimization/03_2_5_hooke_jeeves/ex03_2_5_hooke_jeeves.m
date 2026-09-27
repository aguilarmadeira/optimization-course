% EX03_2_5_HOOKE_JEEVES  Reproduz os exemplos do deck 3.2.5 (Hooke–Jeeves).
%
%   f(x) = 3 x1^2 + x2^2 - 12 x1 - 8 x2  (mínimo em (2,4), f = -28):
%     «Uma exploração à mão»: de xb = (1,1), f = -16, com P = (0.5, 0.5):
%       (1.5,1) -> -18.25, (1.5,1.5) -> -21, com 2 avaliações;
%     «Exemplo completo»: x0 = (1,1), a = 2, P0 = (0.5,0.5), T = (0.1,0.1):
%       padrão (2,2) -> (2,2.5), aceite; (2.5,3.5) -> (2,4), aceite;
%       (2,5.5) -> (2,5), rejeitado; exploração em (2,4) falha com P = 0.5,
%       0.25, 0.125; 0.0625 < T: para em (2,4), f = -28.
%   Rosenbrock a partir de (-1.5, 2), P0 = 0.5, a = 2, T = 1e-6 (como no
%     caderno cap3_comparacao.ipynb): n_f = 353 até f < 1e-4.
%   Custo: n <= n_f <= 2n por exploração; 1 + (n a 2n) por movimento de padrão.
%
%   Contagens: f(x0) conta; f(xb) nunca é reavaliado (fica guardado).
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
xb = [1 1];  fb = quad(xb);  P = [0.5 0.5];
fprintf('\nExploração à mão: xb = (%g,%g), f(xb) = %g, P = (%g,%g)\n', xb, fb, P);
[xe, fe, ie] = Exploratory(quad, xb, fb, P);
for r = 1:size(ie.history, 1)
  h = ie.history(r, :);
  sg = '+-';  res = {'não melhora', 'melhora: guarda-se'};
  fprintf('  x%d %s P%d: (%g, %g) -> f = %g, %s\n', h(1), sg((3 - h(2))/2), h(1), h(3:4), h(5), ...
          res{h(6) + 1});
end
fprintf('xe = (%g, %g), f = %g, com %d avaliações\n', xe, fe, ie.nfev);
c = [confere(fb, -16, 0), confere(ie.history(:, 3:5), [1.5 1 -18.25; 1.5 1.5 -21], 2), ...
     confere(xe, [1.5 1.5], 1), ie.nfev == 2];
fprintf('  f(1,1) = -16: %s | -18.25 e -21: %s | xe = (1.5,1.5): %s | 2 avaliações: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% --------------------------------------------- «Exemplo completo»
x0 = [1 1];  a = 2;  P0 = [0.5 0.5];  T = [0.1 0.1];
fprintf('\nExemplo completo: x0 = (%g,%g), a = %g, P0 = (%g,%g), T = (%g,%g)\n', x0, a, P0, T);
[x, fx, info] = HookeJeeves(quad, x0, a, P0, T, [], true);
fprintf('%s\n', info.message);
fprintf('solução x = (%g, %g), f = %g;  n_f = %d\n', x, fx, info.nfev);
H = info.history;
% linhas 2 a 4 do slide (padrão): xt(1) xt(2) f(xt)  xe'(1) xe'(2) f(xe')  aceite
Tpad = [2.0 2.0 -24.00  2.0 2.5 -25.75  1
        2.5 3.5 -27.00  2.0 4.0 -28.00  1
        2.0 5.5 -25.75  2.0 5.0 -27.00  0];
c = [info.nit == 7 && all(H(:, 2).' == [0 1 1 1 0 0 0]), ...   % explor., 3 padrões, 3 explor.
     confere(H(1, 8:10), [1.5 1.5 -21], 2), ...
     confere(H(2:4, [5 6 7 8 9 10 14]), Tpad, 2), ...
     confere(H(5:7, 12), [0.5; 0.25; 0.125], 4) && all(H(5:7, 14) == 0), ...
     confere(info.P, [0.0625 0.0625], 4), ...
     confere(x, [2 4], 4) && confere(fx, -28, 4)];
fprintf(['  sequência de movimentos: %s | 1.ª exploração -> (1.5,1.5): %s | padrões (2,2), (2.5,3.5), (2,5.5): %s |\n' ...
         '  falha com P = 0.5, 0.25, 0.125: %s | P = 0.0625 < T: %s | (2,4), f = -28: %s\n'], simnao{c + 1});
ok = ok && all(c);

% --------------------------------------------- Rosenbrock
x0 = [-1.5 2];  a = 2;  P0 = [0.5 0.5];  T = [1e-6 1e-6];
fprintf('\nRosenbrock a partir de (%g,%g), a = %g, P0 = (%g,%g), T = (%g,%g)\n', x0, a, P0, T);
[x, fx, info] = HookeJeeves(rosen, x0, a, P0, T);
fprintf('%s\n', info.message);
hit = find(info.ftrace < 1e-4, 1);
fprintf('primeira avaliação com f < 1e-4: n_f = %d\n', hit);
fprintf('no fim: x = (%.6f, %.6f), f = %.2e;  %d movimentos, n_f total = %d\n', ...
        x, fx, info.nit, info.nfev);
% custo de cada movimento: n a 2n (exploração), 1 + (n a 2n) (padrão), n = 2
H = info.history;  n = 2;
dn = diff([1; H(:, end)]);
ce = dn(H(:, 2) == 0);  cp = dn(H(:, 2) == 1);
fprintf('custo por exploração: %d a %d; por movimento de padrão: %d a %d\n', ...
        min(ce), max(ce), min(cp), max(cp));
[~, ~, info5] = HookeJeeves(rosen, x0, a, P0, [1e-5 1e-5]);
fprintf('(com T = 1e-5: n_f = %d até f < 1e-4 e n_f total = %d)\n', ...
        find(info5.ftrace < 1e-4, 1), info5.nfev);
c = [hit == 353, all(ce >= n & ce <= 2*n), all(cp >= 1 + n & cp <= 1 + 2*n)];
fprintf('  n_f = 353 até f < 1e-4: %s | n <= n_f <= 2n: %s | 1 + (n a 2n): %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
