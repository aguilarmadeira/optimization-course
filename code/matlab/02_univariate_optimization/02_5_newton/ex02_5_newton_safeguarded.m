% EX02_5_NEWTON_SAFEGUARDED  Reproduz o exemplo «Newton salvaguardado:
% bisseção + Newton» do deck 2.5.
%
%   f'(x) = arctan x (mínimo de f(x) = x arctan x - ln(1 + x^2)/2 em 0),
%   intervalo [-1, 2], x0 = 1.5 -- o caso em que Newton puro diverge
%   (página «Quando falha»). Um passo de bisseção, depois três de Newton:
%
%     k   x_k       passo de Newton               x_k+1       novo intervalo
%     0   1.5       -1.694: fora de (a,b) -> bis.  0.5         [-1; 0.5]
%     1   0.5       -0.0796: aceite                -0.0796     [-0.0796; 0.5]
%     2   -0.0796   3.4e-4: aceite                 3.4e-4      [-0.0796; 3.4e-4]
%     3   3.4e-4    -2.5e-11: aceite               -2.5e-11
%
%   Paragem |f'| < tolg = 1e-10 (a do script das figuras).
%   Imprime a tabela, as contagens e a comparação com Newton puro, e no fim
%   compara com os slides (vírgula decimal nos slides, ponto aqui).
%
%   Otimização — deck 2.5.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com s mostrado com m algarismos significativos?
confere_sig = @(v, s, m) all(abs(v(:) - s(:)) <= 0.5*10.^(floor(log10(abs(s(:)))) - m(:) + 1)*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.5: Newton salvaguardado (bisseção + Newton)\n');

f   = @(x) x.*atan(x) - 0.5*log(1 + x.^2);
df  = @(x) atan(x);
ddf = @(x) 1./(1 + x.^2);

% ---------------------------------------------- Newton puro diverge
x0 = 1.5;
fprintf('\nNewton puro para f''(x) = arctan x a partir de x0 = %g (3 iterações):\n', x0);
xs = x0;
for j = 1:3                                  % x_{k+1} = x_k - f'(x_k)/f''(x_k)
  xs(end + 1) = xs(end) - df(xs(end))/ddf(xs(end)); %#ok<SAGROW>
end
fprintf('  x_k: %.4f -> %.4f -> %.4f -> %.4f  (|x_k| cresce: diverge)\n', xs);
c = abs(xs(4)) > abs(xs(3)) && abs(xs(3)) > abs(xs(2)) && abs(xs(2)) > abs(xs(1));
fprintf('  diverge: %s\n', simnao{c + 1});
ok = ok && c;

% ---------------------------------------------- Newton salvaguardado
a = -1;  b = 2;  tolg = 1e-10;
fprintf('\nNewton salvaguardado em [%g, %g], x0 = %g, tolg = %g:\n', a, b, x0, tolg);
[x, fx, info] = NewtonSafeguarded(f, df, ddf, a, b, x0, tolg, 50, true);
H = info.history;
fprintf('%s; x* ~ %.1e, f(x*) = %.1e\n', info.message, x, fx);
fprintf('n_g = %d (f''(a), f''(b), f''(x0) e uma por iteração), n_H = %d, n_f = %d\n', ...
        info.ngev, info.nhev, info.nfev);
c = [size(H, 1) == 4, ...
     isequal(H(:, 4)', [1 0 0 0]), ...
     confere_sig(H(:, 2), [1.5; 0.5; -0.0796; 3.4e-4], [2; 1; 3; 2]), ...
     confere_sig(H(:, 3), [-1.694; -0.0796; 3.4e-4; -2.5e-11], [4; 3; 2; 2]), ...
     confere_sig(H(:, 5), [0.5; -0.0796; 3.4e-4; -2.5e-11], [1; 3; 2; 2]), ...
     confere_sig(H(1:3, 6:7), [-1 0.5; -0.0796 0.5; -0.0796 3.4e-4], [1 1; 3 1; 3 2]), ...
     info.nbis == 1, info.ngev == 7, info.nhev == 4];
fprintf(['  4 iterações: %s | 1 bisseção + 3 Newton: %s | x_k: %s | passos de Newton: %s | x_k+1: %s', ...
         ' | intervalos: %s | 1 bisseção: %s | n_g = 7: %s | n_H = 4: %s\n'], simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
