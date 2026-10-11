% EX02_4_BISECTION  Reproduz os exemplos do deck 2.4 (bisseção, paragem |f'(z)| <= eps).
%
%   Exemplo 2 (exemplo-guia, 1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x em [1,2],
%             eps = 0.01: 6 iterações, z = 1.6094 (1.609375), f'(z) = -0.00001,
%             n_g = 6 + 2 = 8, n_f = 0; com eps = 0.03 pararia na it. 3 (z = 1.625).
%             Variante com 15 reduções: x2 - x1 = 1/2^15 = 3.1e-5 (secção áurea: 0.0007).
%   Exemplo 1: f(x) = -2 sin(x) + x^2/16 em [-1,3], eps = 0.03 (e eps = 0.05).
%   «Para resolver», eps = 0.01: x + 31/(x+3) em [0,5] e -3 sin(x) + x^2/8 em [0,4].
%   «Quantas iterações?» (variante): [0,2], tolx = 1e-3 -> k = 11; tolx = 1e-4 -> k = 15.
%
%   Imprime as tabelas das iterações (k = 1, 2, ...), os valores finais e as
%   contagens, e no fim compara com os slides (vírgula decimal nos slides,
%   ponto aqui). A coluna f(z) é só para leitura: o método não a usa.
%
%   Otimização — deck 2.4.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.4: bisseção (paragem |f''(z)| <= eps)\n');

% ------------------------------------------------------------ Exemplo 2
f  = @(x) 8*exp(1 - x) + 7*log(x);
df = @(x) -8*exp(1 - x) + 7./x;
a = 1;  b = 2;  eps_ = 0.01;
fprintf('\nExemplo 2 (exemplo-guia, 1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x em [%g,%g], f''(%g) = %g < 0, f''(%g) = %.4f > 0, eps = %g\n', ...
        a, b, a, df(a), b, df(b), eps_);
[x, fx, info] = Bisection(f, df, a, b, eps_, [], [], true);
fprintf('|f''(z)| = %.5f <= eps na iteração %d: x* ~ z = %.4f, f(z) = %.5f (exato: x* = 1.6093811)\n', ...
        abs(info.history(end, 6)), info.nit, x, fx);
fprintf('n_g = %d + 2 = %d, n_f = %d\n', info.nit, info.ngev, info.nfev);
[x3, ~, i3] = Bisection(f, df, a, b, 0.03);
fprintf('Com eps = 0.03: termina na iteração %d, z = %.4f\n', i3.nit, x3);
[xv, ~, iv] = Bisection(f, df, a, b, 0, [], 1/2^15);
tau = (sqrt(5) - 1)/2;  Lgs = tau^15;
fprintf('Variante, 15 reduções (tolx = 1/2^15): %d iterações, x2 - x1 = %.1e, x = %.4f, n_g = %d\n', ...
        iv.nit, iv.L, xv, iv.ngev);
fprintf('Secção áurea, mesmo [1,2], 15 reduções: L = %.4f, %.0f vezes maior\n', Lgs, Lgs/iv.L);
% tabela do slide: k  x1  x2  z (4 casas)  f(z)  f'(z) (5 casas)
T2 = [1 1.0000 2.0000 1.5000 7.69050 -0.18558
      2 1.5000 2.0000 1.7500 7.69624  0.22107
      3 1.5000 1.7500 1.6250 7.68065  0.02560
      4 1.5000 1.6250 1.5625 7.68227 -0.07826
      5 1.5625 1.6250 1.5938 7.68065 -0.02586
      6 1.5938 1.6250 1.6094 7.68045 -0.00001];
D2 = repmat([0 4 4 4 5 5], 6, 1);
c = [isequal(size(info.history), [6 6]) && confere(info.history, T2, D2), info.flag == 0, ...
     confere(x, 1.6094, 4), info.ngev == 8, info.nfev == 0, i3.nit == 3 && x3 == 1.625, ...
     iv.nit == 15 && iv.flag == 2, confere(iv.L, 3.1e-5, 6), iv.ngev == 17, ...
     confere(Lgs, 0.0007, 4), round(Lgs/iv.L) == 24];
fprintf(['  tabela k = 1..6: %s | parou por |f''| <= eps: %s | z = 1.6094: %s | n_g = 8: %s | n_f = 0: %s | eps = 0.03: 3 it.: %s' ...
         ' | variante 15 it.: %s | 3.1e-5: %s | n_g = 17: %s | L áurea: %s | 24 vezes: %s\n'], simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ Exemplo 1
f  = @(x) -2*sin(x) + x.^2/16;
df = @(x) -2*cos(x) + x/8;
a = -1;  b = 3;  eps_ = 0.03;
fprintf('\nExemplo 1: f(x) = -2 sin(x) + x^2/16 em [%g,%g], f''(%g) = %.2f < 0, f''(%g) = %.2f > 0, eps = %g\n', ...
        a, b, a, df(a), b, df(b), eps_);
[x, fx, info] = Bisection(f, df, a, b, eps_, [], [], true);
xs = Bisection([], df, a, b, 1e-12);                 % referência («exato»)
fprintf('|f''(z)| = %.4f <= eps na iteração %d: x* ~ z = %.4f (exato: %.4f, erro %.4f)\n', ...
        abs(info.history(end, 6)), info.nit, x, xs, abs(x - xs));
[x5, ~, i5] = Bisection(f, df, a, b, 0.05);
fprintf('Com eps = 0.05: termina na iteração %d, z = %.4f, |f''(z)| = %.4f, erro %.3f\n', ...
        i5.nit, x5, abs(i5.history(end, 6)), abs(x5 - xs));
T1 = [1 -1.0000 3.0000 1.0000 -1.6204 -0.9556
      2  1.0000 3.0000 2.0000 -1.5686  1.0823
      3  1.0000 2.0000 1.5000 -1.8544  0.0460
      4  1.0000 1.5000 1.2500 -1.8003 -0.4744
      5  1.2500 1.5000 1.3750 -1.8436 -0.2172
      6  1.3750 1.5000 1.4375 -1.8531 -0.0861
      7  1.4375 1.5000 1.4688 -1.8548 -0.0201];
c = [isequal(size(info.history), [7 6]) && confere(info.history, T1, 4), confere([df(-1) df(3)], [-1.21 2.35], 2), ...
     confere(x, 1.4688, 4), confere(xs, 1.4783, 4), ...
     i5.nit == 3 && confere([x5 abs(i5.history(end, 6)) abs(x5 - xs)], [1.5 0.0460 0.022], [1 4 3])];
fprintf('  tabela k = 1..7: %s | f''(-1), f''(3): %s | z = 1.4688: %s | exato 1.4783: %s | eps = 0.05: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------- Para resolver
fprintf('\nPara resolver (eps = 0.01):\n');
f  = @(x) x + 31./(x + 3);         df = @(x) 1 - 31./(x + 3).^2;
[x3, f3, i3] = Bisection(f, df, 0, 5, 0.01);
fprintf('  x + 31/(x+3) em [0,5]:       %2d iterações, z = %.4f, f(z) = %.4f  (exato: sqrt(31) - 3 = %.4f; erro %.4f)\n', ...
        i3.nit, x3, f3, sqrt(31) - 3, abs(x3 - (sqrt(31) - 3)));
f  = @(x) -3*sin(x) + x.^2/8;      df = @(x) -3*cos(x) + x/4;
[x4, f4, i4] = Bisection(f, df, 0, 4, 0.01);
x4s = Bisection([], df, 0, 4, 1e-12);
[~, ~, i4v] = Bisection([], df, 0, 4, 0, [], 4/2^15);
fprintf('  -3 sin(x) + x^2/8 em [0,4]:  %2d iterações, z = %.4f, f(z) = %.4f  (x* = %.4f; variante, 15 reduções: x2 - x1 = %.1e)\n', ...
        i4.nit, x4, f4, x4s, i4v.L);
c = [i3.nit == 6 && confere([x3 f3], [2.5781 8.1355], 4), confere(sqrt(31) - 3, 2.5678, 4), ...
     confere(abs(x3 - (sqrt(31) - 3)), 0.0104, 4), ...
     i4.nit == 10 && confere([x4 f4], [1.4492 -2.7153], 4), confere(x4s, 1.4497, 4), ...
     i4v.nit == 15 && confere(i4v.L, 1.2e-4, 5)];
fprintf('  1.º: %s | exato: %s | erro 0.0104: %s | 2.º: %s | x* = 1.4497: %s | 1.2e-4: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ Quantas iterações?
kk = ceil(log2(2 ./ [1e-3 1e-4]));
dq = @(x) x - 0.3;                                  % um f' qualquer com zero em [0,2]
[~, ~, q1] = Bisection([], dq, 0, 2, 0, [], 1e-3);
[~, ~, q2] = Bisection([], dq, 0, 2, 0, [], 1e-4);
fprintf(['\nQuantas iterações? (variante) [0,2]: log2(2000) = %.2f -> k = %d (tolx = 1e-3); k = %d (tolx = 1e-4);' ...
         ' a função faz %d e %d\n'], log2(2000), kk, q1.nit, q2.nit);
c = [confere(log2(2000), 10.97, 2), kk(1) == 11, kk(2) == 15, isequal([q1.nit q2.nit], kk)];
fprintf('  log2 2000 = 10.97: %s | k = 11: %s | k = 15: %s | função: %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
