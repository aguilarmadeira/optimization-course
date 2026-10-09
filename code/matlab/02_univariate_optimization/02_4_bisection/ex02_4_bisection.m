% EX02_4_BISECTION  Reproduz os exemplos do deck 2.4 (bisseção, paragem |f'(z)| <= eps).
%
%   Exemplo 2 (exemplo-guia): f(x) = x^2 + 54/x em [1,4], eps = 0.05:
%             7 iterações, z = 2.9922, f'(z) = -0.0470, n_g = 7 + 2 = 9, n_f = 0.
%             Variante com 15 reduções: x2 - x1 = 3/2^15 = 9.2e-5 (secção áurea: 0.0022).
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
f  = @(x) x.^2 + 54./x;
df = @(x) 2*x - 54./x.^2;
a = 1;  b = 4;  eps_ = 0.05;
fprintf('\nExemplo 2: f(x) = x^2 + 54/x em [%g,%g], f''(%g) = %g < 0, f''(%g) = %g > 0, eps = %g\n', ...
        a, b, a, df(a), b, df(b), eps_);
[x, fx, info] = Bisection(f, df, a, b, eps_, [], [], true);
fprintf('|f''(z)| = %.4f <= eps na iteração %d: x* ~ z = %.4f, f(z) = %.4f (exato: x* = 3)\n', ...
        abs(info.history(end, 6)), info.nit, x, fx);
fprintf('n_g = %d + 2 = %d, n_f = %d\n', info.nit, info.ngev, info.nfev);
[xv, ~, iv] = Bisection(f, df, a, b, 0, [], 3/2^15);
tau = (sqrt(5) - 1)/2;  Lgs = 3*tau^15;
fprintf('Variante, 15 reduções (tolx = 3/2^15): %d iterações, x2 - x1 = %.1e, x = %.4f, n_g = %d\n', ...
        iv.nit, iv.L, xv, iv.ngev);
fprintf('Secção áurea, mesmo [1,4], 15 reduções: L = %.4f, %.0f vezes maior\n', Lgs, Lgs/iv.L);
T2 = [1 1.0000 4.0000 2.5000 27.8500 -3.6400
      2 2.5000 4.0000 3.2500 27.1779  1.3876
      3 2.5000 3.2500 2.8750 27.0482 -0.7831
      4 2.8750 3.2500 3.0625 27.0116  0.3674
      5 2.8750 3.0625 2.9688 27.0030 -0.1895
      6 2.9688 3.0625 3.0156 27.0007  0.0933
      7 2.9688 3.0156 2.9922 27.0002 -0.0470];
c = [isequal(size(info.history), [7 6]) && confere(info.history, T2, 4), info.flag == 0, ...
     confere(x, 2.9922, 4), info.ngev == 9, info.nfev == 0, ...
     iv.nit == 15 && iv.flag == 2, confere(iv.L, 9.2e-5, 6), iv.ngev == 17, ...
     confere(Lgs, 0.0022, 4), round(Lgs/iv.L) == 24];
fprintf(['  tabela k = 1..7: %s | parou por |f''| <= eps: %s | z = 2.9922: %s | n_g = 9: %s | n_f = 0: %s' ...
         ' | variante 15 it.: %s | 9.2e-5: %s | n_g = 17: %s | L áurea: %s | 24 vezes: %s\n'], simnao{c + 1});
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
