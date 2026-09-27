% EX02_4_BISECTION  Reproduz os exemplos do deck 2.4 (bisseção).
%
%   Exemplo 2 (exemplo-guia): f(x) = x^2 + 54/x em [1,4], 15 reduções:
%             x ~ 3.0000, b - a = 3/2^15 = 9.2e-5, n_g = N + 2 = 17, n_f = 0.
%   Exemplo 1: f(x) = -2 sin(x) + x^2/16 em [-1,3], 15 reduções.
%   «Para resolver»: x + 31/(x+3) em [0,5] e -3 sin(x) + x^2/8 em [0,4].
%   «Quantas iterações?»: [0,2], tolx = 1e-3 -> k = 11; tolx = 1e-4 -> k = 15.
%
%   Imprime as tabelas das iterações (k = 0..14; os slides mostram k = 0..7),
%   os valores finais e as contagens, e no fim compara com os slides
%   (vírgula decimal nos slides, ponto aqui). A coluna f(x_k) é só para
%   leitura: o método não a usa.
%
%   Otimização — deck 2.4.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.4: bisseção\n');

% ------------------------------------------------------------ Exemplo 2
f  = @(x) x.^2 + 54./x;
df = @(x) 2*x - 54./x.^2;
a = 1;  b = 4;  N = 15;
fprintf('\nExemplo 2: f(x) = x^2 + 54/x em [%g,%g], f''(%g) = %g < 0, f''(%g) = %g > 0, N = %d\n', ...
        a, b, a, df(a), b, df(b), N);
[x, fx, info] = Bisection(f, df, a, b, N, true);
fprintf('Após %d reduções, ponto médio: x* ~ %.4f, f(x*) ~ %.4f, b - a = %.1e (exato: x* = 3)\n', ...
        N, x, fx, info.L);
fprintf('n_g = N + 2 = %d, n_f = %d\n', info.ngev, info.nfev);
tau = (sqrt(5) - 1)/2;  Lgs = 3*tau^15;
fprintf('Secção áurea, mesmo [1,4], 15 reduções: L = %.4f, %.0f vezes maior\n', Lgs, Lgs/info.L);
T2 = [1.0000 4.0000 2.5000 27.8500 -3.6400
      2.5000 4.0000 3.2500 27.1779  1.3876
      2.5000 3.2500 2.8750 27.0482 -0.7831
      2.8750 3.2500 3.0625 27.0116  0.3674
      2.8750 3.0625 2.9688 27.0030 -0.1895
      2.9688 3.0625 3.0156 27.0007  0.0933
      2.9688 3.0156 2.9922 27.0002 -0.0470
      2.9922 3.0156 3.0039 27.0000  0.0234];
c = [confere(info.history(1:8, 2:6), T2, 4), confere(x, 3.0000, 4), ...
     confere(info.L, 9.2e-5, 6), info.ngev == 17, info.nfev == 0, ...
     confere(Lgs, 0.0022, 4), round(Lgs/info.L) == 24];
fprintf('  tabela k = 0..7: %s | x* ~ 3.0000: %s | b - a = 9.2e-5: %s | n_g = 17: %s | n_f = 0: %s | L áurea: %s | 24 vezes: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ Exemplo 1
f  = @(x) -2*sin(x) + x.^2/16;
df = @(x) -2*cos(x) + x/8;
a = -1;  b = 3;  N = 15;
fprintf('\nExemplo 1: f(x) = -2 sin(x) + x^2/16 em [%g,%g], f''(%g) = %.2f < 0, f''(%g) = %.2f > 0, N = %d\n', ...
        a, b, a, df(a), b, df(b), N);
[x, fx, info] = Bisection(f, df, a, b, N, true);
fprintf('Após %d reduções, ponto médio: x* ~ %.4f, f(x*) ~ %.4f, b - a = %.1e\n', N, x, fx, info.L);
T1 = [-1.0000 3.0000 1.0000 -1.6204 -0.9556
       1.0000 3.0000 2.0000 -1.5686  1.0823
       1.0000 2.0000 1.5000 -1.8544  0.0460
       1.0000 1.5000 1.2500 -1.8003 -0.4744
       1.2500 1.5000 1.3750 -1.8436 -0.2172
       1.3750 1.5000 1.4375 -1.8531 -0.0861
       1.4375 1.5000 1.4688 -1.8548 -0.0201
       1.4688 1.5000 1.4844 -1.8548  0.0129];
c = [confere(info.history(1:8, 2:6), T1, 4), confere([df(-1) df(3)], [-1.21 2.35], 2), ...
     confere(x, 1.4783, 4), confere(fx, -1.8549, 4), confere(info.L, 1.2e-4, 5)];
fprintf('  tabela k = 0..7: %s | f''(-1), f''(3): %s | x*: %s | f(x*): %s | b - a = 1.2e-4: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------- Para resolver
fprintf('\nPara resolver (15 reduções, ponto médio final):\n');
f  = @(x) x + 31./(x + 3);         df = @(x) 1 - 31./(x + 3).^2;
[x3, f3] = Bisection(f, df, 0, 5, 15);
fprintf('  x + 31/(x+3) em [0,5]:       x* = %.4f, f(x*) = %.4f  (exato: sqrt(31) - 3 = %.4f)\n', ...
        x3, f3, sqrt(31) - 3);
f  = @(x) -3*sin(x) + x.^2/8;      df = @(x) -3*cos(x) + x/4;
[x4, f4, info4] = Bisection(f, df, 0, 4, 15);
fprintf('  -3 sin(x) + x^2/8 em [0,4]:  x* = %.4f, f(x*) = %.4f, b - a = %.1e\n', x4, f4, info4.L);
c = [confere([x3 f3], [2.5678 8.1355], 4), confere(sqrt(31) - 3, 2.5678, 4), ...
     confere([x4 f4], [1.4496 -2.7153], 4), confere(info4.L, 1.2e-4, 5)];
fprintf('  1.º: %s | exato: %s | 2.º: %s | b - a = 1.2e-4: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ Quantas iterações?
kk = ceil(log2(2 ./ [1e-3 1e-4]));
fprintf('\nQuantas iterações? [0,2]: log2(2000) = %.2f -> k = %d (tolx = 1e-3); k = %d (tolx = 1e-4)\n', ...
        log2(2000), kk);
c = [confere(log2(2000), 10.97, 2), kk(1) == 11, kk(2) == 15];
fprintf('  log2 2000 = 10.97: %s | k = 11: %s | k = 15: %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
