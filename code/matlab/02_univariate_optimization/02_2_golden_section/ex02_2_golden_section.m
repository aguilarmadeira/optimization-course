% EX02_2_GOLDEN_SECTION  Reproduz os exemplos do deck 2.2 (secção áurea).
%
%   Exemplo 2 (exemplo-guia, 1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x em [1,2],
%                             15 reduções (os slides mostram f com 5 casas).
%   Exemplo 1:                f(x) = -3 sin(x) + x^2/8 em [0,4], 15 reduções.
%   Slide «Quantas reduções?»: L0 = 1, tolx = 0.01  =>  N = 10, n_f = 12.
%
%   Imprime as tabelas das iterações (k = 0..15; os slides mostram k = 0..7),
%   os valores finais e as contagens, e no fim compara com os slides.
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 2.2.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.2: método da secção áurea\n');

% ------------------------------------------------------------ Exemplo 2
f = @(x) 8*exp(1 - x) + 7*log(x);
a = 1;  b = 2;  N = 15;
fprintf('\nExemplo 2 (exemplo-guia, 1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x em [%g,%g], N = %d reduções\n', a, b, N);
[x, fx, info] = GoldenSection(f, a, b, N, true);
L15 = info.history(end, 8);
fprintf('Após %d reduções: L_%d = %.4f, x* ~ %.4f, f(x*) ~ %.5f  (exato: x* = 1.6093811, f* = 7.6804459)\n', ...
        N, N, L15, x, fx);
fprintf('n_f = N + 2 = %d para localizar x*; com f no ponto médio final, n_f = N + 3 = %d\n', ...
        info.nfev_loc, info.nfev);

% tabela do slide (k = 0..7): a  b  x1  x2 (4 casas)  f(x1)  f(x2) (5 casas)  L_k (4 casas)
T2 = [1.0000 2.0000 1.3820 1.6180 7.72470 7.68051 1.0000
      1.3820 2.0000 1.6180 1.7639 7.68051 7.69947 0.6180
      1.3820 1.7639 1.5279 1.6180 7.68600 7.68051 0.3820
      1.5279 1.7639 1.6180 1.6738 7.68051 7.68381 0.2361
      1.5279 1.6738 1.5836 1.6180 7.68100 7.68051 0.1459
      1.5836 1.6738 1.6180 1.6393 7.68051 7.68118 0.0902
      1.5836 1.6393 1.6049 1.6180 7.68046 7.68051 0.0557
      1.5836 1.6180 1.5967 1.6049 7.68058 7.68046 0.0344];
D2 = repmat([4 4 4 4 5 5 4], 8, 1);
c = [confere(info.history(1:8, 2:8), T2, D2), confere(L15, 0.0007, 4), ...
     confere(x, 1.6095, 4), confere(fx, 7.68045, 5), ...
     info.nfev_loc == 17, info.nfev == 18];
fprintf('  tabela k = 0..7: %s | L_15: %s | x*: %s | f(x*): %s | n_f = 17: %s | n_f = 18: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ Exemplo 1
f = @(x) -3*sin(x) + x.^2/8;
a = 0;  b = 4;  N = 15;
fprintf('\nExemplo 1: f(x) = -3 sin(x) + x^2/8 em [%g,%g], N = %d reduções\n', a, b, N);
[x, fx, info] = GoldenSection(f, a, b, N, true);
L15 = info.history(end, 8);
fprintf('Após %d reduções: L_%d = %.4f, x* ~ %.4f, f(x*) ~ %.4f;  n_f = %d (+1 = %d)\n', ...
        N, N, L15, x, fx, info.nfev_loc, info.nfev);
T1 = [0.0000 4.0000 1.5279 2.4721 -2.7054 -1.0977 4.0000
      0.0000 2.4721 0.9443 1.5279 -2.3188 -2.7054 2.4721
      0.9443 2.4721 1.5279 1.8885 -2.7054 -2.4040 1.5279
      0.9443 1.8885 1.3050 1.5279 -2.6818 -2.7054 0.9443
      1.3050 1.8885 1.5279 1.6656 -2.7054 -2.6397 0.5836
      1.3050 1.6656 1.4427 1.5279 -2.7152 -2.7054 0.3607
      1.3050 1.5279 1.3901 1.4427 -2.7096 -2.7152 0.2229
      1.3901 1.5279 1.4427 1.4752 -2.7152 -2.7143 0.1378];
c = [confere(info.history(1:8, 2:8), T1, 4), confere(L15, 0.0029, 4), ...
     confere(x, 1.4489, 4), confere(fx, -2.7153, 4)];
fprintf('  tabela k = 0..7: %s | L_15: %s | x*: %s | f(x*): %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ «Quantas reduções?»
tau = (sqrt(5) - 1)/2;
L0 = 1;  tolx = 0.01;
Nmin = log(tolx/L0)/log(tau);
N = ceil(Nmin);
[~, ~, info] = GoldenSection(@(x) (x - 0.3).^2, 0, L0, N);
fprintf('\nQuantas reduções? L0 = %g, tolx = %g: N >= %.2f => N = %d, tau^N = %.4f, n_f = %d\n', ...
        L0, tolx, Nmin, N, tau^N, info.nfev_loc);
c = [confere(Nmin, 9.57, 2), N == 10, confere(tau^N, 0.0081, 4), info.nfev_loc == 12];
fprintf('  N >= 9.57: %s | N = 10: %s | 0.618^10: %s | n_f = 12: %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
