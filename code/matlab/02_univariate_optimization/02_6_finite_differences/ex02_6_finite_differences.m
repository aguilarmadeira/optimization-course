% EX02_6_FINITE_DIFFERENCES  Reproduz os exemplos do deck 2.6.
%
%   Conta à mão: f(x) = x^2 + 54/x em x = 2.5, h = 0.025 (= 0.01|x|):
%       progressiva -3.5295, centrada -3.64086, f'' ~ 8.9127; passo de
%       Newton x1 = 2.9085 (exato 2.9084); com 6 casas f'' daria 8.9136.
%   «Então basta fazer h cada vez menor?»: tabela para h = 0.1 ... 1e-8.
%   Newton numérico de x0 = 1, h = 0.01|x|, tolg = 1e-4: estagna em 3.00010,
%       para em k = 6 e devolve x_7; n_f = 7*3 + 1 = 22.
%
%   Imprime as tabelas como nos slides, os valores finais e as contagens, e
%   no fim compara com os slides (vírgula decimal nos slides, ponto aqui).
%
%   Otimização — deck 2.6.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.6: derivadas numéricas; Newton numérico\n');

f = @(x) x.^2 + 54./x;
dfex = @(x) 2*x - 54./x.^2;          % só para comparar
ddfex = @(x) 2 + 108./x.^3;

% ---------------------------------------------------------- conta à mão
x = 2.5;
[df, ddf, info] = NumDerivs(f, x);
fprintf('\nConta à mão: x = %g, h = %g (= 0.01|x|); exatos f'' = %.2f, f'''' = %.3f\n', ...
        x, info.h, dfex(x), ddfex(x));
fprintf('  f(x-h) = %.8f, f(x) = %.8f, f(x+h) = %.8f   (n_f = %d)\n', info.fm, info.f0, info.fp, info.nfev);
fprintf('  progressiva %.4f (erro %.2f); regressiva %.4f; centrada %.5f (erro %.0e)\n', ...
        info.df_fwd, abs(info.df_fwd - dfex(x)), info.df_bwd, df, abs(df - dfex(x)));
fprintf('  f''''_num = %.4f (erro %.0e)\n', ddf, abs(ddf - ddfex(x)));
x1 = x - df/ddf;
fprintf('  passo de Newton com estes valores: x1 = %.4f (exato: %.4f)\n', x1, x - dfex(x)/ddfex(x));
r6 = @(v) round(v*1e6)/1e6;          % valores com 6 casas
ddf6 = (r6(info.fp) - 2*info.f0 + r6(info.fm)) / info.h^2;
fprintf('  com 6 casas (%.6f; %.6f): f''''_num = %.4f (cancelamento)\n', r6(info.fm), r6(info.fp), ddf6);
c = [confere([info.fm info.f0 info.fp], [27.94380682 27.85 27.76176361], 8), ...
     confere(info.df_fwd, -3.5295, 4), confere(abs(info.df_fwd - dfex(x)), 0.11, 2), ...
     confere(df, -3.64086, 5), confere(abs(df - dfex(x)), 9e-4, 4), ...
     confere(ddf, 8.9127, 4), confere(x1, 2.9085, 4), confere(x - dfex(x)/ddfex(x), 2.9084, 4), ...
     confere(ddf6, 8.9136, 4), info.nfev == 3];
fprintf('  valores de f: %s | progressiva: %s | erro 0.11: %s | centrada: %s | erro 9e-4: %s | f'''': %s | x1: %s | exato: %s | 6 casas: %s | n_f = 3: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------- h cada vez menor?
fprintf('\nh cada vez menor? (dupla precisão, x = 2.5)\n');
fprintf('%8s %12s %9s %12s %9s\n', 'h', 'df central', 'erro', 'ddf central', 'erro');
hh = [0.1 0.025 1e-3 1e-6 1e-8];
R = zeros(numel(hh), 4);
for i = 1:numel(hh)
  [d1, d2] = NumDerivs(f, x, hh(i));
  R(i, :) = [d1, abs(d1 - dfex(x)), d2, abs(d2 - ddfex(x))];
  fprintf('%8.3g %12.5f %9.1e %12.4f %9.1e\n', hh(i), R(i, :));
end
Sd  = [-3.65385 -3.64086 -3.64000 -3.64000 -3.64000];
Se1 = [1.4e-2 9e-4 1.4e-6 1e-9 5e-8];   de1 = [3 4 7 9 8];
Sdd = [8.9231 8.9127 8.9120 8.9067 0.0000];
Se2 = [1e-2 7e-4 1e-6 5e-3 9];         de2 = [2 4 6 3 0];
c = [confere(R(:, 1), Sd, 5), confere(R(:, 2), Se1, de1), confere(R(:, 3), Sdd, 4), ...
     confere(R(:, 4), Se2, de2), R(5, 3) == 0];
fprintf('  df: %s | erros df: %s | ddf: %s | erros ddf: %s | ddf = 0 com h = 1e-8: %s\n', simnao{c + 1});
ok = ok && all(c);
fprintf('  ordens de grandeza do h ótimo: eps^(1/3) = %.0e (df), eps^(1/4) = %.1e (ddf)\n', ...
        eps^(1/3), eps^(1/4));
c = [confere(eps^(1/3), 6e-6, 6), confere(eps^(1/4), 1.2e-4, 5)];
fprintf('  eps^(1/3) ~ 6e-6: %s | eps^(1/4) ~ 1.2e-4: %s\n', simnao{c + 1});
ok = ok && all(c);

% ---------------------------------------------- Newton numérico, x0 = 1
x0 = 1;  tolg = 1e-4;  kmax = 50;
fprintf('\nNewton numérico: f(x) = x^2 + 54/x, x0 = %g, h = max(0.01|x|, 1e-4), tolg = %g\n', x0, tolg);
[x, fx, info] = NewtonNumeric(f, x0, tolg, kmax);
H = info.history;
fprintf('%3s %9s %8s %11s %10s %9s\n', 'k', 'x_k', 'h', 'df_num', 'ddf_num', 'x_k+1');
for r = 1:size(H, 1)
  fprintf('%3d %9.5f %8.4f %+11.4e %10.4f %9.5f\n', H(r, :));
end
fprintf('Para em k = %d e devolve x_%d = %.5f (estagna: x* = 3), f = %.4f\n', ...
        H(end, 1), info.nit, x, fx);
fprintf('n_f = %d linhas x 3 = %d, + 1 (f(x_%d)) = %d\n', info.nit, 3*info.nit, info.nit, info.nfev);
[d3, ~] = NumDerivs(f, 3, 0.03);
fprintf('Com h = 0.03, f''_num(3) = %.1e (não é 0)\n', d3);
[xa, ~, infoa] = NewtonNumeric(f, x, 1e-3, 1, 1e-4);
fprintf('Reduzindo h para 1e-4|x| no fim, uma iteração: x = %.7f (n_f = %d)\n', xa, infoa.nfev);
Tn = [1.00000 0.0100 -52.0054 110.0108 1.47273
      1.47273 0.0147 -21.9541  35.8141 2.08573
      2.08573 0.0209  -8.2428  13.9040 2.67857
      2.67857 0.0268  -2.1700   7.6203 2.96334
      2.96334 0.0296  -0.2233   6.1507 2.99965
      2.99965 0.0300  -0.0027   6.0018 3.00010
      3.00010 0.0300  -4e-7     6.0000 3.00010];
c = [size(H, 1) == 7 && confere(H(:, [2 6]), Tn(:, [1 5]), 5) && confere(H(:, [3 4 5]), Tn(:, [2 3 4]), 4), ...
     confere(H(7, 4), -4e-7, 7), H(end, 1) == 6, info.nit == 7, info.nfev == 22, ...
     confere(x, 3.00010, 5), confere(d3, -6e-4, 4), confere(xa, 3, 7)];
fprintf('  tabela k = 0..6: %s | df_6 = -4e-7: %s | para em k = 6: %s | devolve x_7: %s | n_f = 22: %s | x = 3.00010: %s | f''_num(3) = -6e-4: %s | 3.0000000: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
