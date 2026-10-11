% EX02_5_NEWTON  Reproduz os exemplos do deck 2.5 (Newton-Raphson).
%
%   Exemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x0 = 1:
%       x1 = 2, x2 = 1.5331535, ...; 5 iterações até x* = 1.6093811 (erro 1.3e-13);
%       n_g = n_H = 1 por iteração (+ f'(x0)). Paragem |f'| < tolg com tolg = 1e-6
%       (o slide não fixa tolg; qualquer valor entre 2.2e-13 e 1.08e-6 dá as
%       mesmas 5 iterações). Na pergunta de abertura, x0 = 1.5: x1 = 1.607.
%   «Convergência quadrática»: x^3 - 2x^2 + 4, x0 = 1, tolg = 1e-3 (3 it.).
%   «Newton procura f' = 0»: maximizar 8 + 5x - 3x^4 - 2x^6, x0 = 0.5,
%       tolg = 1e-7 (5 it.).
%   «Os quatro métodos»: o exemplo-guia a partir de x0 = 1 (erro após 5 it.).
%
%   Imprime as tabelas como nos slides, os valores finais e as contagens, e
%   no fim compara com os slides (vírgula decimal nos slides, ponto aqui).
%
%   Otimização — deck 2.5.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.5: Newton-Raphson\n');

% ------------------------------------------------------- exemplo-guia
f   = @(x) 8*exp(1 - x) + 7*log(x);
df  = @(x) -8*exp(1 - x) + 7./x;
ddf = @(x) 8*exp(1 - x) - 7./x.^2;
xs = 1.6093810677230782;                 % minimizante (Newton com muitas casas)
x0 = 1;  tolg = 1e-6;  kmax = 20;
fprintf('\nExemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x0 = %g, tolg = %g\n', x0, tolg);
[x, fx, info] = Newton1D(f, df, ddf, x0, tolg, kmax);
H = info.history;
fprintf('%3s %11s %11s %9s %11s %13s\n', 'k', 'x_k', 'df(x_k)', 'ddf(x_k)', 'x_k+1', '|x_k+1 - x*|');
for r = 1:size(H, 1)
  fprintf('%3d %11.7f %+11.5f %9.4f %11.7f %13.1e\n', H(r, [1 2 4 5 6]), abs(H(r, 6) - xs));
end
fprintf('x* = %.7f, f(x*) = %.7f; %d iterações, n_g = %d, n_H = %d, n_f = %d\n', ...
        x, fx, info.nit, info.ngev, info.nhev, info.nfev);
x1q = 1.5 - df(1.5)/ddf(1.5);
fprintf('Pergunta de abertura: em x0 = 1.5, f''(1.5) = %.4f, f''''(1.5) = %.4f, x1 = %.3f\n', df(1.5), ddf(1.5), x1q);
Tg = [1.0000000 -1.00000 1.0000 2.0000000
      2.0000000  0.55696 1.1930 1.5331535
      1.5331535 -0.12826 1.7160 1.6078977
      1.6078977 -0.00244 1.6484 1.6093804
      1.6093804 -0.0000011 1.6469 1.6093811];
eg = abs(H(:, 6) - xs);
c = [size(H, 1) == 5 && confere(H(:, [2 6]), Tg(:, [1 4]), 7) && confere(H(1:4, 4), Tg(1:4, 2), 5) ...
       && confere(H(5, 4), -1.1e-6, 7) && confere(H(:, 5), Tg(:, 3), 4), ...
     confere(eg', [3.9e-1 7.6e-2 1.5e-3 6.6e-7 1.3e-13], [2 3 4 8 14]), ...
     confere([df(1.5) ddf(1.5) x1q], [-0.1856 1.7411 1.607], [4 4 3]), ...
     info.nit == 5, info.ngev == 6, info.nhev == 5, info.nfev == 0];
fprintf('  tabela k = 0..4: %s | erros: %s | abertura x1 = 1.607: %s | 5 it.: %s | n_g = 6: %s | n_H = 5: %s | n_f = 0: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------- convergência quadrática
f   = @(x) x.^3 - 2*x.^2 + 4;
df  = @(x) 3*x.^2 - 4*x;
ddf = @(x) 6*x - 4;
x0 = 1;  tolg = 1e-3;
fprintf('\nConvergência quadrática: f(x) = x^3 - 2x^2 + 4, x0 = %g, tolg = %g\n', x0, tolg);
[x, ~, info] = Newton1D(f, df, ddf, x0, tolg, 20);
H = info.history;
fprintf('%3s %9s %10s %9s %9s %11s\n', 'k', 'x_k', 'df(x_k)', 'ddf(x_k)', 'x_k+1', 'df(x_k+1)');
for r = 1:size(H, 1)
  fprintf('%3d %9.5f %+10.5f %9.4f %9.5f %+11.5f\n', H(r, [1 2 4 5 6 7]));
end
e = abs([H(:, 2); x] - 4/3);
Hx = ddf(x);                         % classificar: 1 avaliação extra de f''
fprintf('|f''(x_3)| = %.1e < tolg: para em x_3 = %.5f (x* = 4/3); f''''(x_3) = %.4f > 0 (f''''(4/3) = %g): mínimo\n', ...
        abs(H(end, 7)), x, Hx, ddf(4/3));
fprintf('Erros e_k: %.4f -> %.4f -> %.4f -> %.4f;  e_3/e_2^2 = %.2f (C = 0.75)\n', e, e(4)/e(3)^2);
Tc = [1.00000 -1.00000 2.0000 1.50000  0.75000
      1.50000  0.75000 5.0000 1.35000  0.06750
      1.35000  0.06750 4.1000 1.33354  0.00081];
c = [size(H, 1) == 3 && confere(H(:, [2 4 6 7]), Tc(:, [1 2 4 5]), 5) && confere(H(:, 5), Tc(:, 3), 4), ...
     confere(abs(H(end, 7)), 8.1e-4, 5), confere(e', [0.33 0.17 0.017 0.0002], [2 2 3 4]), ...
     ddf(4/3) == 4, Hx > 0, info.ngev == 4, info.nhev == 3];
fprintf('  tabela k = 0..2: %s | |f''(x_3)| = 8.1e-4: %s | erros: %s | f''''(4/3) = 4: %s | mínimo: %s | n_g = 4: %s | n_H = 3: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);
[x03, ~, i03] = Newton1D(f, df, ddf, 0.3, 1e-8, 20);
fprintf('De x0 = 0.3: x = %.6f, f''''(x) = %g < 0: máximo (ponto estacionário errado)\n', x03, ddf(x03));
c = [abs(x03) < 1e-6, ddf(x03) < 0];
fprintf('  converge para 0: %s | máximo: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------ Newton procura f' = 0 (maximizar)
f   = @(x) 8 + 5*x - 3*x.^4 - 2*x.^6;
df  = @(x) 5 - 12*(x.^3 + x.^5);
ddf = @(x) -12*(3*x.^2 + 5*x.^4);
x0 = 0.5;  tolg = 1e-7;
fprintf('\nMaximizar f(x) = 8 + 5x - 3x^4 - 2x^6, x0 = %g, tolg = %g\n', x0, tolg);
[x, fx, info] = Newton1D(f, df, ddf, x0, tolg, 20);
H = info.history;
fprintf('%3s %11s %11s %11s %11s %11s\n', 'k', 'x_k', 'f(x_k)', 'df(x_k)', 'ddf(x_k)', 'x_k+1');
for r = 1:size(H, 1)
  fprintf('%3d %11.7f %11.6f %+11.6f %11.5f %11.7f\n', H(r, 1:6));
end
Hx = ddf(x);                         % classificar: 1 avaliação extra de f''
fprintf('Após %d iterações |f''| < 1e-7: x* = %.7f, f(x*) = %.6f, f''''(x*) = %.4f < 0: máximo\n', ...
        info.nit, x, fx, Hx);
Tm = [0.5000000 10.281250  3.125000 -12.75000 0.7450980
      0.7450980 10.458621 -2.719687 -38.47906 0.6744184
      0.6744184 10.563259 -0.355311 -28.78702 0.6620756
      0.6620756 10.565490 -0.009182 -27.30912 0.6617394
      0.6617394 10.565491 -0.000007 -27.26970 0.6617391];
c = [size(H, 1) == 5 && confere(H(:, [2 6]), Tm(:, [1 5]), 7) && confere(H(:, 3:4), Tm(:, 2:3), 6) ...
       && confere(H(:, 5), Tm(:, 4), 5), ...
     info.nit == 5, confere(x, 0.6617391, 7), confere(fx, 10.565491, 6), Hx < 0];
fprintf('  tabela k = 0..4: %s | 5 it.: %s | x*: %s | f(x*): %s | máximo: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ os quatro métodos
f   = @(x) 8*exp(1 - x) + 7*log(x);
df  = @(x) -8*exp(1 - x) + 7./x;
ddf = @(x) 8*exp(1 - x) - 7./x.^2;
[x5, ~, i5] = Newton1D(f, df, ddf, 1, 1e-6, 20);      % para após 5 it.
fprintf('\nOs quatro métodos, Newton de x0 = 1: erro após %d it. = %.1e (n_g = %d, n_H = %d)\n', ...
        i5.nit, abs(x5 - xs), i5.ngev, i5.nhev);
% o slide escreve 1,3e-13
c = [i5.nit == 5, round(abs(x5 - xs)/1e-14) == 13, i5.ngev == 6, i5.nhev == 5];
fprintf('  5 it.: %s | erro 1,3e-13: %s | n_g = 6: %s | n_H = 5: %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
