% EX02_3_POWELL  Reproduz os exemplos do deck 2.3 (interpolação quadrática).
%
%   Exemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x1 = 1, Delta = 0.5,
%            tolx = tolf = 1e-3: 4 iterações (k = 0..3), sem reflexões,
%            n_f = 3 + 4 = 7, x* ~ 1.6094 (exato 1.6093811).
%   «Quando a parábola não serve»: f(x) = x^4/4 - x^2/2 com os trios
%            (0.3, 0.45, 0.6) -> a2 < 0, xb = -0.46;
%            (0.55, 0.6, 0.65) -> a2 ~ 0+, xb = 5.31, f(xb) = 184;
%            e a salvaguarda (reflexão) a partir do primeiro trio.
%
%   Imprime a tabela das iterações como nos slides, os valores finais e as
%   contagens, e no fim compara com os slides (vírgula decimal nos slides,
%   ponto aqui).
%
%   Otimização — deck 2.3.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.3: interpolação quadrática sucessiva (Powell)\n');

% ------------------------------------------------------------ exemplo-guia
f = @(x) 8*exp(1 - x) + 7*log(x);
xs = 1.6093810677230782;                 % x* (f'(x*) = 0)
x1 = 1;  Delta = 0.5;  tolx = 1e-3;  tolf = 1e-3;  kmax = 50;
fprintf('\nExemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x1 = %g, Delta = %g, tolx = tolf = %g\n', ...
        x1, Delta, tolx);
[x, fx, info] = PowellQuadratic(f, x1, Delta, tolx, tolf, kmax, true);
fprintf('Fim: x* = %.4f, f(x*) = %.4f; %d iterações, %d reflexões, n_f = 3 + %d = %d\n', ...
        x, fx, info.nit, info.nrefl, info.nit + info.nrefl, info.nfev);
fprintf('Erros |xb - x*|:');  fprintf(' %.1e', abs(info.history(:, 10) - xs));  fprintf('\n');

% tabela do slide: x1 x2 x3 (4 casas) f1 f2 f3 (4) a1 a2 (3) xb f(xb) (4)
T = [1.0000 1.5000 2.0000 8.0000 7.6905 7.7951 -0.619 0.828 1.6237 7.6806
     1.5000 1.6237 2.0000 7.6905 7.6806 7.7951 -0.080 0.768 1.6139 7.6805
     1.5000 1.6139 1.6237 7.6905 7.6805 7.6806 -0.088 0.837 1.6096 7.6804
     1.5000 1.6096 1.6139 7.6905 7.6804 7.6805 -0.092 0.840 1.6094 7.6804];
H = info.history;
tabok = size(H, 1) == 4 && confere(H(:, [2:7, 10:11]), T(:, [1:6, 9:10]), 4) ...
        && confere(H(:, 8:9), T(:, 7:8), 3);
erros = abs(H(:, 10) - xs);
c = [tabok, info.nit == 4, info.nrefl == 0, info.nfev == 7, ...
     confere(x, 1.6094, 4), confere(fx, 7.6804, 4), ...
     confere(erros', [0.014 0.0045 0.0002 5e-5], [3 4 4 5]), all(H(2:end, 2) == 1.5)];
fprintf('  tabela k = 0..3: %s | 4 it.: %s | sem reflexões: %s | n_f = 7: %s | x* = 1.6094: %s | f = 7.6804: %s | erros: %s | 1.5 fica no trio: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------ quando a parábola não serve
g = @(x) x.^4/4 - x.^2/2;
fprintf('\nQuando a parábola não serve: f(x) = x^4/4 - x^2/2 (mínimos em -1 e 1)\n');
P = [0.3 0.45 0.6;  0.55 0.6 0.65];
xbP = zeros(2, 1);  a2P = zeros(2, 1);
for r = 1:2
  X = P(r, :);  F = g(X);
  a1 = (F(2)-F(1)) / (X(2)-X(1));
  a2 = ((F(3)-F(1))/(X(3)-X(1)) - a1) / (X(3)-X(2));
  xbP(r) = (X(1)+X(2))/2 - a1/(2*a2);  a2P(r) = a2;
  fprintf('  trio (%.2f, %.2f, %.2f): a1 = %.4f, a2 = %.4f, xb = %.2f, f(xb) = %.1f\n', ...
          X, a1, a2, xbP(r), g(xbP(r)));
end
[xr, fr, inforr] = PowellQuadratic(g, 0.3, 0.15, 1e-3, 1e-3, 50);
fprintf('  com a salvaguarda, de x1 = 0.3, Delta = 0.15: x* = %.4f, f = %.4f, %d reflexão(ões), %d it., n_f = %d\n', ...
        xr, fr, inforr.nrefl, inforr.nit, inforr.nfev);
c = [a2P(1) < 0, confere(xbP(1), -0.46, 2), a2P(2) > 0, confere(xbP(2), 5.31, 2), ...
     confere(g(xbP(2)), 184, 0), inforr.flag == 0, confere(xr, 1, 2)];
fprintf('  a2 < 0: %s | xb = -0.46: %s | a2 > 0: %s | xb = 5.31: %s | f(xb) = 184: %s | reflexão converge: %s | x* ~ 1: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
