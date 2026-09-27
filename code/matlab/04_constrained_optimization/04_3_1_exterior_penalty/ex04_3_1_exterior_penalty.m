% EX04_3_1_EXTERIOR_PENALTY  Reproduz os exemplos do deck 4.3.1 (penalização exterior).
%
%   1. «Ideia»: min x^2 s.a. x - 1 >= 0; x_R = R/(1+R): 1/2, 0.909 (R = 1, 10).
%   2. «O exemplo nas contas» e «As KKT reaparecem»: min x1^2 + x2^2 s.a.
%      x1 + x2 - 1 >= 0, R = 1, 10, 100, 1000: x_R, g(x_R), ||x_R - x*|| = O(1/R), u_R.
%   3. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0:
%      x_R = (1,1) e u_R = 0 para todo o R.
%   4. «A escala das restrições importa»: g1 = 1 - x1, g2 = 100(1 - x2), R = 1, 10, 100.
%   5. «E com igualdades?»: min 2x1^2 + x2^2 s.a. x1 + x2 - 1 = 0: x_R, h(x_R), lambda_R.
%   6. A chamada do slide «No computador» (tolviol = tolx = 1e-6, tmax = 12):
%      com o fminsearch (o do slide; informativo) e com NelderMeadComp.
%   7. A comparação do deck 4.3.2 (lado exterior): de (0,0), 6 ciclos, n_f = 844;
%      de (1,1), 6 ciclos, n_f = 864 (até ||x^(t) - x*|| <= 1e-4, sem a avaliação final).
%
%   Minimizador interno dos exemplos: o Nelder–Mead da UC na versão do script
%   da comparação (NelderMeadComp, tolx = 1e-10, tolf = 1e-12, kmax = 2000;
%   ver help NelderMeadComp e o README). Com ele, os números são iguais aos do
%   exemplo Python. As tabelas 2-5 são uma única corrida com arranque a quente
%   (R^(0) = 1, c = 10), lida em info.history.
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 4.3.1.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;
NM = @(P, x) NelderMeadComp(P, x, 1e-10, 1e-12, 2000);    % interno dos exemplos
sem = @(x) [];                                             % sem restrições deste tipo

fprintf('Otimização — deck 4.3.1: penalização exterior\n');

% ------------------------------------------------------------ 1. Ideia (1D)
[~, ~, info] = PenaltyExterior(@(x) x^2, @(x) x - 1, sem, 0, 1, 10, 0, 0, 2, NM);
xR = info.history(2:end, 6);
fprintf('\n1. min x^2 s.a. x - 1 >= 0: x_R = R/(1+R) (R = 0: x_R = 0, o mínimo livre)\n');
for i = 1:2
  Ri = info.history(i + 1, 2);
  fprintf('   R = %4g: x_R = %.4f  (exato %.4f)\n', Ri, xR(i), Ri/(1 + Ri));
end
c = confere(xR, [0.5; 0.909], 3);
fprintf('  x_R = 1/2, 0.909: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 2. O exemplo nas contas
f = @(x) x(1)^2 + x(2)^2;
g = @(x) x(1) + x(2) - 1;
xs = [0.5 0.5];
[~, ~, info] = PenaltyExterior(f, g, sem, [0 0], 1, 10, 0, 0, 4, NM);
H = info.history(2:end, :);
Rv = H(:, 2);  X = H(:, 6:7);  u = H(:, 8);
gv = sum(X, 2) - 1;  err = sqrt(sum((X - repmat(xs, 4, 1)).^2, 2));
fprintf('\n2. min x1^2 + x2^2 s.a. x1 + x2 - 1 >= 0, de (0,0), R^(0) = 1, c = 10 (x* = (1/2,1/2), u* = 1)\n');
fprintf('  %6s %20s %10s %12s %12s %8s %8s\n', 'R', 'x_R', 'g(x_R)', '||x_R-x*||', 'R*||x_R-x*||', 'u_R', 'n_f');
for i = 1:4
  fprintf('  %6g   (%.4f; %.4f) %10.4f %12.5f %12.4f %8.4f %8d\n', ...
          Rv(i), X(i, 1), X(i, 2), gv(i), err(i), Rv(i)*err(i), u(i), H(i, 3));
end
fprintf('  erro O(1/R): R*||x_R - x*|| -> 1/(2 sqrt 2) = %.4f;  u_R = 2R/(1+2R) -> 1\n', 1/(2*sqrt(2)));
c = [confere(X, repmat([0.333; 0.476; 0.4975; 0.4998], 1, 2), repmat([3; 3; 4; 4], 1, 2)), ...
     confere(gv, [-0.333; -0.0476; -0.0050; -0.0005], [3; 4; 4; 4]), ...
     confere(err, [0.236; 0.0337; 0.0035; 0.00035], [3; 4; 4; 5]), ...
     confere(u, [0.667; 0.952; 0.995; 0.9995], [3; 3; 3; 4])];
fprintf('  x_R: %s | g(x_R): %s | ||x_R - x*||: %s | u_R: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 3. Restrição inativa
[~, ~, info] = PenaltyExterior(@(x) (x(1) - 1)^2 + (x(2) - 1)^2, @(x) 9 - x(1)^2 - x(2)^2, ...
                               sem, [0 0], 1, 10, 0, 0, 4, NM);
H = info.history(2:end, :);
fprintf('\n3. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0\n');
for i = 1:size(H, 1)
  fprintf('   R = %4g: x_R = (%.4f; %.4f), g = %.4f, u_R = %.4f\n', ...
          H(i, 2), H(i, 6), H(i, 7), 9 - H(i, 6)^2 - H(i, 7)^2, H(i, 8) + 0);
end
fprintf('   (%s: x não mudou e a violação é 0)\n', info.message);
c = [confere(H(:, 6:7), ones(size(H, 1), 2), 4), confere(H(:, 8), zeros(size(H, 1), 1), 4)];
fprintf('  x_R = (1,1): %s | u_R = 0: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 4. Escala das restrições
fe = @(x) (x(1) - 2)^2 + (x(2) - 2)^2;
[~, ~, info] = PenaltyExterior(fe, @(x) [1 - x(1); 100*(1 - x(2))], sem, [0 0], 1, 10, 0, 0, 3, NM);
H = info.history(2:end, :);
kap = (2 + 2e4*H(:, 2)) ./ (2 + 2*H(:, 2));    % Hessiana de P em x_R: diag(2+2R, 2+2e4 R)
fprintf('\n4. Escala: min (x1-2)^2 + (x2-2)^2 s.a. 1 - x1 >= 0 e 100(1 - x2) >= 0 (x* = (1,1))\n');
fprintf('  %6s %10s %12s %10s\n', 'R', 'x1_R', 'x2_R', 'kappa');
for i = 1:3
  fprintf('  %6g %10.4f %12.7f %10.0f\n', H(i, 2), H(i, 6), H(i, 7), kap(i));
end
x2 = PenaltyExterior(fe, @(x) [1 - x(1); 1 - x(2)], sem, [0 0], 1, 10, 0, 0, 1, NM);
fprintf('  com g2 = 1 - x2 e R = 1: x_R = (%.4f; %.4f) = ((2+R)/(1+R), (2+R)/(1+R)), kappa = 1\n', x2);
c = [confere(H(:, 6), [1.50; 1.091; 1.0099], [2; 3; 4]), ...
     confere(H(:, 7), [1.0001; 1.00001; 1.000001], [4; 5; 6]), ...
     confere(kap, [5000; 9091; 9901], 0), confere(x2, [1.5 1.5], 4)];
fprintf('  x1_R: %s | x2_R: %s | kappa: %s | com g2: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 5. Igualdades
[~, ~, info] = PenaltyExterior(@(x) 2*x(1)^2 + x(2)^2, sem, @(x) x(1) + x(2) - 1, ...
                               [0 0], 1, 10, 0, 0, 4, NM);
H = info.history(2:end, :);
X = H(:, 6:7);  hv = sum(X, 2) - 1;  lam = H(:, 8);
fprintf('\n5. Igualdade: min 2x1^2 + x2^2 s.a. x1 + x2 - 1 = 0 (x* = (1/3, 2/3), lambda* = 4/3)\n');
for i = 1:4
  fprintf('   R = %4g: x_R = (%.4f; %.4f), h = %.4f, lambda_R = -2Rh = %.3f\n', ...
          H(i, 2), X(i, 1), X(i, 2), hv(i), lam(i));
end
c = [confere(X, [0.2000 0.4000; 0.3125 0.6250; 0.3311 0.6623; 0.3331 0.6662], 4), ...
     confere(hv, [-0.4000; -0.0625; -0.0066; -0.0007], 4), ...
     confere(lam, [0.800; 1.250; 1.325; 1.332], 3)];
fprintf('  x_R: %s | h(x_R): %s | lambda_R: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 6. Chamada do slide
fprintf('\n6. Chamada do slide: PenaltyExterior(f, g, h, [0 0], 0.1, 10, 1e-6, 1e-6, 12)\n');
fprintf('   (a) interno por omissão, o do slide: fminsearch (TolX = TolFun = 1e-10) — informativo\n');
h = @(x) [];
[x, fx, info] = PenaltyExterior(f, g, h, [0 0], 0.1, 10, 1e-6, 1e-6, 12, [], true);
nfc = info.history(2:end, 3);
fprintf('   nit = %d (R = %g), u = %.4f, n_f por ciclo %d-%d, nfev = %d (inclui f(x) final)\n', ...
        info.nit, info.R, info.u(1), min(nfc), max(nfc), info.nfev);
fprintf('   %s\n', info.message);
fprintf('   slide (fminsearch do MATLAB): nit = 10 (R = 1e8); u ~ 1; nfev ~ 130-210 por ciclo\n');
fprintf('   (o fminsearch do Octave não é o do MATLAB e o slide não mostra as opções. Com x_R exatos\n');
fprintf('    o critério vale no ciclo 9: viol = 5e-8 e passo = %.1e <= 1e-6; nit = 10 acontece quando\n', ...
        norm([1 1]*(1e7/(1 + 2e7) - 1e6/(1 + 2e6))));
fprintf('    o erro do minimizador interno ao longo de g = 0 (P mal condicionada) passa 1e-6, como em (b))\n');
fprintf('   (b) com NelderMeadComp (os mesmos números que o exemplo Python):\n');
[x, fx, info] = PenaltyExterior(f, g, h, [0 0], 0.1, 10, 1e-6, 1e-6, 12, NM, true);
nfc = info.history(2:end, 3);
fprintf('   nit = %d (R = %g), u = %.4f, n_f por ciclo %d-%d, nfev = %d (inclui f(x) final)\n', ...
        info.nit, info.R, info.u(1), min(nfc), max(nfc), info.nfev);
fprintf('   %s\n', info.message);
H = info.history;
fprintf('   nota: com x_R exatos o critério valeria no ciclo 9 (viol = 5e-8, passo = %.1e); aqui o passo\n', ...
        sqrt(2)*(1e7/(1 + 2e7) - 1e6/(1 + 2e6)));
fprintf('         no ciclo 9 é %.1e > 1e-6 (erro do interno ao longo de g = 0, P mal condicionada)\n', ...
        norm(H(10, 6:7) - H(9, 6:7))/max(1, norm(H(9, 6:7))));
c = abs(info.u(1) - 1) <= 0.01;
fprintf('  u ~ 1 (a menos de 1 %%): %s | nit = 10: %s (informativo: depende do minimizador interno)\n', ...
        simnao{c + 1}, simnao{(info.nit == 10) + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 7. Comparação (deck 4.3.2)
fprintf('\n7. Comparação exterior vs. interior (deck 4.3.2), lado exterior: R^(0) = 0.1, c = 10,\n');
fprintf('   tolerâncias internas 1e-10; para no 1.º ciclo com ||x^(t) - x*|| <= 1e-4 (lido em history)\n');
x0s = [0 0; 1 1];
res = zeros(2, 3);
for k = 1:2
  [~, ~, info] = PenaltyExterior(f, g, h, x0s(k, :), 0.1, 10, 0, 0, 12, NM);
  H = info.history;
  e = sqrt(sum((H(2:end, 6:7) - repmat(xs, size(H, 1) - 1, 1)).^2, 2));
  t = find(e <= 1e-4, 1);
  res(k, :) = [t, H(t + 1, 4), floor(median(H(2:t + 1, 3)))];
  fprintf('   de (%g,%g):\n', x0s(k, :));
  for i = 2:t + 1
    fprintf('     t = %d  R = %8.1e  x = (%.5f, %.5f)  ||x - x*|| = %.1e  n_f = %3d  acum. = %d\n', ...
            H(i, 1), H(i, 2), H(i, 6), H(i, 7), e(i - 1), H(i, 3), H(i, 4));
  end
  fprintf('     -> %d ciclos, n_f = %d, mediana por ciclo = %d\n', res(k, :));
end
c = [isequal(res(1, 1:2), [6 844]), res(1, 3) == 140, isequal(res(2, 1:2), [6 864])];
fprintf('  (0,0): 6 ciclos, n_f = 844: %s | mediana 140 (~140 por ciclo): %s | (1,1): 6 ciclos, n_f = 864: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% informativo: o mesmo com o NelderMead do deck 3.2.3 (contração exterior aceite com <=)
NM3 = @(P, x) NelderMead(P, x, 1e-10, 1e-12, 2000);
v = zeros(2, 2);
for k = 1:2
  [~, ~, info] = PenaltyExterior(f, g, h, x0s(k, :), 0.1, 10, 0, 0, 12, NM3);
  H = info.history;
  e = sqrt(sum((H(2:end, 6:7) - repmat(xs, size(H, 1) - 1, 1)).^2, 2));
  t = find(e <= 1e-4, 1);  v(k, :) = [t, H(t + 1, 4)];
end
fprintf('   informativo — com NelderMead do deck 3.2.3 (aceita a contração exterior com <=):\n');
fprintf('   (0,0): %d ciclos, n_f = %d;  (1,1): %d ciclos, n_f = %d\n', v(1, :), v(2, :));

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
