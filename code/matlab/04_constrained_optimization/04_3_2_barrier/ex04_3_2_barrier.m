% EX04_3_2_BARRIER  Reproduz os exemplos do deck 4.3.2 (barreira logarítmica).
%
%   1. «Ideia»: min x^2 s.a. x - 1 >= 0; x_R = (1 + sqrt(1+2R))/2: 1.366, 1.048, 1.005
%      (R = 1, 0.1, 0.01).
%   2. «O exemplo nas contas» e «As KKT reaparecem --- perturbadas»: min x1^2 + x2^2
%      s.a. x1 + x2 - 1 >= 0, R = 1, 0.1, 0.01, 0.001: x_R, g(x_R), ||x_R - x*|| = O(R), u_R.
%   3. Barreira inversa R/g, R = 0.01: x_R = (0.548; 0.548), contra (0.5050; 0.5050) da log.
%   4. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0:
%      R = 1 -> (0.882; 0.882), R = 0.1 -> (0.986; 0.986); u_R = R/g -> 0.
%   5. A chamada do slide «No computador» (tolcomp = 1e-8, tolx = 1e-6, tmax = 12):
%      com o fminsearch (o do slide; informativo) e com NelderMeadComp.
%   6. A comparação exterior vs. interior (slide «Exterior vs. interior»): barreira de
%      (1,1) em 5 ciclos, 743 chamadas a P, 26 fora de X, n_f = 717; exterior de (0,0)
%      em 6 ciclos, n_f = 844, e de (1,1), n_f = 864 (até ||x^(t) - x*|| <= 1e-4, sem a
%      avaliação final). Usa PenaltyExterior da pasta 04_3_1_exterior_penalty.
%
%   Minimizador interno dos exemplos: o Nelder–Mead da UC na versão do script
%   da comparação (NelderMeadComp, tolx = 1e-10, tolf = 1e-12, kmax = 2000;
%   ver help NelderMeadComp e o README). Com ele, os números são iguais aos do
%   exemplo Python. As tabelas são uma única corrida com arranque a quente
%   (c = 0.1), lida em info.history.
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 4.3.2.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;
NM = @(P, x) NelderMeadComp(P, x, 1e-10, 1e-12, 2000);    % interno dos exemplos

fprintf('Otimização — deck 4.3.2: barreira logarítmica\n');

% ------------------------------------------------------------ 1. Ideia (1D)
[~, ~, info] = BarrierLog(@(x) x^2, @(x) x - 1, 2, 1, 0.1, 0, 0, 3, NM);
xR = info.history(2:end, 6);
fprintf('\n1. min x^2 s.a. x - 1 >= 0, de x = 2: x_R = (1 + sqrt(1+2R))/2\n');
for i = 1:3
  Ri = info.history(i + 1, 2);
  fprintf('   R = %5g: x_R = %.4f  (exato %.4f)\n', Ri, xR(i), (1 + sqrt(1 + 2*Ri))/2);
end
c = confere(xR, [1.366; 1.048; 1.005], 3);
fprintf('  x_R = 1.366, 1.048, 1.005: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 2. O exemplo nas contas
f = @(x) x(1)^2 + x(2)^2;
g = @(x) x(1) + x(2) - 1;
xs = [0.5 0.5];
[~, ~, info] = BarrierLog(f, g, [1 1], 1, 0.1, 0, 0, 4, NM);
H = info.history(2:end, :);
Rv = H(:, 2);  X = H(:, 6:7);  u = H(:, 8);
gv = sum(X, 2) - 1;  err = sqrt(sum((X - repmat(xs, 4, 1)).^2, 2));
fprintf('\n2. min x1^2 + x2^2 s.a. x1 + x2 - 1 >= 0, de (1,1), R^(0) = 1, c = 0.1 (x* = (1/2,1/2), u* = 1)\n');
fprintf('  %6s %20s %10s %12s %12s %8s %6s %6s\n', 'R', 'x_R', 'g(x_R)', '||x_R-x*||', ...
        '||x_R-x*||/R', 'u_R', 'P', 'fora');
for i = 1:4
  fprintf('  %6g   (%.4f; %.4f) %10.4f %12.5f %12.4f %8.4f %6d %6d\n', ...
          Rv(i), X(i, 1), X(i, 2), gv(i), err(i), err(i)/Rv(i), u(i), H(i, 3), H(i, 4));
end
fprintf(['  erro O(R): ||x_R - x*||/R -> 1/sqrt 2 = %.4f;  u_R = R/g = (1 + sqrt(1+4R))/2 -> 1;', ...
         '  u_R g(x_R) = R\n'], 1/sqrt(2));
c = [confere(X, repmat([0.809; 0.5458; 0.5050; 0.5005], 1, 2), repmat([3; 4; 4; 4], 1, 2)), ...
     confere(gv, [0.618; 0.0916; 0.0099; 0.0010], [3; 4; 4; 4]), ...
     confere(err, [0.437; 0.0648; 0.0070; 0.00071], [3; 4; 4; 5]), ...
     confere(u, [1.618; 1.092; 1.010; 1.001], 3), all(gv > 0)];
fprintf('  x_R: %s | g(x_R): %s | ||x_R - x*||: %s | u_R: %s | todos admissíveis: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 3. Barreira inversa
R = 0.01;
Pinv = @(x) f(x) + R/max(g(x), 0);            % R/0 = Inf fora de X
xi = NelderMeadComp(Pinv, [1 1], 1e-10, 1e-12, 2000);
fprintf('\n3. Barreira inversa P = f + R/g, R = 0.01, de (1,1): x_R = (%.4f; %.4f), ||x_R - x*|| = %.4f\n', ...
        xi(1), xi(2), norm(xi - xs));
fprintf('   log (tabela acima), R = 0.01: x_R = (%.4f; %.4f), ||x_R - x*|| = %.4f  -> inversa O(sqrt R), log O(R)\n', ...
        X(3, 1), X(3, 2), err(3));
c = [confere(xi, [0.548 0.548], 3), confere(X(3, :), [0.505 0.505], 3)];
fprintf('  inversa (0.548; 0.548): %s | log (0.505; 0.505): %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 4. Restrição inativa
[~, ~, info] = BarrierLog(@(x) (x(1) - 1)^2 + (x(2) - 1)^2, @(x) 9 - x(1)^2 - x(2)^2, ...
                          [0 0], 1, 0.1, 0, 0, 2, NM);
H = info.history(2:end, :);
fprintf('\n4. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0, de (0,0) (x* = (1,1), u* = 0)\n');
for i = 1:2
  fprintf('   R = %4g: x_R = (%.4f; %.4f), g = %.4f, u_R = R/g = %.4f\n', ...
          H(i, 2), H(i, 6), H(i, 7), 9 - H(i, 6)^2 - H(i, 7)^2, H(i, 8));
end
c = [confere(H(:, 6:7), [0.882 0.882; 0.986 0.986], 3), H(2, 8) < H(1, 8) && H(1, 8) < 1];
fprintf('  x_R = (0.882; 0.882), (0.986; 0.986): %s | u_R a descer para 0: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 5. Chamada do slide
fprintf('\n5. Chamada do slide: BarrierLog(f, g, [1 1], 1, 0.1, 1e-8, 1e-6, 12)\n');
fprintf('   (a) interno por omissão, o do slide: fminsearch (TolX = TolFun = 1e-10) — informativo\n');
[x, fx, info] = BarrierLog(f, g, [1 1], 1, 0.1, 1e-8, 1e-6, 12, [], true);
nfc = info.history(2:end, 3);
fprintf('   nit = %d (R = %.3g), u = %.4f, chamadas a P por ciclo %d-%d, nfev = chamadas a P + 1 = %d\n', ...
        info.nit, info.R, info.u(1), min(nfc), max(nfc), info.nfev);
fprintf('   fora de X: %d; n_f efetivo = %d (nfev é majorante de n_f); %s\n', ...
        info.nfora, info.nfev_efetivo, info.message);
fprintf('   slide (fminsearch do MATLAB): nit = 10 (R = 1e-9), u ~ 1; nfev ~ 140-190 por ciclo\n');
fprintf('   (o fminsearch do Octave não é o do MATLAB; ver o README)\n');
fprintf('   (b) com NelderMeadComp (os mesmos números que o exemplo Python):\n');
[x, fx, info] = BarrierLog(f, g, [1 1], 1, 0.1, 1e-8, 1e-6, 12, NM, true);
nfc = info.history(2:end, 3);
fprintf('   nit = %d (R = %.3g), u = %.4f, chamadas a P por ciclo %d-%d, nfev = chamadas a P + 1 = %d\n', ...
        info.nit, info.R, info.u(1), min(nfc), max(nfc), info.nfev);
fprintf('   fora de X: %d; n_f efetivo = %d (nfev é majorante de n_f); %s\n', ...
        info.nfora, info.nfev_efetivo, info.message);
fprintf('   nota: R no ciclo 9 = 0.1^8 calculado como %.17g > 1e-8 = tolcomp: J*R <= tolcomp\n', ...
        info.history(10, 2));
fprintf('         só pode valer a partir do ciclo 10 (R ~ 1e-9), qualquer que seja o minimizador interno\n');
c = abs(info.u(1) - 1) <= 0.01;
fprintf('  u ~ 1 (a menos de 1 %%): %s | nit = 10: %s (informativo: depende do minimizador interno)\n', ...
        simnao{c + 1}, simnao{(info.nit == 10) + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 6. Comparação
fprintf('\n6. Comparação exterior vs. interior: o mesmo Nelder–Mead interno (tolerâncias 1e-10),\n');
fprintf('   arranque a quente; para no 1.º ciclo com ||x^(t) - x*|| <= 1e-4 (lido em history)\n');
ate = @(H) find(sqrt(sum((H(2:end, 6:7) - repmat(xs, size(H, 1) - 1, 1)).^2, 2)) <= 1e-4, 1);
[~, ~, info] = BarrierLog(f, g, [1 1], 1, 0.1, 0, 0, 12, NM);
Hb = info.history;  tb = ate(Hb);
nPb = Hb(tb + 1, 5);  forab = sum(Hb(2:tb + 1, 4));
fprintf('   barreira de (1,1), R^(0) = 1, c = 0.1:\n');
for i = 2:tb + 1
  fprintf('     t = %d  R = %8.1e  x = (%.5f, %.5f)  ||x - x*|| = %.1e  P = %3d (fora %2d)  acum. = %d\n', ...
          Hb(i, 1), Hb(i, 2), Hb(i, 6), Hb(i, 7), norm(Hb(i, 6:7) - xs), Hb(i, 3), Hb(i, 4), Hb(i, 5));
end
fprintf('     -> %d ciclos, %d chamadas a P, %d fora de X, n_f = %d (sem a avaliação final)\n', ...
        tb, nPb, forab, nPb - forab);
x0s = [0 0; 1 1];  ext = zeros(2, 3);
for k = 1:2
  [~, ~, info] = PenaltyExterior(f, g, @(x) [], x0s(k, :), 0.1, 10, 0, 0, 12, NM);
  He = info.history;  te = ate(He);
  ext(k, :) = [te, He(te + 1, 4), floor(median(He(2:te + 1, 3)))];
end
fprintf('   exterior de (0,0), R^(0) = 0.1, c = 10: %d ciclos, n_f = %d\n', ext(1, 1:2));
fprintf('   exterior de (1,1), R^(0) = 0.1, c = 10: %d ciclos, n_f = %d\n', ext(2, 1:2));
medb = floor(median(Hb(2:tb + 1, 3)));
fprintf('\n                                                 exterior  barreira\n');
fprintf('   ciclos exteriores t                          %9d %9d\n', ext(1, 1), tb);
fprintf('   por ciclo (mediana; n_f | chamadas a P)      %9d %9d\n', ext(1, 3), medb);
fprintf('   total até ||x^(t) - x*|| <= 1e-4             %9d %9d\n', ext(1, 2), nPb);
c = [isequal([tb nPb forab nPb - forab], [5 743 26 717]), isequal(ext(1, 1:2), [6 844]), ...
     isequal(ext(2, 1:2), [6 864]), isequal([ext(1, 3) medb], [140 146])];
fprintf(['  barreira 5 ciclos, 743 P, 26 fora, n_f = 717: %s | exterior (0,0) 6 ciclos, 844: %s |', ...
         ' exterior (1,1) 6 ciclos, 864: %s | medianas 140 e 146: %s\n'], simnao{c + 1});
ok = ok && all(c);

% informativo: o mesmo com o NelderMead do deck 3.2.3 (contração exterior aceite com <=)
[~, ~, info] = BarrierLog(f, g, [1 1], 1, 0.1, 0, 0, 12, @(P, x) NelderMead(P, x, 1e-10, 1e-12, 2000));
H3 = info.history;  t3 = ate(H3);
fprintf('   informativo — com NelderMead do deck 3.2.3 (aceita a contração exterior com <=):\n');
fprintf('   barreira de (1,1): %d ciclos, %d chamadas a P, %d fora de X, n_f = %d\n', ...
        t3, H3(t3 + 1, 5), sum(H3(2:t3 + 1, 4)), H3(t3 + 1, 5) - sum(H3(2:t3 + 1, 4)));

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
