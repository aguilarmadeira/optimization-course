% EX04_3_3_CONSTRAINT_HANDLING  Reproduz os números do deck 4.3.3 (restrições na prática).
%
%   Um só problema de referência, resolvido pelo MESMO Hooke–Jeeves da UC (3.2.5),
%   com a restrição comunicada ao algoritmo de cinco formas:
%
%     min f(x) = (x1 - 2)^2 + (x2 - 1)^2   s.a.  g(x) = 2 - x1 - x2 >= 0,
%     x* = (1.5, 0.5), f* = 0.5, u* = 1;  (2, 1) é inadmissível.
%   Convenção do capítulo 4: g >= 0, <g> = min(g, 0); violação v = -<g> = max(0, -g).
%
%     A  ignorar a restrição:         f
%     B  penalização quadrática:      f + R <g>^2 = f + R v^2,  R = 1, 10, 100, 1000
%     C  penalização exata (L1):      f + R |<g>| = f + R v,  R = 0.5, 1, 2
%     D  barreira extrema:            f se g >= 0, senão +Inf (f não é avaliada)
%     E  regra de admissibilidade:    compara [v f], v = max(0, -g); nenhuma soma
%        (HookeJeeves com o argumento opcional menor)
%
%   Hooke–Jeeves: x0 = (0, 0), a = 2, P0 = 0.5, T = 1e-6 (todos os casos).
%   Também: (i) B analítico, v(x_R) = 1/(1 + 2R); (ii) restrição oculta: um
%   «simulador» que falha (NaN) sem dar g, tratado com +Inf = D; (iii) E e D a
%   partir de x0 = (3, 3), inadmissível; (iv) o bloqueio em (1, 1) e a correção
%   com as direções rodadas (HJ nas variáveis y, x = Q y, Q = [1 1; 1 -1]/sqrt 2);
%   (v) informativo: o Nelder–Mead do 3.2.3 em C (R = 2) e D.
%   Funções: ExtremeBarrier.m (D e a restrição oculta), simulador.m (a caixa negra);
%   HookeJeeves (3.2.5) e NelderMead (3.2.3), postas no caminho por uc_setup.
%   Os números dos slides saem daqui; os slides usam vírgula decimal.
%
%   Otimização — deck 4.3.3.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

% --- o problema (o mesmo em todo o deck) -----------------------------------
f = @(x) (x(1) - 2)^2 + (x(2) - 1)^2;
g = @(x) 2 - x(1) - x(2);                                 % g(x) >= 0 (convenção do cap. 4)
v = @(x) max(0, -g(x));                                  % violação v = -<g>
xs = [1.5, 0.5];

% --- os cinco tratamentos: só muda a função (ou a comparação) --------------
FA = f;                                                  % A: ignorar
FB = @(R) @(x) f(x) + R*v(x)^2;                          % B: quadrática
FC = @(R) @(x) f(x) + R*v(x);                            % C: exata (L1)
FD = @(x) ExtremeBarrier(f, g, x);                       % D: barreira extrema
FE = @(x) [v(x), f(x)];                                  % E: o par [v f]
menorE = @(a, b) a(1) < b(1) || (a(1) == b(1) && a(2) < b(2));   % E: regra

x0 = [0, 0];  a = 2;  P0 = 0.5;  T = 1e-6;
HJ = @(F, x0) HookeJeeves(F, x0, a, P0, T, [], [], [], '1961');
HJE = @(x0) HookeJeeves(FE, x0, a, P0, T, [], false, menorE, '1961');

fprintf('Otimização — deck 4.3.3: restrições na prática (um problema, cinco tratamentos)\n');
fprintf('min (x1-2)^2 + (x2-1)^2 s.a. g = 2 - x1 - x2 >= 0;  x* = (1.5, 0.5), f* = 0.5, u* = 1\n');
fprintf('Hooke–Jeeves (3.2.5): x0 = (0, 0), a = 2, P0 = 0.5, T = 1e-6\n\n');
linha = @(nome, x, nf) fprintf(['  %-26s x = (%.6f; %.6f)  f = %.6f  g = %+.2e  ', ...
  'admissível: %-3s  n_f = %4d\n'], nome, x(1), x(2), f(x), g(x), simnao{1 + (g(x) >= -1e-12)}, nf);

% ------------------------------------------------------------ 1. a tabela
fprintf('1. A tabela final (slide «O mesmo problema, cinco tratamentos»)\n');
[xA, ~, iA] = HJ(FA, x0);  linha('A ignorar', xA, iA.nfev);
RB = [1 10 100 1000];  xB = zeros(4, 2);  nB = zeros(1, 4);
for i = 1:4
  [xB(i, :), ~, ii] = HJ(FB(RB(i)), x0);  nB(i) = ii.nfev;
  linha(sprintf('B quadrática, R = %g', RB(i)), xB(i, :), nB(i));
end
RC = [0.5 1 2];  xC = zeros(3, 2);
for i = 1:3
  [xC(i, :), ~, ii] = HJ(FC(RC(i)), x0);
  linha(sprintf('C exata, R = %g', RC(i)), xC(i, :), ii.nfev);
end
[xD, ~, iD] = HJ(FD, x0);  linha('D barreira extrema', xD, iD.nfev);
[xE, ~, iE] = HJE(x0);     linha('E regra [v f]', xE, iE.nfev);
vB = max(0, sum(xB, 2) - 2);                             % v(x_R)
c = [confere(xA, [2 1], 4), ...
     confere(xB, [1.667 0.667; 1.524 0.524; 1.502 0.502; 1.500 0.500], 3), ...
     confere(vB, [0.33 0.048 0.0050 0.00050], [2 3 4 5]), ...
     confere(xC, [1.75 0.75; 1.5 0.5; 1 1], 4), ...
     confere([xD; xE], [1 1; 1 1], 4)];
fprintf(['  A = (2; 1) | B -> x* por fora, v = 0.33 ... 0.0005 | C: (1.75; 0.75), x*, (1; 1)', ...
  ' | D = E = (1; 1): %s %s %s %s %s\n'], simnao{1 + c});
ok = ok && all(c);
nfora = sum(isinf(iD.ftrace));
c = nfora == 45;
fprintf('  D: n_F = %d chamadas, %d fora de X (f não avaliada): %s\n', iD.nfev, nfora, simnao{1 + c});
ok = ok && c;
c = isequal(iD.history(:, 3:4), iE.history(:, 3:4)) && iD.nfev == iE.nfev;
fprintf('  de x0 admissível, E passa pelos mesmos pontos base que D, com o mesmo n_F: %s\n', simnao{1 + c});
fprintf('  (as tentativas em torno de pontos inadmissíveis diferem: D compara +Inf com +Inf, E compara v)\n');
ok = ok && c;

% ------------------------------------------------------------ 2. B analítico
fprintf('\n2. B: x_R = (2, 1) - R v_R (1, 1), v_R = 1/(1 + 2R) (slide «Quanto deve valer R?»)\n');
fprintf('  %6s %22s %10s %10s %12s %6s\n', 'R', 'x_R (HJ)', 'v (HJ)', '1/(1+2R)', '||x_R-x*||', 'n_f');
for i = 1:4
  fprintf('  %6g   (%.4f; %.4f) %10.5f %10.5f %12.5f %6d\n', RB(i), xB(i, 1), xB(i, 2), ...
    vB(i), 1/(1 + 2*RB(i)), norm(xB(i, :) - xs), nB(i));
end
c = confere(vB, 1./(1 + 2*RB), 5);
fprintf('  HJ confere com 1/(1+2R) a 5 casas: %s\n', simnao{1 + c});
ok = ok && c;

% ------------------------------------------------------------ 3. restrição oculta
fprintf('\n3. Restrição oculta: o «simulador» falha (NaN) e não diz porquê (não há g)\n');
FH = @(x) ExtremeBarrier(@simulador, [], x);             % falhou -> +Inf
[xH, ~, iH] = HJ(FH, x0);  linha('oculta, +Inf se falha', xH, iH.nfev);
c = isequaln(iH.history, iD.history);
fprintf('  a mesma corrida que D, sem nunca conhecer g: %s\n', simnao{1 + c});
ok = ok && c;

% ------------------------------------------------------------ 4. x0 inadmissível
fprintf('\n4. De x0 = (3, 3), inadmissível (g = -4)\n');
[xD3, fD3, iD3] = HJ(FD, [3 3]);  linha('D barreira extrema', xD3, iD3.nfev);
[xE3, ~, iE3] = HJE([3 3]);       linha('E regra [v f]', xE3, iE3.nfev);
c = [isinf(fD3) && isequal(xD3, [3 3]), confere(xE3, [1 1], 4)];
fprintf('  D não sai de x0 (tudo +Inf) | E chega à fronteira e para em (1; 1): %s %s\n', simnao{1 + c});
ok = ok && all(c);

% ------------------------------------------------------------ 5. o bloqueio em (1, 1)
fprintf('\n5. Porque para em (1, 1)? As 4 tentativas +-P_j (P = 0.5) a partir de (1, 1):\n');
xb = [1 1];  E4 = [1 0; -1 0; 0 1; 0 -1];  nomes = {'+e1', '-e1', '+e2', '-e2'};
FC2 = FC(2);
for k = 1:4
  xp = xb + 0.5*E4(k, :);
  fprintf('  %s: x = (%.1f; %.1f)  f = %.2f  g = %+.1f  C(R=2) = %.2f  D = %.2f\n', ...
    nomes{k}, xp(1), xp(2), f(xp), g(xp), FC2(xp), FD(xp));
end
fprintf(['  em (1, 1): f = %.2f; nenhuma direção coordenada melhora (qualquer P > 0);', ...
  ' a direção (1, -1)/sqrt 2, ao longo da fronteira, melhora.\n'], f(xb));
c1 = true;
for p = [0.5 1e-3 1e-6]
  for k = 1:4, c1 = c1 && FC2(xb + p*E4(k, :)) > f(xb); end
end
c = [c1, FD(xb + 0.25*[1 -1]) < f(xb)];         % (1.25, 0.75): na fronteira
fprintf('  nenhuma tentativa +-P e_j melhora em C (R = 2), P = 0.5, 1e-3, 1e-6 | (1,-1) melhora: %s %s\n', ...
  simnao{1 + c});
ok = ok && all(c);

% ------------------------------------------------------------ 6. a correção: rodar as direções
fprintf(['\n6. O mesmo Hooke–Jeeves nas variáveis y, x = Q y, Q = [1 1; 1 -1]/sqrt 2', ...
  ' (direções ao longo e na normal à fronteira)\n']);
Q = [1 1; 1 -1]/sqrt(2);
rot = @(F) @(y) F((Q*y(:)).');
y0 = (Q.'*x0(:)).';
[y, ~, i1] = HJ(rot(FC2), y0);  xr(1, :) = (Q*y(:)).';  linha('C exata, R = 2 (rodado)', xr(1, :), i1.nfev);
[y, ~, i2] = HJ(rot(FD), y0);   xr(2, :) = (Q*y(:)).';  linha('D barreira extrema (rodado)', xr(2, :), i2.nfev);
[y, ~, i3] = HookeJeeves(rot(FE), y0, a, P0, T, [], false, menorE, '1961');
xr(3, :) = (Q*y(:)).';  linha('E regra [v f] (rodado)', xr(3, :), i3.nfev);
c = confere(xr, repmat(xs, 3, 1), 4);
fprintf('  com as direções rodadas, C (R = 2), D e E chegam a x* = (1.5; 0.5): %s\n', simnao{1 + c});
ok = ok && c;

% ------------------------------------------------------------ 7. informativo: Nelder–Mead
fprintf(['\n7. Informativo: o Nelder–Mead do 3.2.3 (simplex inicial de x0 = (0, 0), tolx = 1e-8,', ...
  ' tolf = 1e-10) chega a x* aqui, mas sem garantia geral\n']);
[x, ~, in] = NelderMead(FC2, x0, 1e-8, 1e-10, 2000);  linha('C exata, R = 2 (Nelder–Mead)', x, in.nfev);
[x, ~, in] = NelderMead(FD, x0, 1e-8, 1e-10, 2000);   linha('D barreira extrema (Nelder–Mead)', x, in.nfev);

fprintf('\nconfere com os slides: %s\n', simnao{1 + ok});

