% EX05_3_GENETIC_ALGORITHMS  Reproduz os exemplos do deck 5.3 (algoritmos genéticos).
%
%   1. «Uma geração à mão»: min x^2 em [-5, 5], 8 bits, N = 6, roleta com
%      F = 1/(1+x^2), cruzamento a 1 ponto (pc = 0.8), mutação por bit
%      (pm = 0.05). Aqui parte-se dos números sorteados que o slide mostra
%      (população inicial, xi da roleta, cortes 3, 4, 7 e a mutação do bit 6
%      do descendente 2) e refazem-se as contas: x, f, F, p_i, acumulada,
%      pais, descendentes e as médias de f: 11.9 -> 10.4 (o melhor piora:
%      0.111 -> 0.323). É o AG binário do slide, não GeneticAlgorithm (que
%      usa codificação real); por isso está escrito aqui, passo a passo.
%   2. «Rastrigin»: N = 40, G = 60, pc = 0.9, pm = 0.1, sigma = 0.5, 2 elites:
%      n_f = 2440; a chamada do slide «No computador» (30 corridas, rng(s),
%      s = 1..30), sucesso f < 1: 100 % no slide.
%   3. «O que fazem os parâmetros?»: função de referência do cap. 5, sucesso
%      |x - x*| < 0.3: referência 75 %, s = 6 61 %, pc = 0 60 %, N = 6 34 %,
%      N = 60 96 %; pm = 0 76 %, pm = 0.8 78 %; mesmo orçamento n_f = 420:
%      N = 60 (G = 6) 96 %, N = 6 (G = 69) 46 %. O slide usa 300 corridas;
%      aqui, por omissão, 100 (completo = false); com completo = true, 300.
%
%   ATENÇÃO — verificação estatística. Os números dos slides vêm dos scripts
%   Python das figuras (gerador numpy, sementes 5, 2, 21 e 11). O MATLAB/Octave
%   tem outro gerador: com rng(s) as corridas são reprodutíveis, mas não são
%   as do slide. Por isso aqui:
%     - a geração à mão e as contagens n_f = N (G + 1) conferem-se exatamente;
%     - as taxas de sucesso conferem-se estatisticamente. A taxa do slide é
%       também uma estimativa (300 corridas); a diferença entre as duas
%       estimativas tem de ficar abaixo de 3 desvios-padrão binomiais,
%       3*sqrt(p(1-p)*(1/R + 1/Rs)), com p(1-p) >= 1/R (para p = 0 ou 1).
%       Com R = 100 a tolerância é cerca de 15 pontos percentuais para
%       p = 0.75; com R = 300 (completo = true), cerca de 10.6.
%       Nota: em 3000 corridas (Python; Octave semelhante) as taxas «de
%       longo prazo» de pc = 0, N = 6 e N = 6 com G = 69 são cerca de 54, 39
%       e 40 %; os 60, 34 e 46 % do slide são estimativas de 300 corridas a
%       cerca de 2 d.p. (ver o README).
%   O exemplo Python (mesmo gerador e sementes dos scripts) reproduz os
%   valores exatos.
%
%   Demora cerca de 25 s em Octave (completo = true: cerca de 1 min).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 5.3.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

completo = false;     % true: 300 corridas por linha na parte 3, como no slide

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
% taxa p (em [0,1], R corridas) compatível com a taxa ps do slide (Rs corridas)?
dentro = @(p, ps, R, Rs) abs(p - ps) <= 3*sqrt(max(ps*(1 - ps), 1/R)*(1/R + 1/Rs));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 5.3: algoritmos genéticos\n');
fprintf('(verificação estatística das taxas: outro gerador que não o dos scripts Python dos slides)\n');

% ------------------------------------------------ 1. Uma geração à mão
fprintf('\n1. Uma geração à mão: min x^2 em [-5, 5], 8 bits, N = 6, roleta, pc = 0.8, pm = 0.05\n');
L = -5;  U = 5;  nb = 8;  N = 6;
decode = @(b) L + (U - L)*bin2dec(b)/(2^nb - 1);      % b: cadeias de '0'/'1' (linhas)
pop = ['11010110'; '10001000'; '00010110'; '00010111'; '00101110'; '11100100'];  % sorteada (slide)
x = decode(pop);  fx = x.^2;
Fa = 1./(1 + fx);  p = Fa/sum(Fa);  acum = cumsum(p);
fprintf('%2s %9s %7s %7s %6s %6s %6s\n', 'i', 'cromos.', 'x', 'f', 'F', 'p_i', 'acum.');
for i = 1:N
  fprintf('%2d %9s %7.3f %7.3f %6.3f %6.3f %6.3f\n', i, pop(i, :), x(i), fx(i), Fa(i), p(i), acum(i));
end
xi = [0.679 0.870 0.227 0.895 0.872 0.019];            % sorteados (slide)
pais = zeros(1, N);
for k = 1:N, pais(k) = find(acum >= xi(k), 1); end     % roleta: 1.º i com acumulada >= xi
fprintf('roleta: xi =%s -> pais%s\n', sprintf(' %.3f', xi), sprintf(' %d', pais));
cortes = [3 4 7];                                      % sorteados (slide); os 3 pares cruzam
filhos = repmat(' ', N, nb);
for j = 1:N/2
  a = pop(pais(2*j - 1), :);  b = pop(pais(2*j), :);  cj = cortes(j);
  filhos(2*j - 1, :) = [a(1:cj), b(cj+1:end)];         % cruzamento a 1 ponto
  filhos(2*j, :)     = [b(1:cj), a(cj+1:end)];
  fprintf('pais %d, %d: %s, %s  corte após bit %d  -> %s, %s\n', pais(2*j - 1), pais(2*j), ...
          a, b, cj, filhos(2*j - 1, :), filhos(2*j, :));
end
filhos_m = filhos;                                     % mutação sorteada (slide):
filhos_m(2, 6) = char('0' + '1' - filhos(2, 6));       % descendente 2, bit 6
x1 = decode(filhos_m);  fx1 = x1.^2;
for i = 1:N
  fprintf('%2d %9s -> %9s %7.3f %7.3f\n', i, filhos(i, :), filhos_m(i, :), x1(i), fx1(i));
end
fprintf('média de f: %.1f -> %.1f;  melhor: %.3f -> %.3f\n', mean(fx), mean(fx1), min(fx), min(fx1));
c1 = [confere(x, [3.392 0.333 -4.137 -4.098 -3.196 3.941], 3), ...
      confere(fx, [11.507 0.111 17.117 16.794 10.215 15.533], 3), ...
      confere(Fa, [0.080 0.900 0.055 0.056 0.089 0.060], 3), ...
      confere(p, [0.064 0.725 0.044 0.045 0.072 0.049], 3), ...
      confere(acum, [0.064 0.790 0.834 0.879 0.951 1.000], 3), ...
      isequal(pais, [2 4 2 5 4 1]), ...
      isequal(cellstr(filhos_m)', {'10010111', '00001100', '10001110', '00101000', '00010110', '11010111'}), ...
      confere(x1, [0.922 -4.529 0.569 -3.431 -4.137 3.431], 3), ...
      confere(fx1, [0.849 20.516 0.323 11.774 17.117 11.774], 3), ...
      confere([mean(fx) mean(fx1)], [11.9 10.4], 1) && confere([min(fx) min(fx1)], [0.111 0.323], 3)];
fprintf(['  x: %s | f: %s | F: %s | p_i: %s | acumulada: %s | pais 2 4 2 5 4 1: %s |\n', ...
         '  descendentes: %s | x: %s | f: %s | médias 11.9 -> 10.4: %s\n'], simnao{c1 + 1});
ok = ok && all(c1);

% ------------------------------------------------ 2. Rastrigin 2D
f = @(x) 20 + sum(x.^2 - 10*cos(2*pi*x));
fprintf('\n2. Rastrigin 2D em [-5.12, 5.12]^2: N = 40, G = 60, pc = 0.9, pm = 0.1, sigma = 0.5, 2 elites\n');
rng(2);
[x, fx, info] = GeneticAlgorithm(f, 2, -5.12, 5.12, 40, 60, 0.9, 0.1, 0.5, 2);
H = info.history;                                      % [t melhor média x]
for t = [0 5 20 60]
  fprintf('  geração %2d: melhor f = %.3g, média = %.1f\n', t, H(t + 1, 2), H(t + 1, 3));
end
fprintf('rng(2): x = (%.1e; %.1e), f = %.2e, n_f = %d  (slide, outra corrida: f < 1e-12)\n', ...
        x(1), x(2), fx, info.nfev);
c2 = info.nfev == 2440;
fprintf('  n_f = 2440: %s\n', simnao{c2 + 1});
ok = ok && c2;

% a chamada do slide «No computador»
fs = zeros(1, 30);
for s = 1:30
  rng(s);
  [x, fx, info] = GeneticAlgorithm(f, ...
     2, -5.12, 5.12, 40, 60, 0.9, 0.1, 0.5, 2);
  fs(s) = fx;              % info.nfev = 2440
end
% Q1 e Q3 por interpolação linear (como o percentile do numpy); evita
% prctile, que no MATLAB é da Statistics and Machine Learning Toolbox
v = sort(fs);  h = 1 + (numel(v) - 1)*[0.25 0.75];  l = floor(h);
q = v(l) + (h - l).*(v(min(l + 1, numel(v))) - v(l));
fprintf('chamada do slide (rng(s), s = 1..30): mediana f = %.1e, [Q1; Q3] = [%.1e; %.1e], sucesso f < 1: %.0f%%\n', ...
        median(fs), q(1), q(2), 100*mean(fs < 1));
c2 = [dentro(mean(fs < 1), 1.00, 30, 30), info.nfev == 2440];
fprintf('  sucesso compatível com 100 %% (>= 86 %%): %s | info.nfev = 2440: %s\n', simnao{c2 + 1});
ok = ok && all(c2);

% ------------------------------------------------ 3. O que fazem os parâmetros?
f1 = @(x) 0.1*x.^2 + 3*sin(2*x) + 2*cos(3.3*x + 1);
LB = -8;  UB = 8;  XS = 2.4680679685;
if completo, R = 300; else, R = 100; end
fprintf('\n3. Função de referência, %d corridas (rng(11) em cada linha), sucesso |x - x*| < 0.3\n', R);
%          N   G   pc   pm   s  slide(%) n_f(slide)
linhas = [20  20  0.9  0.1  2  75   420
          20  20  0.9  0.1  6  61   420
          20  20  0.0  0.1  2  60   420
           6  20  0.9  0.1  2  34   126
          60  20  0.9  0.1  2  96  1260
          20  20  0.9  0.0  2  76   420
          20  20  0.9  0.8  2  78   420
          60   6  0.9  0.1  2  96   420
           6  69  0.9  0.1  2  46   420];
nomes = {'referência: N = 20, s = 2, pm = 0.1', 'pressão forte: s = 6', 'sem cruzamento: pc = 0', ...
         'população pequena: N = 6', 'população maior: N = 60', 'sem mutação: pm = 0', ...
         'mutação alta: pm = 0.8', 'mesmo orçamento: N = 60, G = 6', 'mesmo orçamento: N = 6, G = 69'};
c3t = true(1, 9);  c3n = true(1, 9);
for i = 1:9
  Ni = linhas(i, 1);  Gi = linhas(i, 2);  ps = linhas(i, 6)/100;
  rng(11);
  suc = 0;
  for r = 1:R
    [xr, ~, inf_] = GeneticAlgorithm(f1, 1, LB, UB, Ni, Gi, linhas(i, 3), linhas(i, 4), 0.5, 2, linhas(i, 5));
    suc = suc + (abs(xr - XS) < 0.3);
  end
  c3t(i) = dentro(suc/R, ps, R, 300);
  c3n(i) = inf_.nfev == linhas(i, 7);
  fprintf('  sucesso %5.1f%%  (slide %2d%%, 3 d.p.: %s)  n_f = %4d (slide %4d)  %s\n', ...
          100*suc/R, linhas(i, 6), simnao{c3t(i) + 1}, inf_.nfev, linhas(i, 7), nomes{i});
end
fprintf('  taxas dentro de 3 d.p. das do slide: %s | n_f = N(G+1): %s\n', simnao{all(c3t) + 1}, ...
        simnao{all(c3n) + 1});
ok = ok && all(c3t) && all(c3n);

fprintf('\nconfere com os slides (valores determinísticos exatos; taxas: estatisticamente, %d corridas): %s\n', ...
        R, simnao{ok + 1});
