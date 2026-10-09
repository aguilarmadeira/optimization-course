% EX05_2_SIMULATED_ANNEALING  Reproduz os exemplos do deck 5.2 (simulated annealing, versão das aulas).
%
%   Função de referência do cap. 5: f(x) = 0.1 x^2 + 3 sin 2x + 2 cos(3.3x + 1)
%   em X = [-8, 8], mínimo global x* = 2.4681 (f = -4.2381).
%
%   1. «Uma execução passo a passo»: x0 = 6, T0 = 3, alpha = 0.8, n = 10,
%      Ts = 0.05, sigma = 0.5: termina com t = 191 aceites e
%      T = 3*0.8^19 = 0.043 < Ts; e «SA -> método local»: de x_t = 2.457,
%      uma pesquisa direta chega a |x - x*| < 1e-6 com 38 avaliações.
%   2. «Uma execução não chega: 100 corridas»: a tabela de 6 linhas (sempre
%      de x0 = 6): sucesso |x - x*| < 0.3 do ponto final e do melhor
%      visitado, e a mediana de n_f.
%   3. «Em duas dimensões: Rastrigin»: T0 = 20, alpha = 0.5, n = 10, Ts = 0.1,
%      sigma = 0.5 (81 aceites, T final = 20*0.5^8 = 0.078); a chamada do
%      slide «No computador» (30 corridas, rng(s), s = 1..30), sucesso f < 1:
%      90 % no slide. Com completo = true, também alpha = 0.3 e 0.8.
%
%   ATENÇÃO — verificação estatística. Os números dos slides vêm do script
%   Python das figuras (gerador numpy, sementes 3, 99 e 11). O MATLAB/Octave
%   tem outro gerador: com rng(s) as corridas são reprodutíveis, mas não são
%   as do slide. Por isso aqui:
%     - os valores determinísticos (aceites e temperatura final, a pesquisa
%       direta final) conferem-se exatamente;
%     - as taxas de sucesso conferem-se estatisticamente. A taxa do slide é
%       também uma estimativa (em Rs corridas); a diferença entre as duas
%       estimativas tem de ficar abaixo de 3 desvios-padrão binomiais,
%       3*sqrt(p(1-p)*(1/R + 1/Rs)), com p(1-p) >= 1/R (para p = 0 ou 1).
%       Com R = Rs = 100 a tolerância é larga (cerca de 20 pontos
%       percentuais para p = 0.6): verifica a ordem de grandeza, não o
%       valor exato; a mediana de n_f tem de cair no intervalo
%       interquartil [Q1; Q3] do slide;
%     - a corrida isolada (passo a passo) é outra corrida: mostra-se, mas o
%       seu resultado (escapar ou não) não entra na verificação.
%   O exemplo Python (mesmo gerador e sementes do script) reproduz os valores
%   exatos.
%
%   Demora cerca de um minuto em Octave (completo = true: mais um minuto).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 5.2.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

completo = false;     % true: acrescenta alpha = 0.3 e 0.8 no Rastrigin (2 x 30 corridas)

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
% taxa p (em [0,1], R corridas) compatível com a taxa ps do slide (Rs corridas)?
dentro = @(p, ps, R, Rs) abs(p - ps) <= 3*sqrt(max(ps*(1 - ps), 1/R)*(1/R + 1/Rs));
simnao = {'não', 'sim'};
ok = true;

f1 = @(x) 0.1*x.^2 + 3*sin(2*x) + 2*cos(3.3*x + 1);
LB = -8;  UB = 8;
XS = 2.4680679685;               % mínimo global da função de referência

fprintf('Otimização — deck 5.2: simulated annealing (versão das aulas)\n');
fprintf('(verificação estatística: outro gerador que não o do script Python dos slides)\n');

% ------------------------------------------------ 1. Execução passo a passo
x0 = 6;  T0 = 3;  alpha = 0.8;  n = 10;  Ts = 0.05;  sigma = 0.5;
fprintf('\n1. Execução passo a passo: x0 = %g, T0 = %g, alpha = %g, n = %d, Ts = %g, sigma = %g, rng(3)\n', ...
        x0, T0, alpha, n, Ts, sigma);
rng(3);
[x, fx, info] = SimulatedAnnealing(f1, x0, T0, alpha, n, Ts, sigma, LB, UB);
fprintf('ponto final x_t = %.4f (f = %.4f), melhor visitado %.4f;  n_f = %d, aceites = %d, T final = %.4f\n', ...
        x, fx, info.xbest, info.nfev, info.nacc, info.T);
txt = {'não terminou', 'terminou'};
fprintf('  (outra corrida que não a do slide: aqui %s na bacia do global; no slide, sim, com o salto em n_f = 325)\n', ...
        txt{(abs(x - XS) < 0.3) + 1});
c1 = [info.flag == 0, info.nacc == 191, confere(info.T, 0.043, 3), confere(3*0.8^19, 0.043, 3)];
fprintf('  terminou por T < Ts: %s | 191 aceites: %s | T final = 0.043: %s | 3*0.8^19 = 0.043: %s\n', simnao{c1 + 1});
ok = ok && all(c1);

% SA -> método local: pesquisa direta 1D (poll {+1, -1}), passo inicial 0.01,
% mantido no sucesso e reduzido a metade no insucesso, até ser <= 1e-6; não
% conta f(x_SA). Parte do x_t do slide (corrida Python, semente 3).
x = 2.45737074;  fx = f1(x);  nf = 0;  a = 0.01;
while a > 1e-6
  sucesso = false;
  for d = [1 -1]
    y = x + a*d;  fy = f1(y);  nf = nf + 1;
    if fy < fx, x = y; fx = fy; sucesso = true; break; end
  end
  if ~sucesso, a = a/2; end
end
fprintf('SA -> método local: de 2.4574, pesquisa direta: x = %.7f, |x - x*| = %.1e, %d avaliações\n', ...
        x, abs(x - XS), nf);
c1 = [abs(x - XS) < 1e-6, nf == 38];
fprintf('  |x - x*| < 1e-6: %s | 38 avaliações: %s\n', simnao{c1 + 1});
ok = ok && all(c1);

% ------------------------------------------------ 2. 100 corridas
R = 100;
fprintf('\n2. %d corridas de x0 = 6 (rng(99)), Ts = 0.05, sucesso |x - x*| < 0.3\n', R);
fprintf('%6s %5s %6s %4s %9s %9s %8s %8s %16s %8s\n', 'sigma', 'T0', 'alpha', 'n', 'n_f med', 'slide', ...
        'final', 'melhor', 'slide', '3 d.p.');
linhas = [0.5 3.0 0.80 10  693 60 63  525  828
          0.5 3.0 0.50 10  148 11 11  136  168
          0.5 3.0 0.95 10 3318 69 93 3187 3470
          0.5 0.3 0.80 10  404 15 15  368  440
          1.0 3.0 0.80 10 1518 60 92 1445 1578
          0.5 3.0 0.80 50 3942 64 93 3669 4125];       % ... mediana, sucessos, Q1 e Q3 de n_f do slide
cp = true(1, 6);  cn = true(1, 6);
for i = 1:6
  sg = linhas(i, 1);  T0_ = linhas(i, 2);  a_ = linhas(i, 3);  n_ = linhas(i, 4);
  rng(99);
  suc = 0;  sucb = 0;  nfs = zeros(1, R);
  for r = 1:R
    [xr, ~, inf_] = SimulatedAnnealing(f1, 6, T0_, a_, n_, 0.05, sg, LB, UB);
    suc = suc + (abs(xr - XS) < 0.3);  sucb = sucb + (abs(inf_.xbest - XS) < 0.3);  nfs(r) = inf_.nfev;
  end
  cp(i) = dentro(suc/R, linhas(i, 6)/100, R, 100) && dentro(sucb/R, linhas(i, 7)/100, R, 100);
  cn(i) = median(nfs) >= linhas(i, 8) && median(nfs) <= linhas(i, 9);
  fprintf('%6.1f %5.1f %6.2f %4d %9.1f %9d %7.0f%% %7.0f%%   %3d%% / %3d%% %8s\n', sg, T0_, a_, n_, ...
          median(nfs), linhas(i, 5), suc, sucb, linhas(i, 6), linhas(i, 7), simnao{(cp(i) && cn(i)) + 1});
end
fprintf('  taxas dentro de 3 d.p. das do slide: %s | medianas de n_f dentro de [Q1; Q3] do slide: %s\n', ...
        simnao{all(cp) + 1}, simnao{all(cn) + 1});
ok = ok && all(cp) && all(cn);

% ------------------------------------------------ 3. Rastrigin 2D
f = @(x) 20 + sum(x.^2 - 10*cos(2*pi*x));
fprintf('\n3. Rastrigin 2D em [-5.12, 5.12]^2: T0 = 20, alpha = 0.5, n = 10, Ts = 0.1, sigma = 0.5\n');
rng(3);
[x, fx, info] = SimulatedAnnealing(f, [4 -4], 20, 0.5, 10, 0.1, 0.5, -5.12, 5.12);
fprintf('de (4, -4), rng(3): x_t = (%.4f; %.4f), f = %.4f, n_f = %d, aceites = %d, T final = %.4f\n', ...
        x(1), x(2), fx, info.nfev, info.nacc, info.T);
fprintf('  (slide, outra corrida: f = 0.18, n_f = 5421)\n');
c3 = [info.nacc == 81, confere(info.T, 0.078, 3)];
fprintf('  81 aceites: %s | T final = 20*0.5^8 = 0.078: %s\n', simnao{c3 + 1});
ok = ok && all(c3);

% a chamada do slide «No computador»
fs = zeros(1, 30);  nf = zeros(1, 30);
for s = 1:30                    % 30 corridas
    rng(s);  x0 = -5.12 + 10.24*rand(1,2);
    [x, fx, info] = SimulatedAnnealing(f, x0, ...
              20, 0.5, 10, 0.1, 0.5, -5.12, 5.12);
    fs(s) = fx;  nf(s) = info.nfev;
end
% Q1 e Q3 por interpolação linear (como o percentile do numpy); evita
% prctile, que no MATLAB é da Statistics and Machine Learning Toolbox
v = sort(fs);  h = 1 + (numel(v) - 1)*[0.25 0.75];  l = floor(h);
q = v(l) + (h - l).*(v(min(l + 1, numel(v))) - v(l));
fprintf('chamada do slide (rng(s), s = 1..30): mediana f = %.3f, [Q1; Q3] = [%.3f; %.3f], sucesso f < 1: %.0f%%, n_f mediana %.0f\n', ...
        median(fs), q(1), q(2), 100*mean(fs < 1), median(nf));
c3 = [dentro(mean(fs < 1), 0.90, 30, 30), median(nf) >= 4166 && median(nf) <= 5521];
fprintf('  sucesso compatível com 90 %%: %s | n_f mediana em [4166; 5521]: %s\n', simnao{c3 + 1});
ok = ok && all(c3);

if completo
  fprintf('30 corridas de pontos aleatórios (rng(11)), sucesso f < 1 (ponto final):\n');
  cc = [0.3 97 3630 3279 4395; 0.8 97 18572 17010 19879];
  c3 = true(1, 2);
  for i = 1:2
    rng(11);
    suc = 0;  nf = zeros(1, 30);
    for r = 1:30
      x0 = -5.12 + 10.24*rand(1, 2);
      [~, fx, info] = SimulatedAnnealing(f, x0, 20, cc(i, 1), 10, 0.1, 0.5, -5.12, 5.12);
      suc = suc + (fx < 1);  nf(r) = info.nfev;
    end
    c3(i) = dentro(suc/30, cc(i, 2)/100, 30, 30);     % com 30 corridas, a mediana de n_f é só informativa
    fprintf('  alpha = %.1f: sucesso %3.0f%%, n_f mediana %.0f  (slide %d%%, %d)  sucesso compatível: %s\n', ...
            cc(i, 1), 100*suc/30, median(nf), cc(i, 2), cc(i, 3), simnao{c3(i) + 1});
  end
  ok = ok && all(c3);
else
  fprintf('(alpha = 0.3 e 0.8: 97 %% com n_f ~ 3600 e ~ 18 600 — só com completo = true)\n');
end

fprintf('\nconfere com os slides (valores determinísticos exatos; taxas: estatisticamente): %s\n', ...
        simnao{ok + 1});
