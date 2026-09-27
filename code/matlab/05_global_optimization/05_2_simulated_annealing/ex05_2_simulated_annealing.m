% EX05_2_SIMULATED_ANNEALING  Reproduz os exemplos do deck 5.2 (simulated annealing).
%
%   Função de referência do cap. 5: f(x) = 0.1 x^2 + 3 sin 2x + 2 cos(3.3x + 1)
%   em X = [-8, 8], mínimo global x* = 2.4681 (f = -4.2381).
%
%   1. «Uma execução passo a passo»: x0 = 6, T0 = 3, c = 0.98, sigma = 0.5,
%      Nmax = 400: T em k = 30, 100, 200 e n_f = 401; e «SA -> método local»:
%      de x_SA = 2.464, uma pesquisa direta chega a |x - x*| < 1e-6 com 37
%      avaliações.
%   2. «Uma execução não chega: 100 corridas»: a tabela de 6 linhas (sucesso
%      |x - x*| < 0.3, sempre de x0 = 6): 35, 19, 74, 20, 88, 73 %.
%   3. «Em duas dimensões: Rastrigin»: T0 = 20, c = 0.995, sigma = 0.5,
%      Nmax = 3000, n_f = 3001; a chamada do slide «No computador» (30
%      corridas, rng(s), s = 1..30), sucesso f < 1: 100 % no slide.
%      Com completo = true, também a variação de c (97/97/97/100/90 %).
%
%   ATENÇÃO — verificação estatística. Os números dos slides vêm do script
%   Python das figuras (gerador numpy, sementes 3, 99 e 11). O MATLAB/Octave
%   tem outro gerador: com rng(s) as corridas são reprodutíveis, mas não são
%   as do slide. Por isso aqui:
%     - os valores determinísticos (T_k, T_Nmax, n_f, a pesquisa direta final)
%       conferem-se exatamente;
%     - as taxas de sucesso conferem-se estatisticamente. A taxa do slide é
%       também uma estimativa (em Rs corridas); a diferença entre as duas
%       estimativas tem de ficar abaixo de 3 desvios-padrão binomiais,
%       3*sqrt(p(1-p)*(1/R + 1/Rs)), com p(1-p) >= 1/R (para p = 0 ou 1).
%       Com R = Rs = 100 a tolerância é larga (cerca de 20 pontos
%       percentuais para p = 0.35): verifica a ordem de grandeza, não o
%       valor exato;
%     - a corrida isolada (passo a passo) é outra corrida: mostra-se, mas o
%       seu resultado (escapar ou não) não entra na verificação.
%   O exemplo Python (mesmo gerador e sementes do script) reproduz os valores
%   exatos.
%
%   Demora cerca de 25 s em Octave (completo = true: cerca de 50 s).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 5.2.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

completo = false;     % true: acrescenta a variação de c no Rastrigin (5 x 30 corridas)

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
% taxa p (em [0,1], R corridas) compatível com a taxa ps do slide (Rs corridas)?
dentro = @(p, ps, R, Rs) abs(p - ps) <= 3*sqrt(max(ps*(1 - ps), 1/R)*(1/R + 1/Rs));
simnao = {'não', 'sim'};
ok = true;

f1 = @(x) 0.1*x.^2 + 3*sin(2*x) + 2*cos(3.3*x + 1);
LB = -8;  UB = 8;
XS = 2.4680679685;               % mínimo global da função de referência

fprintf('Otimização — deck 5.2: simulated annealing\n');
fprintf('(verificação estatística: outro gerador que não o do script Python dos slides)\n');

% ------------------------------------------------ 1. Execução passo a passo
x0 = 6;  T0 = 3;  c = 0.98;  sigma = 0.5;  Nmax = 400;
fprintf('\n1. Execução passo a passo: x0 = %g, T0 = %g, c = %g, sigma = %g, Nmax = %d, rng(3)\n', ...
        x0, T0, c, sigma, Nmax);
rng(3);
[x, fx, info] = SimulatedAnnealing(f1, x0, T0, c, sigma, Nmax, LB, UB);
H = info.history;                % [k T aceite f(x') f(x_k) f_best x' x_k]
fprintf('%5s %9s %9s %9s %9s\n', 'k', 'T', 'x_k', 'f(x_k)', 'f_best');
for k = [0 30 100 200 400]
  fprintf('%5d %9.4f %9.4f %9.4f %9.4f\n', k, H(k+1, 2), H(k+1, 8), H(k+1, 5), H(k+1, 6));
end
fprintf('x = %.4f, f(x) = %.4f;  n_f = %d, aceites = %d, T final = %.2e\n', ...
        x, fx, info.nfev, info.nacc, info.T);
txt = {'não terminou', 'terminou'};
fprintf('  (outra corrida que não a do slide: aqui %s na bacia do global; no slide, sim, em k = 315)\n', ...
        txt{(abs(x - XS) < 0.3) + 1});
c1 = [confere(H(31, 2), 1.7, 1), confere(H(101, 2), 0.4, 1), confere(H(201, 2), 0.05, 2), ...
      info.nfev == 401];
fprintf('  T(30) = 1.7: %s | T(100) = 0.4: %s | T(200) = 0.05: %s | n_f = 401: %s\n', simnao{c1 + 1});
ok = ok && all(c1);

% SA -> método local: pesquisa direta 1D (poll {+1, -1}), passo inicial 0.01,
% mantido no sucesso e reduzido a metade no insucesso, até ser <= 1e-6; não
% conta f(x_SA). Parte do x_SA do slide (corrida Python, semente 3).
x = 2.46421979;  fx = f1(x);  nf = 0;  a = 0.01;
while a > 1e-6
  sucesso = false;
  for d = [1 -1]
    y = x + a*d;  fy = f1(y);  nf = nf + 1;
    if fy < fx, x = y; fx = fy; sucesso = true; break; end
  end
  if ~sucesso, a = a/2; end
end
fprintf('SA -> método local: de 2.4642, pesquisa direta: x = %.7f, |x - x*| = %.1e, %d avaliações\n', ...
        x, abs(x - XS), nf);
c1 = [abs(x - XS) < 1e-6, nf == 37];
fprintf('  |x - x*| < 1e-6: %s | 37 avaliações: %s\n', simnao{c1 + 1});
ok = ok && all(c1);

% ------------------------------------------------ 2. 100 corridas
R = 100;
fprintf('\n2. %d corridas de x0 = 6 (rng(99)), sucesso |x - x*| < 0.3\n', R);
fprintf('%6s %6s %7s %6s %10s %8s %8s %10s\n', 'sigma', 'T0', 'c', 'Nmax', 'T_Nmax', ...
        'sucesso', 'slide', '3 d.p.');
linhas = [0.5 3.0 0.980  400 35
          0.5 3.0 0.900  400 19
          0.5 3.0 0.995  400 74
          0.5 0.3 0.980  400 20
          1.0 3.0 0.980  400 88
          0.5 3.0 0.980 2000 73];
Tslide = [9.3e-4 1.5e-18 4.0e-1 9.3e-5 9.3e-4 8.5e-18];
cT = true(1, 6);  cp = true(1, 6);
for i = 1:6
  sg = linhas(i, 1);  T0_ = linhas(i, 2);  c_ = linhas(i, 3);  N_ = linhas(i, 4);
  ps = linhas(i, 5)/100;
  rng(99);
  suc = 0;
  for r = 1:R
    [xr, ~, inf_] = SimulatedAnnealing(f1, 6, T0_, c_, sg, N_, LB, UB);
    suc = suc + (abs(xr - XS) < 0.3);
  end
  p = suc/R;
  cT(i) = strcmp(sprintf('%.1e', inf_.T), sprintf('%.1e', Tslide(i)));
  cp(i) = dentro(p, ps, R, 100);
  fprintf('%6.1f %6.1f %7.3f %6d %10.1e %7.0f%% %7.0f%% %9s\n', sg, T0_, c_, N_, inf_.T, ...
          100*p, 100*ps, simnao{cp(i) + 1});
end
fprintf('  T_Nmax: %s | taxas dentro de 3 d.p. das do slide (35/19/74/20/88/73 %%): %s\n', ...
        simnao{all(cT) + 1}, simnao{all(cp) + 1});
ok = ok && all(cT) && all(cp);

% ------------------------------------------------ 3. Rastrigin 2D
f = @(x) 20 + sum(x.^2 - 10*cos(2*pi*x));
fprintf('\n3. Rastrigin 2D em [-5.12, 5.12]^2: T0 = 20, c = 0.995, sigma = 0.5, Nmax = 3000\n');
rng(3);
[x, fx, info] = SimulatedAnnealing(f, [4 -4], 20, 0.995, 0.5, 3000, -5.12, 5.12);
fprintf('de (4, -4), rng(3): x = (%.4f; %.4f), f = %.4f, |x| = %.4f, n_f = %d, aceites = %d\n', ...
        x(1), x(2), fx, norm(x), info.nfev, info.nacc);
fprintf('  (slide, outra corrida: f = 0.068, |x| = 0.02)\n');
c3 = info.nfev == 3001;
fprintf('  n_f = 3001: %s\n', simnao{c3 + 1});
ok = ok && c3;

% a chamada do slide «No computador»
fs = zeros(1, 30);
for s = 1:30                    % 30 corridas
    rng(s);  x0 = -5.12 + 10.24*rand(1,2);
    [x, fx, info] = SimulatedAnnealing(f, x0, ...
              20, 0.995, 0.5, 3000, -5.12, 5.12);
    fs(s) = fx;                 % info.nfev = 3001
end
% Q1 e Q3 por interpolação linear (como o percentile do numpy); evita
% prctile, que no MATLAB é da Statistics and Machine Learning Toolbox
v = sort(fs);  h = 1 + (numel(v) - 1)*[0.25 0.75];  l = floor(h);
q = v(l) + (h - l).*(v(min(l + 1, numel(v))) - v(l));
fprintf('chamada do slide (rng(s), s = 1..30): mediana f = %.3f, [Q1; Q3] = [%.3f; %.3f], sucesso f < 1: %.0f%%\n', ...
        median(fs), q(1), q(2), 100*mean(fs < 1));
c3 = [dentro(mean(fs < 1), 1.00, 30, 30), info.nfev == 3001];
fprintf('  sucesso compatível com 100 %% (>= 86 %%): %s | info.nfev = 3001: %s\n', simnao{c3 + 1});
ok = ok && all(c3);

if completo
  fprintf('30 corridas de pontos aleatórios (rng(11)), sucesso f < 1:\n');
  cc = [0.90 97; 0.95 97; 0.99 97; 0.995 100; 0.999 90];
  c3 = true(1, 5);
  for i = 1:5
    rng(11);
    suc = 0;
    for r = 1:30
      x0 = -5.12 + 10.24*rand(1, 2);
      [~, fx] = SimulatedAnnealing(f, x0, 20, cc(i, 1), 0.5, 3000, -5.12, 5.12);
      suc = suc + (fx < 1);
    end
    c3(i) = dentro(suc/30, cc(i, 2)/100, 30, 30);
    fprintf('  c = %.3f: sucesso %3.0f%%  (slide %d%%)  dentro de 3 d.p.: %s\n', ...
            cc(i, 1), 100*suc/30, cc(i, 2), simnao{c3(i) + 1});
  end
  ok = ok && all(c3);
else
  fprintf('(variação de c: 97/97/97/100/90 %% — só com completo = true)\n');
end

fprintf('\nconfere com os slides (valores determinísticos exatos; taxas: estatisticamente): %s\n', ...
        simnao{ok + 1});
