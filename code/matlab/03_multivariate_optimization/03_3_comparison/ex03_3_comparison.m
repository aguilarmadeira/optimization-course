% EX03_3_COMPARISON  Reproduz a tabela «Os métodos estudados no mesmo
% problema, com o custo em avaliações» do fim do capítulo 3 (deck 3.3.3),
% com as funções da UC.
%
%   Regras (folha de convenções da UC):
%   * mesmo problema (Rosenbrock, fg.m), mesmo ponto inicial x0 = (-1.5; 2),
%     mesmo critério: custo até à PRIMEIRA avaliação de f com f < 1e-4 (f nos
%     iterandos dos métodos com derivadas também conta em n_f);
%   * custo em avaliações contadas à parte, n_f, n_g, n_H (contador.m);
%     equivalente em f com derivadas centrais (n = 2): n_f + 4 n_g + 4 n_H
%     (a Hessiana, 2n^2+1 = 9 pontos, reaproveita os 4 do gradiente e f(x_k));
%   * métodos de 3.3 com a pesquisa em linha precisa (Brent, tol 1e-10);
%   * Newton amortecido: o n_f até f < 1e-4 depende do último bit da solução
%     de H*d = -g (a pesquisa em linha de Brent com tol 1e-10 é sensível a
%     isso): o numpy (LU) dá 281, como os slides; H\g do Octave (Cholesky)
%     dá 280. Aqui aceita-se n_f a menos de 2 % (n_g, n_H exatos).
%   * pesquisa aleatória localizada (estocástica): 30 corridas, rng(0..29).
%     O gerador do MATLAB/Octave não é o do numpy usado nos slides, por isso
%     aqui a verificação é ESTATÍSTICA: mediana a menos de 15 % de 975 e
%     sucesso (f < 1e-4) nas 30 corridas. A versão Python reproduz
%     exatamente 975 [867; 1052].
%
%   Valores dos slides (n_f / n_g / n_H / equiv.):
%     aleatória localizada 975 [867; 1052]; Nelder-Mead 166; Box 13 445;
%     Hooke-Jeeves 353; gradiente 60 900 / 3000 / - / 72 900 (não atinge
%     f < 1e-4 em 3000 it.); FR (reinício n) 718 / 34 / - / 854;
%     Newton puro 6 / 5 / 5 / 46; Newton amortecido 281 / 12 / 12 / 377.
%
%   Usa as funções das pastas 03_2_*, 03_3_* e common/ (postas no caminho
%   por uc_setup). Em Octave demora cerca de 20 s.
%
%   Otimização — deck 3.3.3 (comparação do cap. 3).  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

simnao = {'não', 'sim'};
equiv = @(nf, ng, nH) nf + 4*ng + 4*nH;        % n = 2: 2n = 4, 2n^2 - 2n = 4
pad = @(s, w) [s, blanks(w - sum(s < 128 | s >= 192))];   % alinhar texto com acentos (UTF-8)
F = @(x) contador('f', x);
G = @(x) contador('g', x);
H = @(x) contador('H', x);
X0 = [-1.5; 2];
FTOL = 1e-4;
ok = true;

fprintf('Otimização — capítulo 3: os métodos estudados no mesmo problema (Rosenbrock de (-1.5; 2))\n');
fprintf('Critério: primeira avaliação com f < %g; custo equivalente em f = n_f + 4 n_g + 4 n_H\n', FTOL);

% ------------------------------------------ pesquisa aleatória localizada
nhit = zeros(30, 1);
for s = 0:29
  contador('inicio', FTOL);
  LocalRandomSearch(F, X0', 1.0, 20, 0.9, 1e-8, 5000, -Inf, s);
  C = contador('estado');
  if isempty(C.hit), nhit(s + 1) = Inf; else, nhit(s + 1) = C.hit(1); end
end
fin = sort(nhit(isfinite(nhit)));
med = median(fin);
q = interp1(0:numel(fin) - 1, fin, [0.25 0.75]*(numel(fin) - 1));   % quartis como no numpy
nfalhas = sum(~isfinite(nhit));
fprintf('\n%s %6s %18s %6s %5s %8s\n', pad('Método', 30), 'ordem', 'n_f', 'n_g', 'n_H', 'equiv.');
fprintf('%s %6d %18s %6s %5s %8.0f\n', pad('Aleatória localizada (30x)', 30), 0, ...
        sprintf('%.0f [%.0f; %.0f]', med, q(1), q(2)), '--', '--', med);
linhas_ok = abs(med - 975) <= 0.15*975 && nfalhas == 0;

% ------------------------------------------------ os outros seis métodos
nomes = {'Nelder-Mead', 'Box (Delta_0 = 1)', 'Hooke-Jeeves (P_0 = 0.5)', 'Gradiente', ...
         'Grad. conj. (FR, reinício n)', 'Newton puro', 'Newton amortecido'};
ordem = [0 0 0 1 1 2 2];
slides = [166 0 0 166; 13445 0 0 13445; 353 0 0 353; 60900 3000 0 72900; ...
          718 34 0 854; 6 5 5 46; 281 12 12 377];
res = zeros(7, 4);
for j = 1:7
  contador('inicio', FTOL);
  switch j
    case 1, [x, fx, info] = NelderMead(F, [X0'; X0' + [0.5 0]; X0' + [0 0.5]], 1e-8, 1e-10, 1000);
    case 2, [x, fx, info] = BoxEvo(F, X0', [1 1], 1e-6, 10000);
    case 3, [x, fx, info] = HookeJeeves(F, X0', 2, 0.5, 1e-6, 100000);
    case 4, [x, fx, info] = SteepestDescent(F, G, X0, 1e-8, 3000);
    case 5, [x, fx, info] = ConjGrad(F, G, X0, 1e-8, 3000);
    case 6, [x, fx, info] = NewtonND(F, G, H, X0, 1e-8, 200);
    case 7, [x, fx, info] = NewtonND(F, G, H, X0, 1e-8, 200, struct('damped', true, 'modify', true));
  end
  C = contador('estado');
  % as contagens do contador coincidem com as da própria função?
  if info.nfev ~= C.nf || info.ngev ~= C.ng || info.nhev ~= C.nH
    fprintf('  AVISO: %s conta (%d, %d, %d), o contador (%d, %d, %d)\n', nomes{j}, ...
            info.nfev, info.ngev, info.nhev, C.nf, C.ng, C.nH);
    ok = false;
  end
  if isempty(C.hit)                            % não atinge: custo total da corrida
    h = [C.nf, C.ng, C.nH];  marca = '‡';
    if j == 4, fgrad = fx; end
  else
    h = C.hit;  marca = '';
  end
  res(j, :) = [h, equiv(h(1), h(2), h(3))];
  sg = '--'; if h(2) > 0, sg = sprintf('%d', h(2)); end
  sH = '--'; if h(3) > 0, sH = sprintf('%d', h(3)); end
  fprintf('%s %6d %18d %6s %5s %8d\n', pad([nomes{j} marca], 30), ordem(j), h(1), sg, sH, res(j, 4));
  if j == 7
    % Newton amortecido: n_f até f < 1e-4 é sensível ao último bit de H\g
    % (numpy/LU dá 281; Octave, H\g por Cholesky, 280): tolerância de 2 %
    linhas_ok(end + 1) = isequal(res(j, 2:3), slides(j, 2:3)) && abs(res(j, 1) - slides(j, 1)) <= 0.02*slides(j, 1); %#ok<SAGROW>
  else
    linhas_ok(end + 1) = isequal(res(j, :), slides(j, :)); %#ok<SAGROW>
  end
end
fprintf('‡ não atinge f < 1e-4 em 3000 iterações (f = %.1e); custo total da corrida\n', fgrad);

nd = res(7, :);  fr = res(5, :);
fprintf('\nDeck 3.3.3: até f < 1e-4, Newton amortecido n_g = n_H = %d, n_f = %d na pesquisa em linha (%d com f nos iterandos);\n', ...
        nd(2), nd(1) - nd(2), nd(1));
fprintf('            FR n_g = %d, n_f = %d na pesquisa em linha (%d).\n', fr(2), fr(1) - fr(2), fr(1));
fprintf('            (a primeira f < 1e-4 surge numa pesquisa em linha; até aí os iterandos avaliados são tantos quantos os gradientes)\n');
c2 = [nd(2) == 12 && abs(nd(1) - 281) <= 0.02*281, fr(2) == 34 && fr(1) - fr(2) == 684];
fprintf('  Newton: 269 + 12 = 281 (a menos de 2 %%): %s | FR: 684 + 34 = 718: %s\n', simnao{c2 + 1});

fprintf(['\n  aleatória: mediana a menos de 15%% de 975 e 30/30 (verificação estatística): %s | NM 166: %s', ...
         ' | Box 13 445: %s | HJ 353: %s | gradiente 72 900: %s | FR 854: %s | Newton puro 46: %s', ...
         ' | amortecido 377 (n_f a menos de 2 %%): %s\n'], simnao{linhas_ok + 1});
fprintf(['  Nota: o n_f do Newton amortecido até f < 1e-4 depende do arredondamento de H\\g (Brent com tol 1e-10\n', ...
         '  num ponto quase ótimo): numpy (LU) dá 281, como nos slides; H\\g do Octave (Cholesky) %d.\n'], nd(1));
ok = ok && all(linhas_ok) && all(c2);

% ------------------------------------------------ o mesmo com o fminunc (informativo)
if exist('fminunc') > 0 %#ok<EXIST>
  try
    opts = optimset('GradObj', 'on', 'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'off');
    [xu, fu, ~, out] = fminunc(@fg, X0, opts);
    fprintf('\nO mesmo com o fminunc (informativo: quasi-Newton, critério e pesquisa em linha do solver): f = %.1e, funcCount = %d, iterations = %d\n', ...
            fu, out.funcCount, out.iterations);
  catch err
    fprintf('\n(fminunc não correu: %s)\n', err.message);
  end
end

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
