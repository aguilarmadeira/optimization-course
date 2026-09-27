% EX03_2_4_BOX  Reproduz os exemplos do deck 3.2.4 (método de Box, EVOP simplificado).
%
%   «Uma iteração à mão»: Himmelblau em x = (0,0), Delta = (2,2): f nos 4 vértices.
%   Himmelblau a partir de (0,0), Delta = (2,2), tolx = 1e-4: tabela k = 0..9,
%     20 iterações, n_f = 4*20 + 1 = 81, termina em (3,2) com f = 0.
%   Rosenbrock a partir de (-1.5, 2), Delta_0 = (1,1): n_f = 13 445 até f < 1e-4
%     (tolx = 1e-6, como no caderno cap3_comparacao.ipynb; o deck não indica tolx);
%     a caixa tem de encolher até Delta ~ 6e-5.
%   «Porque falha em vales curvos?»: em (0,0), passo 0.25, as 4 diagonais dão
%     f entre 4.08 e 11.33 (> 1), mas f(0.25, 0) = 0.953.
%
%   Contagens: f(x0) conta e em cada iteração avaliam-se SEMPRE os 2^n vértices
%   (nota do deck: «nas contagens deste deck avaliam-se sempre os 2^n vértices»).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 3.2.4.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

himmel = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
rosen  = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;

fprintf('Otimização — deck 3.2.4: método de Box (EVOP simplificado)\n');

% ------------------------------------------------------------ Himmelblau
x0 = [0 0];  Delta = [2 2];  tolx = 1e-4;
fprintf('\nHimmelblau a partir de (%g,%g), Delta = (%g,%g), tolx = %g\n', x0, Delta, tolx);
[x, fx, info] = BoxEvo(himmel, x0, Delta, tolx, [], true);
fprintf('%s\n', info.message);
fprintf('x = (%.4f, %.4f), f(x) = %.4g;  %d iterações, n_f = %d (= 4 x %d + 1)\n', ...
        x, fx, info.nit, info.nfev, info.nit);

% «Uma iteração à mão»: vértices (-1,-1), (1,-1), (-1,1), (1,1), pela ordem do slide
% (a ordem de avaliação de BoxEvo é (-1,-1), (-1,1), (1,-1), (1,1))
fv = info.ftrace(1 + [1 3 2 4]);
fprintf('Iteração à mão: f(0,0) = %g; vértices (-1,-1), (1,-1), (-1,1), (1,1): %g, %g, %g, %g\n', ...
        info.ftrace(1), fv);

% tabela do slide (k = 0..9): x_k(1) x_k(2)  Delta  f(x_k)  decisão (1 = move)  f novo
T = [0.0 0.0 2.000 170.000 1 106.000
     1.0 1.0 2.000 106.000 1  26.000
     2.0 2.0 2.000  26.000 1  10.000
     3.0 1.0 2.000  10.000 0  10.000
     3.0 1.0 1.000  10.000 1   9.125
     3.5 1.5 1.000   9.125 1   0.000
     3.0 2.0 1.000   0.000 0   0.000
     3.0 2.0 0.500   0.000 0   0.000
     3.0 2.0 0.250   0.000 0   0.000
     3.0 2.0 0.125   0.000 0   0.000];
H = info.history(1:10, [2 3 4 6 7 8]);
c = [confere(info.ftrace(1), 170, 0), confere(fv, [170 146 130 106], 0), ...
     confere(H, T, 3), info.nit == 20, info.nfev == 81, ...
     confere(x, [3 2], 3), confere(fx, 0, 3)];
fprintf(['  f(0,0) = 170: %s | vértices 170/146/130/106: %s | tabela k = 0..9: %s | ' ...
         '20 it.: %s | n_f = 81: %s | x = (3,2): %s | f = 0: %s\n'], simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ Rosenbrock
x0 = [-1.5 2];  Delta = [1 1];  tolx = 1e-6;
fprintf('\nRosenbrock a partir de (%g,%g), Delta_0 = (%g,%g), tolx = %g\n', x0, Delta, tolx);
[x, fx, info] = BoxEvo(rosen, x0, Delta, tolx);
fprintf('%s\n', info.message);
hit = find(info.ftrace < 1e-4, 1);
kh = ceil((hit - 1)/4);                % linha de info.history (iteração k = kh - 1)
Dh = info.history(kh, 4);              % Delta_1 nessa iteração
fprintf('primeira avaliação com f < 1e-4: n_f = %d (na iteração k = %d, com Delta = %.2e)\n', ...
        hit, kh - 1, Dh);
fprintf('no fim: x = (%.6f, %.6f), f = %.2e;  %d iterações, n_f total = %d\n', ...
        x, fx, info.nit, info.nfev);
c = [hit == 13445, confere(Dh, 6e-5, 5)];
fprintf('  n_f = 13 445 até f < 1e-4: %s | Delta ~ 6e-5: %s\n', simnao{c + 1});
ok = ok && all(c);

% --------------------------------------- «Porque falha em vales curvos?»
[~, ~, info] = BoxEvo(rosen, [0 0], [0.5 0.5], 0, 1);    % uma só iteração
fd = info.ftrace(2:5);  fe1 = rosen([0.25 0]);
nomes = {'encolhe', 'move'};
fprintf('\nRosenbrock em (0,0), f = %g, passo 0.25: diagonais f = %s -> %s\n', ...
        info.ftrace(1), sprintf('%.2f ', fd), nomes{info.history(1, 7) + 1});
fprintf('eixo +e1: f(0.25, 0) = %.3f < 1\n', fe1);
c = [all(fd > 1), confere(min(fd), 4.08, 2), confere(max(fd), 11.33, 2), ...
     info.history(1, 7) == 0, confere(fe1, 0.953, 3)];
fprintf('  nenhuma diagonal melhora: %s | 4.08: %s | 11.33: %s | encolhe: %s | 0.953: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
