% EX06_6_GAMULTIOBJ  O exemplo convexo da UC com gamultiobj, como no slide.
%
%   min f1(x) = x1^2 + x2^2,  min f2(x) = (x1 - 1)^2 + x2^2
%   s.a. x1^2 + x2^2 <= 1,  0 <= x1, x2 <= 1
%
%   Frame «No computador: MATLAB e Python» do deck 6.6:
%       [x,fval] = gamultiobj(fun, nvars, [],[],[],[], lb, ub, nonlcon);
%   gamultiobj (Global Optimization Toolbox) é uma VARIANTE da família
%   NSGA-II (NSGA-II controlado, com distância de aglomeração e fração de
%   Pareto); não dá os mesmos números que NSGA2.m nem que o pymoo.
%   Aqui: população 40, 50 gerações, rng(1); o HV usa Hypervolume2D com
%   r = (1.1; 1.65), o mesmo dos slides (frente exata: 1.648).
%   Por omissão, gamultiobj devolve só ParetoFraction = 35 % da população
%   (cerca de 14 pontos); com 'ParetoFraction', 1 devolve até 40.
%
%   NÃO TESTADO NO OCTAVE: o GNU Octave não tem gamultiobj. Escrito para
%   MATLAB com a Global Optimization Toolbox; no Octave o script só avisa.
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 6.6.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

if ~exist('gamultiobj', 'file')
  fprintf('gamultiobj não disponível (Global Optimization Toolbox; não existe no Octave).\n');
  fprintf('Use NSGA2.m (ver ex06_6_nsga2.m), que não depende de toolboxes.\n');
  return
end

fun = @(x) [x(1)^2 + x(2)^2, (x(1) - 1)^2 + x(2)^2];
nvars = 2;
lb = [0 0];  ub = [1 1];
nonlcon = @(x) deal(x(1)^2 + x(2)^2 - 1, []);   % c(x) <= 0, sem igualdades

rng(1);
opts = optimoptions('gamultiobj', 'PopulationSize', 40, 'MaxGenerations', 50);
% opts = optimoptions(opts, 'ParetoFraction', 1);   % devolver até 40 pontos
[x, fval, exitflag, output] = gamultiobj(fun, nvars, [], [], [], [], lb, ub, nonlcon, opts);

hv = Hypervolume2D(fval, [1.1 1.65]);
fprintf('gamultiobj: %d pontos, %d avaliações (output.funccount), %d gerações, HV = %.3f (frente exata: 1.648)\n', ...
        size(fval, 1), output.funccount, output.generations, hv);
fprintf('exitflag = %d: %s\n', exitflag, output.message);

t = linspace(0, 1, 200);
plot(t.^2, (1 - t).^2, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);  hold on
plot(fval(:, 1), fval(:, 2), 'o', 'Color', [190 30 45]/255);  hold off
axis equal;  xlabel('f_1');  ylabel('f_2');
legend('frente exata', 'gamultiobj');
title(sprintf('gamultiobj, população 40, 50 gerações: HV = %.3f', hv));
