% EX06_4_EPSILON_CONSTRAINT  Reproduz os exemplos do deck 6.4 (epsilon-constrangimento).
%
%   Exemplo não convexo da UC (6.2):  min (f1, f2) = (x1, 1 - x1^2 + x2),
%   x em [0,1] x [0, 0.6].  Subproblema:  min f2  s.a.  f1 <= e.
%     1. «O mesmo exemplo do 6.3, à mão»: e = 0.25, 0.50, 0.75.
%     2. «Em MATLAB: varrer e com fmincon»: e = 0.05:0.1:0.95, x0 = (0.5, 0.3),
%        arranque a quente; 10 pontos da frente f2 = 1 - f1^2, u = 2e, n_f.
%     3. «O multiplicador: a taxa de troca»: e = 0.5 -> u = 1; f2(0.51) = 0.7399;
%        e = 0.9 -> u = 1.8, a mesma folga compra ~0.018.
%     4. «Restrição inativa»: e = 1.05 e 1.2 -> a mesma solução (1, 0), u = 0.
%     5. «Nantes–Lille»: o mais barato com duração <= e (e = 4, 5, 10 h).
%
%   Contagens: as do slide (60 avaliações de f2, 90 de f1, 150 chamadas) são
%   do SLSQP do SciPy e só se verificam na versão Python. Aqui o solver é o
%   fmincon (MATLAB) ou o sqp (GNU Octave): os pontos e os multiplicadores
%   coincidem; as contagens são as desse solver e imprimem-se só a título
%   informativo (não entram na verificação).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 6.4.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
perto = @(v, s, tol) all(abs(v(:) - s(:)) <= tol);
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 6.4: método do epsilon-constrangimento\n');

f1 = @(x) x(1);  f2 = @(x) 1 - x(1)^2 + x(2);
lb = [0 0];  ub = [1 0.6];  x0 = [0.5 0.3];

% ------------------------------------------------ 1. À mão: três valores
fprintf('\n1. min f2 = 1 - x1^2 + x2  s.a.  f1 = x1 <= e  (x* = (e, 0))\n');
[~, F, ~] = EpsilonConstraint(f1, f2, x0, lb, ub, [0.25 0.50 0.75]);
for k = 1:3
  fprintf('   e = %.2f  =>  (f1, f2) = (%.4f; %.4f)\n', F(k, 1), F(k, 1), F(k, 2));
end
c = confere(F, [0.25 0.9375; 0.50 0.75; 0.75 0.4375], 4);
fprintf('  (0.25; 0.9375), (0.50; 0.75), (0.75; 0.4375): %s\n', simnao{c + 1});
ok = ok && c;

% ------------------------------------------------ 2. Varrimento do slide
fprintf('\n2. Varrimento e = 0.05:0.1:0.95 a partir de x0 = (%g, %g), arranque a quente\n', x0);
e = 0.05:0.1:0.95;
[x, F, info] = EpsilonConstraint(f1, f2, x0, lb, ub, e, true);
fprintf('   solver: %s (%s)\n', info.solver, info.message);
nd = size(unique(round(F*1e6), 'rows'), 1);
fprintf('   %d valores de e -> %d pontos distintos, todos na frente f2 = 1 - f1^2\n', numel(e), nd);
fprintf('   contagens (%s): n_f = %d; chamadas a f2 = %d, a f1 = %d (total %d); %d iterações\n', ...
        info.solver, info.nfev, info.nf2, info.nf1, info.ncalls, info.nit);
fprintf('   slide (SLSQP do SciPy): 60 avaliações de f2, 90 de f1 (150 chamadas) — só se\n');
fprintf('   verificam em Python; com o %s as contagens são outras (informativo).\n', info.solver);
c = [numel(e) == 10, nd == 10, perto(x, [e' zeros(10, 1)], 1e-5), ...
     perto(F(:, 2), 1 - e'.^2, 1e-5), perto(info.u, 2*e', 1e-4), info.flag == 0];
fprintf('  10 valores: %s | 10 pontos: %s | x* = (e, 0): %s | f2 = 1 - e^2: %s | u = 2e: %s | sem falhas: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 3. O multiplicador
fprintf('\n3. Multiplicador de f1 <= e (g = e - f1 >= 0): u = 2e = -d f2*/de\n');
[x, F, info] = EpsilonConstraint(f1, f2, x0, lb, ub, [0.50 0.51 0.90 0.91]);
d1 = F(2, 2) - F(1, 2);  d2 = F(3, 2) - F(4, 2);
fprintf('   e = 0.50: x* = (%.4f; %.4f), f2* = %.4f, u = %.4f\n', x(1, :), F(1, 2), info.u(1));
fprintf('   e = 0.51: f2* = %.4f; diferença %.4f  (aproximação -u*0.01 = %.4f)\n', ...
        F(2, 2), d1, -info.u(1)*0.01);
fprintf('   e = 0.90: u = %.4f; a folga 0.90 -> 0.91 compra %.4f em f2\n', info.u(3), d2);
c = [confere(x(1, 1), 0.5, 1), confere(F(1, 2), 0.75, 2), confere(info.u(1), 1, 3), ...
     confere(F(2, 2), 0.7399, 4), confere(d1, -0.0101, 4), confere(info.u(3), 1.8, 3), ...
     confere(d2, 0.018, 3)];
fprintf('  x* = 0.5: %s | f2* = 0.75: %s | u = 1: %s | f2(0.51) = 0.7399: %s | -0.0101: %s | u(0.9) = 1.8: %s | ~0.018: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 4. Restrição inativa
fprintf('\n4. Restrição inativa: e > 1\n');
[x, F, info] = EpsilonConstraint(f1, f2, x0, lb, ub, [1.05 1.2]);
for k = 1:2
  fprintf('   e = %.2f: x* = (%.4f; %.4f), f1 = %.4f < e, u = %.4f\n', ...
          info.history(k, 1), x(k, :), F(k, 1), info.u(k));
end
c = perto(x, [1 0; 1 0], 1e-6) && perto(info.u, [0; 0], 1e-6);
fprintf('  mesma solução (1, 0) com u = 0: %s\n', simnao{c + 1});
ok = ok && c;

% ------------------------------------------------ 5. Nantes–Lille (6.1)
fprintf('\n5. Nantes–Lille: o mais barato com duração <= e\n');
nomes = {'Carro via Paris', 'Carro via Rouen', 'Carro sem portagens', 'Autocarro 1', ...
         'Autocarro 2', 'TGV', 'TGV low-cost', 'Avião'};
dur = [5+33/60, 5+49/60, 7+18/60, 9.5, 9.5, 4+2/60, 4+20/60, 1+5/60];
custo = [230 251 204 38 28 84 40 134];
esc = zeros(1, 3);  E = [4 5 10];
for k = 1:3
  adm = find(dur <= E(k));
  [~, j] = min(custo(adm));                   % empate: fica o primeiro
  esc(k) = adm(j);
  fprintf('   e = %2d h: %s (%d €)\n', E(k), nomes{esc(k)}, custo(esc(k)));
end
c = isequal(esc, [8 7 5]) && isequal(custo(esc), [134 40 28]);
fprintf('  avião (134), TGV low-cost (40), autocarro 2 (28): %s\n', simnao{c + 1});
ok = ok && c;

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
