% EX03_2_2_GRID_SEARCH  Reproduz os exemplos do deck 3.2.2 (pesquisa em grelha).
%
%   Himmelblau em [-5,5]^2:
%   1. Grelha 10 x 10 (100 avaliações): melhor f = 4.51; os 4 melhores nós,
%      f = 4.51, 4.70, 6.15, 11.88, ficam um em cada bacia.
%   2. Grosseira -> fina: mais 100 avaliações (grelha 10 x 10 em
%      [x - h, x + h], h = 10/9) à volta do melhor nó: f = 0.32.
%   3. Grelha + Nelder–Mead (3.2.3) a partir de cada um dos 4 melhores nós
%      (simplex do deck, tolerâncias 1e-6): n_f = 89, 94, 94, 87; os 4
%      mínimos globais, todos com f < 1e-11, em 464 avaliações no total;
%      o 5.º melhor nó leva outra vez a (3,2), com mais 95 avaliações.
%
%   Usa NelderMead da pasta 03_2_3_nelder_mead (posta no caminho por
%   uc_setup; não há cópia da função nesta pasta).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 3.2.2.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 3.2.2: pesquisa em grelha\n');

himmel = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
MIN = [3 2; -2.805118 3.131312; -3.779310 -3.283186; 3.584428 -1.848126];   % os 4 mínimos
bacia = @(p) find(sqrt(sum((MIN - repmat(p, 4, 1)).^2, 2)) == ...
                  min(sqrt(sum((MIN - repmat(p, 4, 1)).^2, 2))), 1);   % mínimo mais próximo

% ------------------------------------------------------------ 1. grelha 10 x 10
fprintf('\n1. Himmelblau em [-5,5]^2, grelha 10 x 10\n');
[x, fx, info] = GridSearch(himmel, [-5 -5], [5 5], 10, true);
top4 = info.history(1:4, :);
b4 = zeros(1, 4);  for i = 1:4, b4(i) = bacia(top4(i, 1:2)); end
fprintf('melhor nó (%.4f, %.4f), f = %.2f;  n_f = %d;  bacias dos 4 melhores: %d %d %d %d\n', ...
        x, fx, info.nfev, b4);
c = [confere(fx, 4.51, 2), info.nfev == 100, confere(top4(:, 3), [4.51 4.70 6.15 11.88], 2), ...
     numel(unique(b4)) == 4];
fprintf('  f = 4.51: %s | n_f = 100: %s | 4 melhores: %s | um em cada bacia: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 2. grosseira -> fina
h = info.h;
[x2, f2, info2] = GridSearch(himmel, x - h, x + h, 10);
fprintf('\n2. Grosseira -> fina: grelha 10 x 10 em [x - h, x + h], h = 10/9 = %.4f\n', h(1));
fprintf('melhor nó (%.4f, %.4f), f = %.2f;  mais %d avaliações\n', x2, f2, info2.nfev);
c = [confere(f2, 0.32, 2), info2.nfev == 100];
fprintf('  f = 0.32: %s | +100 avaliações: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 3. grelha + Nelder–Mead
fprintf('\n3. Grelha + Nelder–Mead a partir dos 5 melhores nós (tolx = tolf = 1e-6)\n');
fprintf('  nó      ponto inicial         f                 mínimo    f final   n_f\n');
tot = info.nfev;  nf = zeros(1, 5);  ff = zeros(1, 5);  bm = zeros(1, 5);
for i = 1:5
  p = info.history(i, 1:2);
  [xm, fm, im] = NelderMead(himmel, p, 1e-6, 1e-6);   % simplex do deck a partir de p
  nf(i) = im.nfev;  ff(i) = fm;  bm(i) = bacia(xm);
  fprintf('%3d  (%7.4f, %7.4f) %9.2f  (%9.6f, %9.6f) %10.1e %5d\n', i, p, info.history(i, 3), xm, fm, im.nfev);
  if i <= 4, tot = tot + im.nfev; end
end
fprintf('total grelha + 4 Nelder–Mead: %d avaliações;  o 5.º nó leva ao mínimo %d, com mais %d\n', ...
        tot, bm(5), nf(5));
c = [isequal(nf(1:4), [89 94 94 87]), all(ff(1:4) < 1e-11), numel(unique(bm(1:4))) == 4, ...
     tot == 464, bm(5) == 1, nf(5) == 95];
fprintf(['  n_f = 89, 94, 94, 87: %s | f < 1e-11: %s | 4 mínimos diferentes: %s | total 464: %s\n' ...
         '  5.º nó -> (3,2): %s | mais 95: %s\n'], simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
