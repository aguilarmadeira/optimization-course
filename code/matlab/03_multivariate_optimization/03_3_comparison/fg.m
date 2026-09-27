function [f, g, H] = fg(x)
%FG  Função de Rosenbrock, com gradiente e Hessiana só quando são pedidos.
%
%   f         = fg(x)      só f (1 avaliação de f)
%   [f, g]    = fg(x)      f e o gradiente (a forma que o fminunc espera com
%                          'SpecifyObjectiveGradient' / 'GradObj' = 'on')
%   [f, g, H] = fg(x)      f, gradiente e Hessiana
%
%   f(x) = (1 - x1)^2 + 100 (x2 - x1^2)^2, mínimo f = 0 em (1, 1).
%   O gradiente só é calculado quando nargout > 1 e a Hessiana quando
%   nargout > 2: quem só precisa de f (p. ex. a pesquisa em linha) não paga
%   o resto.
%
%   Uso no deck 3.3 («Os contadores já feitos»):
%     opts = optimoptions('fminunc', 'SpecifyObjectiveGradient', true);
%     [x, fx, ~, out] = fminunc(@fg, [-1.5; 2], opts);
%     out.funcCount, out.iterations
%   (em Octave: opts = optimset('GradObj', 'on')).
%
%   Otimização — deck 3.3 (Métodos com derivadas) e comparação do cap. 3.
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

f = (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
if nargout > 1                     % gradiente só se for pedido
  g = [-2*(1 - x(1)) - 400*x(1)*(x(2) - x(1)^2)
        200*(x(2) - x(1)^2)];
end
if nargout > 2                     % Hessiana só se for pedida
  H = [2 - 400*(x(2) - 3*x(1)^2), -400*x(1)
       -400*x(1),                  200];
end
end
