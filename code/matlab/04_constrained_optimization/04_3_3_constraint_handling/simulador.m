function y = simulador(x)
%SIMULADOR  «Simulador» do exemplo 4.3.3: uma caixa negra que falha.
%
%   y = simulador(x) devolve (x1 - 2)^2 + (x2 - 1)^2, ou NaN quando «a
%   simulação não converge» (aqui, quando x1 + x2 > 2), sem dizer porquê:
%   quem o usa não conhece g. Serve para ilustrar uma restrição oculta.
%
%   Otimização — deck 4.3.3.  J. F. A. Madeira — Licença MIT.

if x(1) + x(2) > 2
  y = NaN;
else
  y = (x(1) - 2)^2 + (x(2) - 1)^2;
end
end
