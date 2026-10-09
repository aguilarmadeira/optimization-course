# 2.4 — Bisseção (MATLAB)

Estado: **disponível**.

- Função: `Bisection.m` — `[x, fx, info] = Bisection(f, df, a, b, tolg, kmax, tolx, verbose)` (paragem |f′(z)| ≤ `tolg`; `kmax` = 100, `tolx` = 0 e `verbose` = false por omissão; `f` só para leitura, pode ser `[]`)
- Exemplo: `ex02_4_bisection.m`, que reproduz os números do deck 2.4.
- Apontamentos: `notes/pt/` (deck 2.4).

## O algoritmo (como nas aulas e no Deb)

1. **Passo 1.** Escolher a < b com f′(a) < 0 < f′(b) e ε > 0; fazer x₁ = a e x₂ = b.
2. **Passo 2.** Fazer z = (x₁ + x₂)/2 e calcular f′(z).
3. **Passo 3.** Se |f′(z)| ≤ ε, terminar: x* ≈ z. Senão, se f′(z) < 0, fazer x₁ = z; se f′(z) > 0, fazer x₂ = z. Voltar ao passo 2.

|f′(z)| ≤ ε não garante que z esteja a menos de ε de x*. Variante (`tolx` > 0): parar também quando x₂ − x₁ ≤ tolx e devolver o ponto médio, com |x − x*| ≤ tolx/2. Com `tolg` = 0 é a bisseção com N = ⌈log₂((b − a)/tolx)⌉ reduções.

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex02_4_bisection
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex02_4_bisection"`). A ajuda da função: `help Bisection`.

## O que o exemplo imprime

As tabelas das iterações (k = 1, 2, …) dos exemplos 2 (x² + 54/x em [1,4], ε = 0,05) e 1 (−2 sin x + x²/16 em [−1,3], ε = 0,03), os valores finais e as contagens, a variante com 15 reduções, os dois exercícios «Para resolver» (ε = 0,01), o exemplo «Quantas iterações?» e a verificação com os slides. A coluna f(z) é só para leitura. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exemplo 2: 7 iterações, z = 2.9922 (f′(z) = −0.0470), n_g = 7 + 2 = 9, n_f = 0. Variante com 15 reduções: x₂ − x₁ = 3/2¹⁵ = 9.2e-05; a secção áurea (L = 0.0022) fica 24 vezes maior.
- Exemplo 1: 7 iterações, z = 1.4688 (f′(z) = −0.0201; exato 1.4783). Com ε = 0,05: 3 iterações, z = 1.5.
- Para resolver: 6 iterações, 2.5781 / 8.1355 (exato 2.5678: erro 0.0104 > ε); 10 iterações, 1.4492 / −2.7153.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev` (avaliações de f), `ngev` (de f'), `nhev` (de f''), `nit`, `flag` (0: |f′(z)| ≤ tolg; 1: kmax; 2: variante), e a tabela das iterações em `info.history`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
