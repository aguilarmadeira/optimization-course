function uc_setup()
%UC_SETUP  Acrescenta ao caminho do MATLAB/Octave todas as pastas do código da UC.
%
%   uc_setup
%
%   Faz addpath(genpath(<pasta deste ficheiro>)), isto é, de code/matlab e
%   de todas as subpastas (common/ e as pastas dos métodos), sem as pastas
%   escondidas (nome começado por '.'; o genpath do Octave incluí-las-ia).
%   Pode chamar-se várias vezes (addpath não duplica entradas).
%
%   Os exemplos ex*.m chamam-no na primeira linha de código. Para usar as
%   funções da UC nos seus próprios scripts:
%       addpath('<repo>/code/matlab'); uc_setup
%
%   Otimização — código da UC.  J. F. A. Madeira — Licença MIT.

raiz = fileparts(mfilename('fullpath'));
anterior = cd(raiz);            % forma canónica do caminho (sem '..')
raiz = pwd;
cd(anterior);

pastas = strsplit(genpath(raiz), pathsep);
manter = true(size(pastas));
for j = 1:numel(pastas)
  resto = pastas{j}(numel(raiz) + 1:end);          % parte abaixo de raiz
  manter(j) = ~isempty(pastas{j}) && isempty(regexp(resto, '[\\/]\.', 'once'));
end
addpath(strjoin(pastas(manter), pathsep));
end
