# Simulação de Monte Carlo no Ox

## Descrição:

Código produzido no Ox para realizar a simulação de Monte Carlo, tendo como objetivo estimar os parâmetros da distribuição Gaussiana Inversa. Tem-se como método de  maximização da distribuição o processo de Máxima Verossimilhança. Vale ressaltar que a simulação se utiliza de geradores de números aleatórios para recriar processos estocásticos e assim prever o valor ótimo dos parâmetros, o fazem repetidas vezes (R = 10000), para cada passagem da simulação ocorre a otimização e para este código é considerado método não-linear Quasi-Newton com algoritmo BFGS. Para a análise da otimização é calculado para cada parâmetro estimado: o Viés, Viés Relativo, e o Erro Quadrático Médio.

## O que é Simulação de Monte Carlo:

A simulação serve para estimar os possíveis resultados de um evento incerto. Pode ser utilizada para avaliar o impacto de risco em muitos cenários da vida real, como em  inteligência artificial, preços de ações, previsão de vendas, gerenciamento de projetos e precificação.

## Bibliotecas

 - oxstd.oxh
 - oxdraw.oxh
 - oxprob.oxh 
 - maximize
 - oxfloat.oxh
