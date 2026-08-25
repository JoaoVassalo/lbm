[🇺🇸 English](README.md)

# LBM-CUDA

Solver bidimensional do método Lattice Boltzmann (LBM) desenvolvido em CUDA, com armazenamento baseado em momentos e execução paralela na GPU. Este é meu repositório mais recente e o principal foco do meu desenvolvimento atual, além de outros projetos nos quais colaboro.

## Estrutura

- `src/config/` — parâmetros de geometria, stencil D2Q9, momentos, física e tempo.
- `src/lbm/init/` — inicialização do domínio e da malha.
- `src/lbm/simulation/` — streaming, reconstrução das distribuições e colisão.
- `src/io/` — escrita dos resultados em VTK e cálculo de desempenho.
- `src/core/` — indexação e tipos auxiliares.
- `run.sh` — compila e executa o solver.

## Método

O código utiliza o stencil **D2Q9** e armazena seis variáveis por nó:

rho, ux, uy, mxx, mxy, myy

As distribuições são reconstruídas a partir desses momentos, seguidas pela etapa de streaming e colisão. A evolução utiliza dois buffers (`momA`/`momB`) em esquema ping-pong.

## Execução

É necessário ter o `nvcc` disponível no `PATH`.

```bash
bash run.sh
```

O script compila os arquivos `.cu` com C++20 e gera o executável em `build/debug`.

Os resultados são escritos em arquivos `.vti` em:

```text
results/plot/vtk/
```

Eles podem ser visualizados no ParaView.

## Configuração

Os principais parâmetros atualmente estão em `src/config/`, incluindo:

- tamanho do domínio;
- número de Reynolds;
- velocidade característica;
- viscosidade, `tau` e `omega`;
- parâmetros do D2Q9;
- intervalo e número de iterações.

Este repositório corresponde ao estado atual de desenvolvimento do solver, portanto a estrutura e os modelos podem evoluir conforme novas formulações e condições de contorno forem adicionadas.
