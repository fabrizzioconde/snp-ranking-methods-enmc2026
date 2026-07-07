# Ordenamento de SNPs por Métodos Estatísticos e de Aprendizado de Máquina

Artigo submetido ao **XXIX ENMC -- Encontro Nacional de Modelagem Computacional** e **XVII ECTM -- Encontro de Ciência e Tecnologia de Materiais** (06 a 09 de outubro de 2026, Bento Gonçalves -- RS).

## Resumo

Comparação de seis métodos de ordenamento de SNPs (*Single Nucleotide Polymorphisms*) em duas bases de dados simuladas com cenários de efeitos genéticos independentes e interações epistáticas. Os métodos avaliados são: Correlação Marginal, Correlação Parcial, CAR Score, I Score, Importância por Permutação e Importância por Impureza.

## Estrutura do Repositório

```
.
├── artigo/                   # Artigo em LaTeX
│   ├── main.tex              # Documento principal
│   ├── bib/                  # Referências bibliográficas
│   │   └── Referencias.bib
│   ├── figures/              # Figuras do artigo
│   └── setup/                # Estilo do evento
│       └── config.sty
├── R/                        # Código R das análises
│   ├── simulacao_base_dados_1.R   # Simulação da Base de Dados 1
│   ├── simulacao_base_dados_2.R   # Simulação da Base de Dados 2
│   ├── scoring_methods.R          # Métodos de ordenamento
│   ├── gerar_tabelas_artigo.R     # Geração das tabelas e figuras
│   ├── utils.R                    # Funções auxiliares e instalação de pacotes
│   └── resultados/                # CSVs e tabelas geradas
├── referencias/              # PDFs de referência
└── .gitignore
```

## Dependências

### R

- **R >= 4.1** (testado com R 4.4.x)

### Pacotes R

| Pacote         | Finalidade                                                        |
| -------------- | ----------------------------------------------------------------- |
| `scrime`       | Simulação de dados genotípicos (`simulateSNPglm`)                 |
| `randomForest` | Ajuste da Random Forest e cálculo das importâncias                |
| `care`         | Cálculo do CAR Score e da Correlação Marginal (`carscore`)        |
| `corpcor`      | Correlação Parcial shrinkage (`pcor.shrink`, `estimate.lambda`)   |

Os pacotes são instalados automaticamente pelo script `utils.R` caso não estejam presentes.

### LaTeX

- Distribuição TeX com `pdflatex` e `bibtex` (ex.: MiKTeX, TeX Live)

## Como Reproduzir os Resultados

Execute os scripts R na ordem abaixo a partir da raiz do repositório. Todos usam `set.seed(123)` para garantir reprodutibilidade.

```r
# 1. Carregar funções auxiliares e instalar pacotes (se necessário)
source("R/utils.R")

# 2. Simular as bases de dados
source("R/simulacao_base_dados_1.R")   # gera a Base de Dados 1 (Simulação 1)
source("R/simulacao_base_dados_2.R")   # gera a Base de Dados 2 (Simulação 4)

# 3. Calcular os scores e rankings por todos os métodos
source("R/scoring_methods.R")

# 4. Gerar as tabelas LaTeX e figuras do artigo
source("R/gerar_tabelas_artigo.R")     # salva em R/resultados/ e artigo/figures/
```

Os resultados (CSVs e tabelas `.tex`) serão gravados em `R/resultados/`.

## Como Compilar o Artigo

```bash
cd artigo
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Autores

- João Vitor Barreto Baptista (UFF)
- Thiago Jordem Pereira (UFF)
- Fabrízzio Condé de Oliveira (UFF)

## Licença

Este repositório contém material acadêmico. Consulte os autores antes de reutilizar.
