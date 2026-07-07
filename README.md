# Ordenamento de SNPs por Metodos Estatisticos e de Aprendizado de Maquina

Artigo submetido ao **XXIX ENMC -- Encontro Nacional de Modelagem Computacional** e **XVII ECTM -- Encontro de Ciencia e Tecnologia de Materiais** (06 a 09 de outubro de 2026, Bento Goncalves -- RS).

## Resumo

Comparacao de seis metodos de ordenamento de SNPs (*Single Nucleotide Polymorphisms*) em duas bases de dados simuladas com cenarios de efeitos geneticos independentes e interacoes epistaticas. Os metodos avaliados sao: Correlacao Marginal, Correlacao Parcial, CAR Score, I Score, Importancia por Permutacao e Importancia por Impureza.

## Estrutura do Repositorio

```
.
├── artigo/                   # Artigo em LaTeX
│   ├── main.tex              # Documento principal
│   ├── bib/                  # Referencias bibliograficas
│   │   └── Referencias.bib
│   ├── figures/              # Figuras do artigo
│   └── setup/                # Estilo do evento
│       └── config.sty
├── R/                        # Codigo R das analises
│   ├── simulacao_base_dados_1.R   # Simulacao da Base de Dados 1
│   ├── simulacao_base_dados_2.R   # Simulacao da Base de Dados 2
│   ├── scoring_methods.R          # Metodos de ordenamento
│   ├── gerar_tabelas_artigo.R     # Geracao das tabelas e figuras
│   ├── utils.R                    # Funcoes auxiliares e instalacao de pacotes
│   └── resultados/                # CSVs e tabelas geradas
├── referencias/              # PDFs de referencia
└── .gitignore
```

## Dependencias

### R

- **R >= 4.1** (testado com R 4.4.x)

### Pacotes R

| Pacote | Finalidade |
|--------|-----------|
| `scrime` | Simulacao de dados genotipicos (`simulateSNPglm`) |
| `randomForest` | Ajuste da Random Forest e calculo das importancias |
| `care` | Calculo do CAR Score e da Correlacao Marginal (`carscore`) |
| `corpcor` | Correlacao Parcial shrinkage (`pcor.shrink`, `estimate.lambda`) |

Os pacotes sao instalados automaticamente pelo script `utils.R` caso nao estejam presentes.

### LaTeX

- Distribuicao TeX com `pdflatex` e `bibtex` (ex.: MiKTeX, TeX Live)

## Como Reproduzir os Resultados

Execute os scripts R na ordem abaixo a partir da raiz do repositorio. Todos usam `set.seed(123)` para garantir reprodutibilidade.

```r
# 1. Carregar funcoes auxiliares e instalar pacotes (se necessario)
source("R/utils.R")

# 2. Simular as bases de dados
source("R/simulacao_base_dados_1.R")   # gera a Base de Dados 1 (Simulacao 1)
source("R/simulacao_base_dados_2.R")   # gera a Base de Dados 2 (Simulacao 4)

# 3. Calcular os scores e rankings por todos os metodos
source("R/scoring_methods.R")

# 4. Gerar as tabelas LaTeX e figuras do artigo
source("R/gerar_tabelas_artigo.R")     # salva em R/resultados/ e artigo/figures/
```

Os resultados (CSVs e tabelas `.tex`) serao gravados em `R/resultados/`.

## Como Compilar o Artigo

```bash
cd artigo
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Autores

- Joao Vitor Barreto Baptista (UFF)
- Thiago Jordem Pereira (UFF)
- Fabrizzio Conde de Oliveira (UFF)

## Licenca

Este repositorio contem material academico. Consulte os autores antes de reutilizar.
