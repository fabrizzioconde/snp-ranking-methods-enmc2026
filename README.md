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
│   ├── gerar_tabelas_artigo.R     # Geracao das tabelas LaTeX
│   ├── utils.R                    # Funcoes auxiliares
│   └── resultados/                # CSVs e tabelas geradas
├── referencias/              # PDFs de referencia
└── .gitignore
```

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
