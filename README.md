# rsa-cryptanalysis-demo

Repositório de código para o **Relatório Técnico RSA** da disciplina
*Introdução à Cibersegurança* — Pós-graduação *Lato Sensu* em
Comunicação Quântica — **SENAI CIMATEC**.

> **Objetivo educacional:** demonstrar o funcionamento interno do RSA,
> seus ataques clássicos conhecidos e a ameaça imposta pelo algoritmo
> de Shor executado em simulador quântico (Qiskit Aer).

---

## Estrutura do repositório

```
rsa-cryptanalysis-demo/
├── rsa_core.py          # Módulo 1 — RSA puro: chaves, cifragem, assinatura
├── rsa_attacks.py       # Módulo 2 — Ataques clássicos ao RSA
├── shor_simulator.py    # Módulo 3 — Algoritmo de Shor com Qiskit
└── README.md
```

---

## Pré-requisitos

| Requisito | Versão mínima |
|---|---|
| Python | 3.10 |
| qiskit | 2.0 |
| qiskit-aer | 0.17 |

> Os módulos `rsa_core.py` e `rsa_attacks.py` **não requerem pacotes externos** —
> usam apenas a biblioteca padrão do Python.
> Somente `shor_simulator.py` depende do Qiskit.

### Instalação das dependências

```bash
# Clone o repositório
git clone https://github.com/jullyanolino/rsa-cryptanalysis-demo.git
cd rsa-cryptanalysis-demo

# (Opcional) Crie um ambiente virtual
python -m venv .venv
source .venv/bin/activate          # Linux/macOS
# .venv\Scripts\activate           # Windows

# Instale as dependências
pip install qiskit qiskit-aer

# Para execução em hardware IBM Quantum real (opcional)
pip install qiskit-ibm-runtime
```

---

## Módulo 1 — `rsa_core.py`

Implementação educacional do RSA **sem dependências externas**.

### Funcionalidades

- Teste de primalidade probabilístico (Miller-Rabin, 40 rodadas)
- Geração de primos criptograficamente adequados (`secrets.randbits`)
- Geração de par de chaves RSA (`gerar_chaves`)
- Cifragem e decifragem *textbook* (`cifrar_textbook`, `decifrar_textbook`)
- Cifragem e decifragem com OAEP simplificado (`cifrar`, `decifrar`)
- Assinatura digital e verificação SHA-256 + RSA (`assinar`, `verificar_assinatura`)
- Modo demonstrativo: reproduz o **exemplo numérico do relatório** (p=61, q=53)

### Execução

```bash
# Demonstração completa (exemplo numérico + OAEP com chave de 1024 bits)
python rsa_core.py

# Especificar tamanho da chave (em bits)
python rsa_core.py 2048
```

### Saída esperada (trecho)

```
EXEMPLO NUMÉRICO: p=61, q=53
  p            = 61
  q            = 53
  n  = p·q     = 3233
  φ(n)         = 3120
  e            = 17   [mdc(17,3120) = 1  ✓]
  d  = e⁻¹ mod φ(n) = 2753

  C = M^e mod n = 65^17 mod 3233 = 2790
  M = C^d mod n = 2790^2753 mod 3233 = 65  ✓
```

### API pública

```python
from rsa_core import gerar_chaves, cifrar, decifrar, assinar, verificar_assinatura

# Gerar par de chaves (2048 bits recomendado em produção)
par = gerar_chaves(bits=2048)

# Cifrar com RSA-OAEP
cifrado   = cifrar(par.publica, b"Mensagem secreta")
recuperado = decifrar(par.privada, cifrado)

# Assinatura digital
assinatura = assinar(par.privada, b"Documento")
valida     = verificar_assinatura(par.publica, b"Documento", assinatura)
```

---

## Módulo 2 — `rsa_attacks.py`

Simulação dos principais **ataques clássicos** ao RSA.

### Ataques implementados

| # | Ataque | Pré-condição | Complexidade |
|---|--------|-------------|-------------|
| 1 | **Fatoração por Força Bruta** | n pequeno | O(√n) |
| 2 | **Wiener** | d < n^{1/4}/3 | O(log² n) |
| 3 | **Módulo Comum** | n compartilhado, e₁ ≠ e₂ | O(log n) |
| 4 | **Hastad** | e=3 sem padding, 3 destinatários | O(log n) + raiz cúbica |

### Execução

```bash
python rsa_attacks.py
```

### Saída esperada (trecho — Ataque de Wiener)

```
ATAQUE 2: Wiener – Expoente Privado Pequeno
  Gerando parâmetros RSA com d pequeno (256 bits)... concluído.

  [Segredo do experimento] d real = 1362452641891029727
  Executando ataque de Wiener...
    *** SUCESSO ***
    d encontrado : 1362452641891029727
    Tempo total  : 7.394 ms  ✓
```

### API pública

```python
from rsa_attacks import (
    fatorar_forca_bruta,
    wiener_attack,
    ataque_modulo_comum,
    hastad_attack,
    crt_geral,
)

# Fatorar n pequeno
p, q = fatorar_forca_bruta(3233)

# Ataque de Wiener
d = wiener_attack(e, n)

# Ataque de módulo comum
M = ataque_modulo_comum(C1, C2, e1, e2, n)
```

---

## Módulo 3 — `shor_simulator.py`

Simulação do **Algoritmo de Shor** com Qiskit Aer (local) ou IBM Quantum (real).

### Casos suportados

| N | a | Período r | Fatores | Circuito |
|---|---|-----------|---------|---------|
| 15 | 7 | 4 | 3 × 5 | 8 qubits, profundidade ~48 |
| 21 | 2 | 6 | 3 × 7 | 11 qubits, profundidade ~256 |

### Execução

```bash
# Ambos os casos (N=15 e N=21) no simulador local
python shor_simulator.py

# Caso específico
python shor_simulator.py --n 15 --a 7

# Mais medições para maior precisão estatística
python shor_simulator.py --shots 4096

# Sem diagrama do circuito (execução mais limpa)
python shor_simulator.py --no-circuit

# Hardware IBM Quantum real (requer conta em quantum.ibm.com)
python shor_simulator.py --ibmq SEU_TOKEN_AQUI --n 15 --a 7
```

### Saída esperada (N=15, a=7)

```
[FASE 1] Pré-processamento Clássico
  N = 15, a = 7,  mdc(7,15)=1  ✓

[FASE 2] Execução Quântica
  Circuito: 8 qubits | AerSimulator | 2048 shots
  Portas: cswap=38, h=8, cp=6 | Profundidade: 48

[FASE 3] Pós-processamento (frações contínuas)
  r = 2  (candidato) → testando múltiplos → r = 4
  7^4 mod 15 = 1  ✓

[FASE 4] Cálculo dos Fatores
  mdc(7^2-1, 15) = mdc(3, 15) = 3
  mdc(7^2+1, 15) = mdc(5, 15) = 5

*** SUCESSO ***  15 = 3 × 5  ✓
```

### API pública

```python
from shor_simulator import fatorar_shor, construir_circuito_shor

# Fatorar N=15 no simulador local
resultado = fatorar_shor(N=15, a=7, backend="simulator", shots=2048)
if resultado:
    p, q = resultado
    print(f"15 = {p} × {q}")

# Construir apenas o circuito (para inspeção ou exportação)
qc = construir_circuito_shor(N=15, a=7)
print(qc.draw())

# Hardware IBM Quantum real
resultado = fatorar_shor(N=15, a=7,
                         backend="ibmq",
                         ibmq_token="SEU_TOKEN")
```

---

## Arquitetura geral

```
rsa_core.py
│
├── Primitivas matemáticas (mdc, Euclides estendido, inverso modular,
│   exponenciação modular, Miller-Rabin, geração de primos)
│
├── Geração de chaves (gerar_chaves)
├── Cifragem textbook e OAEP (cifrar_textbook, cifrar)
├── Decifragem (decifrar_textbook, decifrar)
└── Assinatura digital (assinar, verificar_assinatura)

rsa_attacks.py
│
├── Importa primitivas de rsa_core
├── Fatoração por força bruta (fatorar_forca_bruta)
├── Ataque de Wiener / frações contínuas (wiener_attack)
├── Ataque de módulo comum (ataque_modulo_comum)
└── Ataque de Hastad / CRT (hastad_attack, crt_geral)

shor_simulator.py
│
├── Construção do circuito quântico (construir_circuito_shor)
├── QFT inversa (qft_inversa)
├── Execução local: Qiskit Aer (executar_simulador)
├── Execução remota: IBM Quantum (executar_ibmq)
├── Extração do período (frações contínuas) (_extrair_periodo)
└── Fatoração completa (fatorar_shor)
```

---

## Referências

- RIVEST, R. L.; SHAMIR, A.; ADLEMAN, L. A Method for Obtaining Digital Signatures and
  Public-Key Cryptosystems. *Communications of the ACM*, v. 21, n. 2, 1978.
- SHOR, P. W. Algorithms for Quantum Computation: Discrete Logarithms and Factoring.
  *Proceedings of FOCS 1994*, pp. 124-134. IEEE, 1994.
- STAMP, M. *Information Security: Principles and Practice*. 2. ed. Wiley, 2011.
- PAAR, C.; PELZL, J. *Understanding Cryptography*. Springer, 2010.
- MENEZES, A. J.; VAN OORSCHOT, P. C.; VANSTONE, S. A.
  *Handbook of Applied Cryptography*. CRC Press, 2018.
- WEBBER, M. et al. The Impact of Hardware Specifications on Reaching Quantum Advantage
  in the Fault-Tolerant Regime. *AVS Quantum Science*, v. 4, n. 1, 2022.
- NIST FIPS 203: Module-Lattice-Based Key-Encapsulation Mechanism Standard, 2024.

---

## Licença

Código disponibilizado para fins educacionais.
Não deve ser utilizado em sistemas de produção.
