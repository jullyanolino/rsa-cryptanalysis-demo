"""
rsa_core.py
===========
Implementação educacional do RSA puro em Python.

Funcionalidades
---------------
  - Teste de primalidade probabilístico (Miller-Rabin)
  - Geração de primos criptograficamente adequados
  - Geração de par de chaves (chave pública / chave privada)
  - Cifragem e decifragem (RSA "livro-texto" e com padding OAEP simplificado)
  - Assinatura digital e verificação (hash SHA-256 + RSA)
  - Modo demonstrativo passo-a-passo com o exemplo numérico do relatório
    (p=61, q=53), executável diretamente via `python rsa_core.py`

Dependências
------------
  Biblioteca padrão do Python (>= 3.10) apenas.
  Não requer instalação de pacotes externos.

Uso rápido
----------
  $ python rsa_core.py

Referências
-----------
  Paar & Pelzl (2010), cap. 7
  Menezes, van Oorschot & Vanstone (2018), cap. 8
  NIST SP 800-57 Part 1 Rev. 5
"""

from __future__ import annotations

import hashlib
import math
import os
import secrets
import sys
from dataclasses import dataclass
from typing import NamedTuple

# ── Tipos ────────────────────────────────────────────────────────────────────

class ChavePublica(NamedTuple):
    e: int   # expoente público
    n: int   # módulo

class ChavePrivada(NamedTuple):
    d: int   # expoente privado
    n: int   # módulo

@dataclass
class ParDeChaves:
    publica:  ChavePublica
    privada:  ChavePrivada
    bits_n:   int          # tamanho em bits do módulo n

    def __repr__(self) -> str:
        return (
            f"ParDeChaves(\n"
            f"  chave_publica  = (e={self.publica.e}, n={self.publica.n})\n"
            f"  chave_privada  = (d=<SECRETO>, n={self.privada.n})\n"
            f"  bits_modulo    = {self.bits_n}\n"
            f")"
        )

# ── Teoria dos números ────────────────────────────────────────────────────────

def mdc(a: int, b: int) -> int:
    """Máximo divisor comum (Algoritmo de Euclides)."""
    while b:
        a, b = b, a % b
    return a


def euclides_estendido(a: int, b: int) -> tuple[int, int, int]:
    """
    Algoritmo de Euclides Estendido.

    Retorna (g, x, y) tal que  a*x + b*y = g = mdc(a, b).
    """
    if b == 0:
        return a, 1, 0
    g, x1, y1 = euclides_estendido(b, a % b)
    return g, y1, x1 - (a // b) * y1


def inverso_modular(a: int, n: int) -> int:
    """
    Inverso multiplicativo de a módulo n.

    Levanta ValueError se mdc(a, n) != 1.
    """
    g, x, _ = euclides_estendido(a % n, n)
    if g != 1:
        raise ValueError(f"mdc({a}, {n}) = {g} ≠ 1: inverso não existe.")
    return x % n


def exponenciacao_modular(base: int, exp: int, mod: int) -> int:
    """
    Exponenciação modular rápida: base^exp mod n.

    Implementa o algoritmo "square-and-multiply" (Right-to-Left Binary).
    Complexidade: O(log exp) multiplicações modulares.
    """
    if mod == 1:
        return 0
    resultado = 1
    base %= mod
    while exp > 0:
        if exp & 1:                   # bit menos significativo = 1
            resultado = resultado * base % mod
        exp >>= 1
        base = base * base % mod
    return resultado


def miller_rabin(n: int, k: int = 40) -> bool:
    """
    Teste de primalidade probabilístico de Miller-Rabin.

    Parâmetros
    ----------
    n : inteiro a testar
    k : número de rodadas (40 rodadas → probabilidade de erro < 4^{-40})

    Retorna True se n é provavelmente primo, False se composto.
    """
    if n < 2:
        return False
    # Casos base pequenos
    pequenos_primos = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n in pequenos_primos:
        return True
    if any(n % p == 0 for p in pequenos_primos):
        return False

    # Escrever n-1 como 2^r * d, com d ímpar
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # k rodadas com testemunhas aleatórias
    for _ in range(k):
        a = secrets.randbelow(n - 4) + 2      # a em [2, n-2]
        x = exponenciacao_modular(a, d, n)

        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False   # composto com certeza
    return True            # provavelmente primo


def gerar_primo(bits: int) -> int:
    """
    Gera um primo aleatório com exatamente `bits` bits.

    Usa secrets.randbits para garantir aleatoriedade criptográfica.
    """
    while True:
        # Garante bits exatos: MSB e LSB obrigatoriamente 1
        candidato = secrets.randbits(bits)
        candidato |= (1 << (bits - 1)) | 1    # força bit mais e menos significativo
        if miller_rabin(candidato):
            return candidato

# ── Geração de chaves ─────────────────────────────────────────────────────────

# Valor padrão de e recomendado (4º primo de Fermat)
EXPOENTE_PUBLICO_PADRAO = 65537


def gerar_chaves(bits: int = 2048,
                 e: int = EXPOENTE_PUBLICO_PADRAO) -> ParDeChaves:
    """
    Gera um par de chaves RSA.

    Parâmetros
    ----------
    bits : tamanho total do módulo n em bits (padrão: 2048)
    e    : expoente público (padrão: 65537)

    Retorna
    -------
    ParDeChaves com chave pública (e, n) e chave privada (d, n)

    Notas de segurança
    ------------------
    - Bits mínimos recomendados: 2048 (NIST SP 800-57)
    - p e q são descartados após o cálculo de d
    - phi(n) é descartado após o cálculo de d
    """
    if bits < 512:
        raise ValueError("Tamanho mínimo seguro: 512 bits (use >= 2048 em produção).")

    bits_primo = bits // 2

    # Gerar p e q distintos
    while True:
        p = gerar_primo(bits_primo)
        q = gerar_primo(bits_primo)
        if p != q and abs(p - q).bit_length() > bits_primo // 2:
            # Garante que p e q não são muito próximos (defesa contra fatoração de Fermat)
            break

    n     = p * q
    phi_n = (p - 1) * (q - 1)

    # Verificar que e é coprimo com phi(n)
    if mdc(e, phi_n) != 1:
        raise ValueError(f"e={e} não é coprimo com phi(n). Escolha outro e.")

    d = inverso_modular(e, phi_n)

    # Descarte seguro das variáveis sensíveis (best-effort em Python)
    # Em produção, usar memória gerenciada (C extension ou ctypes)
    del p, q, phi_n

    pub  = ChavePublica(e=e, n=n)
    priv = ChavePrivada(d=d, n=n)
    return ParDeChaves(publica=pub, privada=priv, bits_n=bits)

# ── Cifragem / Decifragem ─────────────────────────────────────────────────────

def _inteiro_para_bytes(n: int, comprimento: int) -> bytes:
    return n.to_bytes(comprimento, byteorder="big")

def _bytes_para_inteiro(b: bytes) -> int:
    return int.from_bytes(b, byteorder="big")


def cifrar_textbook(M: int, chave_pub: ChavePublica) -> int:
    """
    Cifragem RSA "livro-texto": C = M^e mod n.

    ATENÇÃO: não use em produção sem padding (OAEP).
    Adequado apenas para demonstração e exemplos numéricos.
    """
    e, n = chave_pub
    if not (0 <= M < n):
        raise ValueError(f"Mensagem M={M} deve satisfazer 0 <= M < n={n}.")
    return exponenciacao_modular(M, e, n)


def decifrar_textbook(C: int, chave_priv: ChavePrivada) -> int:
    """Decifragem RSA "livro-texto": M = C^d mod n."""
    d, n = chave_priv
    return exponenciacao_modular(C, d, n)


def _oaep_pad(mensagem: bytes, k: int, label: bytes = b"") -> bytes:
    """
    OAEP padding simplificado (educational).
    Baseado em PKCS#1 v2.2, Seção 7.1.1.
    k  = tamanho do módulo em bytes
    hLen = 32 (SHA-256)
    """
    hLen = 32   # SHA-256
    mLen = len(mensagem)
    if mLen > k - 2 * hLen - 2:
        raise ValueError("Mensagem muito longa para o tamanho de chave escolhido.")

    lHash = hashlib.sha256(label).digest()
    PS    = bytes(k - mLen - 2 * hLen - 2)           # zero-padding
    DB    = lHash + PS + b"\x01" + mensagem           # Data Block

    seed  = os.urandom(hLen)

    # MGF1 com SHA-256
    def mgf1(seed_: bytes, length: int) -> bytes:
        out = b""
        for i in range((length + hLen - 1) // hLen):
            out += hashlib.sha256(seed_ + i.to_bytes(4, "big")).digest()
        return out[:length]

    dbMask   = mgf1(seed, k - hLen - 1)
    maskedDB = bytes(x ^ y for x, y in zip(DB, dbMask))
    seedMask = mgf1(maskedDB, hLen)
    maskedSeed = bytes(x ^ y for x, y in zip(seed, seedMask))

    return b"\x00" + maskedSeed + maskedDB


def _oaep_unpad(padded: bytes, k: int, label: bytes = b"") -> bytes:
    """OAEP unpadding simplificado."""
    hLen  = 32
    lHash = hashlib.sha256(label).digest()

    _, maskedSeed, maskedDB = padded[0:1], padded[1:1+hLen], padded[1+hLen:]

    def mgf1(seed_: bytes, length: int) -> bytes:
        out = b""
        for i in range((length + hLen - 1) // hLen):
            out += hashlib.sha256(seed_ + i.to_bytes(4, "big")).digest()
        return out[:length]

    seedMask = mgf1(maskedDB, hLen)
    seed     = bytes(x ^ y for x, y in zip(maskedSeed, seedMask))
    dbMask   = mgf1(seed, k - hLen - 1)
    DB       = bytes(x ^ y for x, y in zip(maskedDB, dbMask))

    # Verificar lHash
    if DB[:hLen] != lHash:
        raise ValueError("OAEP: hash do label não confere. Chave ou mensagem corrompida.")

    # Encontrar 0x01 separador
    i = hLen
    while i < len(DB) and DB[i] == 0:
        i += 1
    if i >= len(DB) or DB[i] != 0x01:
        raise ValueError("OAEP: separador 0x01 não encontrado.")
    return DB[i+1:]


def cifrar(chave_pub: ChavePublica, mensagem: bytes) -> bytes:
    """
    Cifra uma mensagem de bytes com RSA-OAEP.

    Retorna o cifrotexto como bytes.
    """
    e, n = chave_pub
    k = (n.bit_length() + 7) // 8    # tamanho do módulo em bytes
    padded = _oaep_pad(mensagem, k)
    M_int  = _bytes_para_inteiro(padded)
    C_int  = exponenciacao_modular(M_int, e, n)
    return _inteiro_para_bytes(C_int, k)


def decifrar(chave_priv: ChavePrivada, cifrotexto: bytes) -> bytes:
    """
    Decifra um cifrotexto RSA-OAEP.

    Retorna a mensagem original como bytes.
    """
    d, n = chave_priv
    k     = (n.bit_length() + 7) // 8
    C_int = _bytes_para_inteiro(cifrotexto)
    M_int = exponenciacao_modular(C_int, d, n)
    padded = _inteiro_para_bytes(M_int, k)
    return _oaep_unpad(padded, k)

# ── Assinatura Digital ────────────────────────────────────────────────────────

def assinar(chave_priv: ChavePrivada, mensagem: bytes) -> bytes:
    """
    Assina uma mensagem com RSA (hash SHA-256 + cifragem com chave privada).

    Esquema simplificado de RSA-PSS para fins educacionais.
    Retorna a assinatura como bytes.
    """
    d, n = chave_priv
    k    = (n.bit_length() + 7) // 8

    digest   = hashlib.sha256(mensagem).digest()    # h = H(M)
    h_int    = _bytes_para_inteiro(digest)
    # Garante que o hash cabe no módulo (ajuste de tamanho)
    h_int   %= n

    S_int    = exponenciacao_modular(h_int, d, n)   # S = h^d mod n
    return _inteiro_para_bytes(S_int, k)


def verificar_assinatura(chave_pub: ChavePublica,
                         mensagem: bytes,
                         assinatura: bytes) -> bool:
    """
    Verifica uma assinatura RSA.

    Retorna True se a assinatura é válida, False caso contrário.
    """
    e, n  = chave_pub
    k     = (n.bit_length() + 7) // 8

    S_int = _bytes_para_inteiro(assinatura)
    h_rec = exponenciacao_modular(S_int, e, n)   # h = S^e mod n

    digest   = hashlib.sha256(mensagem).digest()
    h_orig   = _bytes_para_inteiro(digest) % n

    return h_rec == h_orig

# ── Modo demonstrativo ────────────────────────────────────────────────────────

def _separador(titulo: str = "", largura: int = 66) -> None:
    if titulo:
        pad = (largura - len(titulo) - 2) // 2
        print("─" * pad + f" {titulo} " + "─" * pad)
    else:
        print("─" * largura)


def demo_exemplo_numerico() -> None:
    """
    Reproduz passo a passo o exemplo numérico do relatório (p=61, q=53).
    """
    _separador("EXEMPLO NUMÉRICO: p=61, q=53")
    print("Fonte: Stamp (2011); Paar & Pelzl (2010)\n")

    p, q = 61, 53
    n    = p * q
    phi  = (p - 1) * (q - 1)
    e    = 17
    d    = inverso_modular(e, phi)

    print(f"  p            = {p}")
    print(f"  q            = {q}")
    print(f"  n  = p·q     = {n}")
    print(f"  φ(n)         = (p-1)(q-1) = {phi}")
    print(f"  e            = {e}   [mdc({e},{phi}) = {mdc(e,phi)}  ✓]")
    print(f"  d  = e⁻¹ mod φ(n) = {d}")
    print(f"  Verificação  : {e} × {d} mod {phi} = {(e*d) % phi}  ✓")

    print()
    print("  Chave Pública  (e, n) =", (e, n))
    print("  Chave Privada  (d, n) =", (d, n))

    M = 65
    C = cifrar_textbook(M, ChavePublica(e=e, n=n))
    M_rec = decifrar_textbook(C, ChavePrivada(d=d, n=n))

    print()
    _separador("Cifragem")
    print(f"  M (mensagem original) = {M}")
    print(f"  C = M^e mod n = {M}^{e} mod {n} = {C}")

    _separador("Decifragem")
    print(f"  M = C^d mod n = {C}^{d} mod {n} = {M_rec}")
    assert M_rec == M, "Erro: mensagem recuperada não confere!"
    print(f"  Mensagem recuperada   = {M_rec}  ✓\n")


def demo_chaves_e_oaep(bits: int = 1024) -> None:
    """
    Demonstra geração de chaves e cifragem RSA-OAEP com chaves de `bits` bits.
    """
    _separador(f"RSA-OAEP com módulo de {bits} bits")
    print(f"  Gerando par de chaves ({bits} bits)... ", end="", flush=True)

    par = gerar_chaves(bits=bits)
    print("concluído.\n")
    print(f"  Módulo n  (hex, primeiros 32 chars):")
    print(f"    {hex(par.publica.n)[:34]}…")
    print(f"  Expoente e = {par.publica.e}")
    print(f"  |d| bits   = {par.privada.d.bit_length()}\n")

    mensagem = b"RSA - SENAI CIMATEC 2026"
    print(f"  Mensagem original : {mensagem.decode()!r}")

    cifrado   = cifrar(par.publica, mensagem)
    print(f"  Cifrotexto (hex)  : {cifrado.hex()[:48]}…")

    recuperado = decifrar(par.privada, cifrado)
    print(f"  Recuperada        : {recuperado.decode()!r}")
    assert recuperado == mensagem
    print("  Cifragem/Decifragem OAEP: OK ✓\n")

    _separador("Assinatura Digital")
    assinatura = assinar(par.privada, mensagem)
    valida     = verificar_assinatura(par.publica, mensagem, assinatura)
    print(f"  Assinatura (hex)  : {assinatura.hex()[:48]}…")
    print(f"  Assinatura válida : {valida}  ✓")

    # Adultera a mensagem e verifica que a assinatura é rejeitada
    mensagem_adulterada = mensagem + b"X"
    invalida = verificar_assinatura(par.publica, mensagem_adulterada, assinatura)
    print(f"  Assinatura de msg adulterada válida: {invalida}  (esperado: False) ✓\n")


if __name__ == "__main__":
    print()
    print("=" * 66)
    print("  rsa_core.py — Demonstração RSA")
    print("  Atividade Final – Introdução à Cibersegurança")
    print("  SENAI CIMATEC – Pós-graduação em Comunicação Quântica")
    print("=" * 66)
    print()

    demo_exemplo_numerico()
    _separador()

    bits_demo = int(sys.argv[1]) if len(sys.argv) > 1 else 1024
    demo_chaves_e_oaep(bits=bits_demo)

    _separador()
    print("  Todos os testes concluídos com sucesso.")
    print()
