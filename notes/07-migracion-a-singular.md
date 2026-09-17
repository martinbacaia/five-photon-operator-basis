# Fase 2 — Decisión: migrar a Singular para el cálculo del anillo de invariantes de S₅

## Por qué

El método artesanal (sympy + operador de Reynolds + adivinar semillas grado por grado, ver
`04-generadores-explicitos-s5.md` y `06-generadores-grado5-6.md`) encontró generadores hasta grado 6,
pero mostró 2 problemas reales a medida que crece el grado:

1. **Velocidad**: los cálculos de grado 6 ya necesitaban varios minutos corriendo en segundo plano
   (sumas de Reynolds sobre 120 elementos, con polinomios de cientos de términos en sympy puro).
   El usuario indicó explícitamente (en esta y otras sesiones) que si Python se vuelve lento para
   este tipo de cómputo, hay que migrar a una herramienta/lenguaje más rápido en vez de seguir
   esperando u optimizando micro-detalles en Python.
2. **Fragilidad**: se encontraron y corrigieron 3 bugs reales en el método artesanal esta misma
   sesión (convención de homomorfismo, mal uso de `nsimplify`, semilla de grado equivocado) — cada
   uno originado en un paso manual (elegir semillas, evaluar en puntos, convertir tipos) que un
   algoritmo de teoría de invariantes ya implementado y probado por la comunidad no necesita.

**Decisión**: migrar el cálculo del anillo de invariantes a **Singular** (Computer Algebra System
especializado en álgebra conmutativa/geometría algebraica, con la librería `finvar.lib` dedicada
específicamente a anillos de invariantes de grupos finitos — algoritmos de Heydtmann,
`invariant_ring` / `invariant_ring_random`). Es estándar de facto en la comunidad para este tipo de
cálculo (más citable en caso de publicar) y está escrito en C (mucho más rápido que sympy puro para
álgebra simbólica pesada).

## Cómo se instaló

Vía WSL (ya se tenía Ubuntu 24.04 instalado), con el paquete de Ubuntu/Debian `singular`:

```bash
wsl -d Ubuntu-24.04 -- sudo apt update
wsl -d Ubuntu-24.04 -- sudo apt install -y singular
```

Se investigó el tamaño real ANTES de instalar (vía `apt show`, sin necesitar privilegios de root):

| paquete | tamaño instalado | descarga |
|---|---|---|
| singular-ui (interfaz de línea de comandos) | 46 KB | 11 KB |
| singular-data (datos compartidos, arquitectura-independiente) | 15.2 MB | 4.9 MB |
| singular-modules (módulos básicos) | 4.95 MB | 1.0 MB |
| libsingular4m3n0t64 (librería núcleo) | 8.0 MB | 3.0 MB |
| **subtotal Singular** | **~28 MB** | **~9 MB** |
| + recomendados (graphviz, 4ti2, normaliz, surf-alggeo, topcom — no esenciales para `finvar.lib`) | ~4.3 MB | ~0.85 MB |
| **total con recomendados** | **~32 MB** | **~9.7 MB** |

Instalación liviana — nada comparable a instalar SageMath completo (que pesa varios GB). Versión
instalada: `1:4.3.2-p10+ds-1.1build1` (empaquetado por Ubuntu/Debian, no la última versión upstream
4.4.1 de enero 2025, pero suficiente para `finvar.lib`).

## Próximo paso concreto (no hecho todavía)

Traducir la representación de S₅ ya calculada (las matrices de 5×5 para los generadores (12) y
(12345), ya verificadas exhaustivamente en `n5_mandelstam_representation.py` y
`reverificacion_generadores.py`) al formato que espera `finvar.lib`, y correr
`invariant_ring`/`invariant_ring_random` para obtener de una sola vez el conjunto completo y mínimo
de generadores, más las relaciones entre ellos (vía bases de Gröbner) — reemplazando el método de
"adivinar semillas" por un cálculo sistemático y confiable.

## Referencia para futuras sesiones

Ver también la memoria de Claude: `feedback_slow_python_switch_language.md` (el usuario pidió
explícitamente cambiar de lenguaje/herramienta cuando Python se vuelve lento para cómputo pesado).
