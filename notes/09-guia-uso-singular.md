# Guía práctica: usar Singular en este proyecto (vía WSL)

Documento de referencia operativa — no narra un hallazgo, sino cómo usar la herramienta día a día.
Para el porqué de la migración y los resultados encontrados, ver `07-migracion-a-singular.md` y
`08-anillo-completo-via-singular.md`.

## Instalación (ya hecha, dejar constancia por si hay que reinstalar)

Singular vive dentro de WSL (distribución `Ubuntu-24.04`), NO como binario nativo de Windows:

```bash
wsl -d Ubuntu-24.04 -- sudo apt update
wsl -d Ubuntu-24.04 -- sudo apt install -y singular
```

`sudo` pide contraseña interactiva — Claude no puede escribirla, así que este paso lo tiene que
correr el usuario (con el prefijo `!` en el prompt de Claude Code). Versión instalada: `4.3.2`
(paquete de Ubuntu 24.04, no la última versión upstream pero suficiente). Tamaño real: ~9-10 MB de
descarga, ~30-32 MB instalado (ver `07-migracion-a-singular.md` para el desglose completo).

Verificar que quedó bien instalado:
```bash
wsl -d Ubuntu-24.04 -- Singular --version
```

## Cómo correr un script desde Windows/Claude Code

Los archivos del proyecto viven en Windows (`C:\Users\marti\Documents\repositorios\cuerdas\...`),
pero Singular corre dentro de Linux (WSL). WSL monta el disco de Windows en `/mnt/c/...`, así que
se puede acceder a los mismos archivos sin copiarlos:

```bash
wsl -d Ubuntu-24.04 -- bash -c "cd '/mnt/c/Users/marti/Documents/repositorios/cuerdas/fase2-objetivo-tecnico/scripts' && Singular -q nombre_del_script.sing"
```

La flag `-q` suprime el banner de bienvenida (que si no, contamina la salida). Los scripts de
Singular en este proyecto usan extensión `.sing` por convención (no es obligatorio, Singular no
exige extensión).

## Lección aprendida esta sesión: separar cálculo completo de resumen

Un script que calcula el anillo de invariantes completo (`invariant_ring(...)`) y hace `print()` de
todos los resultados puede generar una salida de **decenas de miles de caracteres por línea** (un
solo invariante secundario de grado alto puede tener miles de términos). Si esa salida se lee
directamente con la herramienta `Read` de Claude, puede consumir la ventana de contexto completa de
un tirón (~110.000 tokens en un caso real de esta sesión).

**Patrón recomendado**: escribir SIEMPRE 2 versiones de cada script de cómputo pesado:
1. Un script "completo" (ej. `s5_invariant_ring.sing`) que calcula y guarda todo — correrlo con la
   salida redirigida a un archivo (`> salida.txt`), NUNCA leído directamente con `Read`.
2. Un script "resumen" (ej. `s5_invariant_ring_summary.sing`) que solo imprime metadatos: cuántos
   invariantes hay, de qué grado (`deg(...)`), cuántos términos tiene cada uno (`size(...)`) — sin
   imprimir los polinomios en sí. Este es el que se lee directamente.

Para inspeccionar la salida completa sin cargarla entera, usar `Bash` con `wc`, `awk`, `grep`,
`head -c`/`tail -c` (nunca la herramienta `Read` sobre un archivo de ese tamaño).

## Sintaxis mínima usada en este proyecto (referencia rápida)

```singular
LIB "finvar.lib";                          // cargar la libreria de anillos de invariantes

ring R = 0, (a1,a2,a3,a4,a5), dp;          // cuerpo de caracteristica 0 (racionales),
                                             // 5 variables, orden "dp" (grado, luego lexicografico)

matrix G1[5][5] = <15 entradas por fila>;   // matriz generadora del grupo (entera o racional)
matrix G2[5][5] = ...;

list result = invariant_ring(G1, G2);       // G1, G2 deben generar un grupo FINITO
matrix primaries = result[1];               // invariantes primarios (Noether normalization)
matrix secondaries = result[2];             // invariantes secundarios

deg(primaries[1,i]);                        // grado del i-esimo invariante primario
size(secondaries[1,i]);                     // numero de terminos del i-esimo secundario
```

`invariant_ring` requiere que las matrices generen un grupo **finito** — no lo verifica por sí
mismo; si se le dan matrices que generan un grupo infinito, corre indefinidamente. Por eso conviene
verificar antes (como se hizo en Python) que el grupo generado tiene el orden esperado.

## Dónde está todo en este proyecto

- `scripts/s5_invariant_ring.sing` — cálculo completo (genera los polinomios enteros).
- `scripts/s5_invariant_ring_summary.sing` — resumen (grados y tamaños).
- `scripts/s5_invariant_ring_output.txt` — salida guardada del cálculo completo.
- `07-migracion-a-singular.md` — por qué se migró, tamaños de instalación.
- `08-anillo-completo-via-singular.md` — el resultado matemático encontrado y su verificación.
- Memoria de Claude: `reference_singular_via_wsl.md` (para que otras sesiones/proyectos sepan que
  Singular está disponible sin tener que re-investigar la instalación).
