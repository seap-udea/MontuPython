# MontuPython — Guía del desarrollador

Este documento describe los comandos del `makefile` del repositorio: entorno local, pruebas, documentación, releases del paquete Python (`montu`) y releases de la aplicación de escritorio (`montu_gui`).

Para ver un resumen rápido en terminal:

```bash
make help
```

---

## Variables del Makefile

| Variable | Valor por defecto | Uso |
|----------|-------------------|-----|
| `PYTHON` | `python3` | Intérprete Python |
| `PIP` | `pip3` | Instalador de paquetes |
| `COMMIT_MSG` | `[MAN] Maintenance` | Mensaje para `make commit` / `make push` |
| `RELMODE` | `release` | Modo de release del paquete: `release` (PyPI) o `test` (TestPyPI) |
| `VERSION` | — | Versión nueva (obligatoria en releases) |
| `TAG` | — | Tag de Git para disparar CI del Desktop (p. ej. `desktop-v0.1.2`) |
| `SPHINXOPTS` | vacío | Opciones extra para Sphinx |

Versiones que el Makefile lee automáticamente:

| Componente | Archivo |
|------------|---------|
| Paquete `montu` | `setup.py`, `montu/version.py`, historial en `.versions` |
| Desktop | `montu_gui/version.py`, historial en `.desktop-versions` |

El paquete y la app de escritorio tienen **versiones independientes**.

---

## Entorno y paquete Python

| Comando | Descripción |
|---------|-------------|
| `make env` | Crea el entorno virtual `.montuenv` e instala `montu` en modo editable |
| `make install` | Instala el paquete (`pip install .`) |
| `make install-dev` | Instala en modo editable (`pip install -e .`) + `requirements.txt` |
| `make build` | Limpia artefactos viejos y construye sdist/wheel en `dist/` (`python -m build`) |
| `make import` | Verifica la instalación importando `montu` e imprimiendo la versión |
| `make show` | Muestra versión (desde `.versions`) y rama Git actual |
| `make status` | Ejecuta `git status` |

Activar el entorno de desarrollo:

```bash
make env
source .montuenv/bin/activate
```

Lanzar la app de escritorio en desarrollo:

```bash
./bin/montu-gui
./bin/montu-gui --debug
```

---

## Limpieza

| Comando | Descripción |
|---------|-------------|
| `make clean` | Elimina archivos temporales, `__pycache__`, `.DS_Store`, etc. |
| `make cleanall` | `clean` + objetos compilados + `dist/`, `build/`, artefactos del Desktop |
| `make cleancrap` | Solo basura de editor y sistema |
| `make cleanout` | Objetos compilados y checkpoints de Jupyter |
| `make cleandist` | `dist/`, `build/`, bundles del Desktop |

---

## Git

| Comando | Descripción |
|---------|-------------|
| `make addall` | `cleanall` + `git add -A` |
| `make commit` | Commit con `COMMIT_MSG` y push a la rama actual |
| `make pull` | `git pull` de la rama actual |
| `make push` | Commit (si hay cambios) y `git push -u origin HEAD` |

Ejemplo:

```bash
make push COMMIT_MSG="Fix calendar conversion for BCE dates"
```

---

## Pruebas

| Comando | Descripción |
|---------|-------------|
| `make test-install` | Instala dependencias de `requirements-test.txt` |
| `make test` | Suite completa de `pytest` |
| `make test-docstrings` | Tests derivados de ejemplos en docstrings |
| `make test-notebooks` | Tests derivados de notebooks de ejemplo |
| `make test-structure` | Valida la estructura de los notebooks de ejemplo |

---

## Documentación

| Comando | Descripción |
|---------|-------------|
| `make docs-install` | Instala dependencias de Sphinx (`docs/requirements.txt`) |
| `make docs-prepare` | Copia notebooks de ejemplo a `docs/examples/` |
| `make docs-build` | Genera HTML en `docs/_build/html/` |
| `make docs` | `docs-prepare` + `docs-build` |
| `make docs-clean` | Borra `docs/_build/` |
| `make readme` | Convierte `README.ipynb` a `README.md` |

---

## Release del paquete Python (`montu` → PyPI)

El release del **paquete** publica en PyPI (o TestPyPI) la librería instalable con `pip install montu`. No incluye el empaquetado de la app de escritorio.

### Comando principal

```bash
make release RELMODE=release VERSION=0.21.0
```

### Modos

| `RELMODE` | Destino |
|-----------|---------|
| `release` | [PyPI](https://pypi.org/project/montu/) |
| `test` | [TestPyPI](https://test.pypi.org/) |

### Qué hace `bin/release.sh`

1. Exige un working tree limpio (sin cambios sin commitear).
2. Actualiza la versión en `setup.py` y `montu/version.py`.
3. Añade la versión a `.versions`.
4. Construye sdist y wheel (`python -m build`).
5. Valida con `twine check`.
6. Sube con `twine upload`.
7. Si algo falla, **revierte** automáticamente los archivos de versión.

### Dry-run (sin subir a PyPI)

```bash
bash bin/release.sh release 0.21.0 --dry-run
```

### Flujo recomendado

```bash
make test
make release RELMODE=test VERSION=0.21.0      # probar en TestPyPI
make release RELMODE=release VERSION=0.21.0   # publicar en PyPI
git add setup.py montu/version.py .versions WHATSNEW.md
git commit -m "Release montu 0.21.0"
git push origin main
```

Requisitos: credenciales de PyPI configuradas para `twine` (p. ej. `~/.pypirc` o variables de entorno).

---

## Release de MontuPython Desktop (app de escritorio)

El release del **Desktop** genera ejecutables para macOS y Windows. La versión vive en `montu_gui/version.py` y es independiente de la del paquete.

### Comandos del Makefile

| Comando | Descripción |
|---------|-------------|
| `make desktop-show` | Muestra la versión actual del Desktop |
| `make desktop-install-build` | Prepara el venv `.desktop-build` con PyInstaller y dependencias |
| `make desktop-build` | Construye la app (`.app` en macOS, carpeta onedir en Windows/Linux) |
| `make desktop-package` | Build + empaquetado (`.zip`/`.dmg` en Mac, `.zip` en Windows) |
| `make desktop-clean` | Borra artefactos del Desktop en `dist/` y `build/` |
| `make desktop-release VERSION=x.y.z` | Sube versión en `montu_gui/version.py` y construye localmente |
| `make desktop-ci TAG=desktop-vx.y.z` | Empuja un tag a GitHub para disparar CI |

### Build local

```bash
make desktop-package
```

Salida en macOS:

```
dist/MontuPython Desktop.app
dist/desktop/MontuPython-Desktop-<versión>-<fecha>-macos.zip
dist/desktop/MontuPython-Desktop-<versión>-<fecha>-macos.dmg
```

En Windows (PowerShell):

```powershell
.\bin\build-desktop.ps1
```

Salida:

```
dist/MontuPython-Desktop/MontuPython-Desktop.exe
dist/desktop/MontuPython-Desktop-<versión>-<fecha>-windows.zip
```

> **Nota:** el `.dmg`/`.zip` de macOS no sirve en Windows. Hay que generar (o descargar del CI) el artefacto `*-windows.zip`.

### Release local completo

```bash
make desktop-release VERSION=0.1.2
```

Equivale a:

1. Actualizar `montu_gui/version.py` y `.desktop-versions`.
2. Ejecutar PyInstaller (`bin/build-desktop.sh --clean`).
3. Crear archivos en `dist/desktop/`.

Opciones del script:

```bash
./bin/desktop-release.sh 0.1.2 --no-bump   # reconstruir sin cambiar versión
./bin/desktop-release.sh 0.1.2 --dry-run # solo actualizar versión, sin build
```

### Release con CI en GitHub (Mac + Windows)

El workflow `.github/workflows/desktop-release.yml` construye en paralelo en `macos-latest` y `windows-latest`.

**Disparadores:**

| Evento | Resultado |
|--------|-----------|
| Push de tag `desktop-v*` | Build Mac + Windows + [GitHub Release](https://github.com/seap-udea/MontuPython/releases) con archivos adjuntos |
| *Actions → Desktop release → Run workflow* | Build Mac + Windows; artefactos en la pestaña *Artifacts* |

**Flujo recomendado:**

```bash
# 1. Actualizar versión y notas
make desktop-release VERSION=0.1.2   # o editar montu_gui/version.py a mano

# 2. Commit y push
git add montu_gui/version.py montu_gui/WHATSNEW.md .desktop-versions
git commit -m "Desktop 0.1.2"
git push origin main

# 3. Tag y CI
git tag desktop-v0.1.2
make desktop-ci TAG=desktop-v0.1.2
```

Cuando termine el workflow (~15–25 min), descarga los instaladores desde **Releases** o desde **Actions → Artifacts**.

### Qué entregar a los usuarios finales

| Sistema | Archivo |
|---------|---------|
| macOS | `.dmg` o `.zip` (`*-macos.*`) |
| Windows | `.zip` (`*-windows.zip`) → ejecutar `MontuPython-Desktop.exe` |

Sin firma de código, macOS puede pedir clic derecho → **Abrir**; Windows puede mostrar advertencia de SmartScreen.

---

## Resumen: paquete vs. Desktop

| | Paquete `montu` | Desktop `montu_gui` |
|--|-----------------|---------------------|
| **Versión** | `setup.py` / `montu/version.py` | `montu_gui/version.py` |
| **Comando release** | `make release VERSION=x.y.z` | `make desktop-release VERSION=x.y.z` |
| **CI** | — (release manual a PyPI) | Tag `desktop-v*` → GitHub Actions |
| **Salida** | Wheel/sdist en PyPI | `.app`, `.dmg`, `.zip`, `.exe` |
| **Instalación usuario** | `pip install montu` | Ejecutar app empaquetada |

---

## Archivos relevantes

| Archivo | Rol |
|---------|-----|
| `makefile` | Punto de entrada de todos los comandos |
| `bin/release.sh` | Release del paquete a PyPI |
| `bin/build-desktop.sh` | Build local del Desktop (macOS/Linux) |
| `bin/build-desktop.ps1` | Build local del Desktop (Windows) |
| `bin/desktop-release.sh` | Bump de versión + build del Desktop |
| `montu_gui/montu-desktop.spec` | Configuración de PyInstaller |
| `requirements-desktop-build.txt` | PyInstaller y hooks |
| `.github/workflows/desktop-release.yml` | CI Mac + Windows |
