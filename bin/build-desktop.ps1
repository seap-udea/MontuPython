# Build MontuPython Desktop distributable with PyInstaller (Windows).
#
# Usage:
#   powershell -ExecutionPolicy Bypass -File .\bin\build-desktop.ps1
#   .\bin\build-desktop.ps1 -Clean
#   .\bin\build-desktop.ps1 -NoPackage
#
param(
    [switch]$Clean,
    [switch]$NoPackage
)

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $Root

$Venv = Join-Path $Root ".desktop-build"
$Spec = Join-Path $Root "montu_gui\montu-desktop.spec"
$Dist = Join-Path $Root "dist"
$Build = Join-Path $Root "build"

function Get-DesktopVersion {
    $text = Get-Content (Join-Path $Root "montu_gui\version.py") -Raw
    if ($text -match "version\s*=\s*['""]([^'""]+)['""]") {
        return $Matches[1]
    }
    return "unknown"
}

$DesktopVersion = Get-DesktopVersion
Write-Host "MontuPython Desktop build"
Write-Host "  version : $DesktopVersion"
Write-Host "  platform: Windows"
Write-Host "  root    : $Root"

if (-not (Test-Path $Venv)) {
    Write-Host "Creating build venv at .desktop-build ..."
    python -m venv $Venv
}

$Activate = Join-Path $Venv "Scripts\Activate.ps1"
. $Activate

Write-Host "Installing/upgrading build dependencies ..."
python -m pip install -q --upgrade pip
python -m pip install -q -e .
python -m pip install -q -r requirements.txt
python -m pip install -q -r requirements-desktop-build.txt

if ($Clean) {
    Write-Host "Cleaning previous desktop build artifacts ..."
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue `
        $Dist, $Build, (Join-Path $Dist "MontuPython-Desktop"), (Join-Path $Dist "desktop")
}

Write-Host "Running PyInstaller ..."
pyinstaller $Spec --noconfirm

$OneDir = Join-Path $Dist "MontuPython-Desktop"
$Exe = Join-Path $OneDir "MontuPython-Desktop.exe"
if (-not (Test-Path $Exe)) {
    throw "Expected executable not found: $Exe"
}

Write-Host "Built: $Exe"

if (-not $NoPackage) {
    $OutDir = Join-Path $Dist "desktop"
    New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
    $Stamp = (Get-Date).ToUniversalTime().ToString("yyyyMMdd")
    $ArchiveBase = "MontuPython-Desktop-$DesktopVersion-$Stamp-windows"
    $Zip = Join-Path $OutDir "$ArchiveBase.zip"
    Write-Host "Packaging zip: $Zip"
    Compress-Archive -Path $OneDir -DestinationPath $Zip -Force
    Write-Host "Packaged: $Zip"
}

Write-Host "Done."
