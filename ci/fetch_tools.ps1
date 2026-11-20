# Downloads the external programs into the folders that build_portable.R / build_exe.py look for.
# CHECK the URLs and versions once (the project pages are in THIRD_PARTY_NOTICES.md) and keep them pinned.
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
function Get-Zip($url, $dest) {
  $zip = Join-Path $env:RUNNER_TEMP ([IO.Path]::GetFileName($url))
  Invoke-WebRequest -Uri $url -OutFile $zip
  New-Item -ItemType Directory -Force $dest | Out-Null
  Expand-Archive -Path $zip -DestinationPath $dest -Force
}
# IQ-TREE 3 (GPL-2.0) -> iqtree-3.0.1-Windows\bin\iqtree3.exe
Get-Zip "https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Windows.zip" $root
# trimAl (GPL-3.0) -> trimAL\trimal.exe          (TODO: put the URL of the Windows build you use)
if ($env:TRIMAL_URL) { Get-Zip $env:TRIMAL_URL (Join-Path $root "trimAL") } else { Write-Warning "TRIMAL_URL not set: the build will have no trimAl" }
# minimap2 (MIT) -> minimap2\minimap2.exe        (TODO: optional)
if ($env:MINIMAP2_URL) { Get-Zip $env:MINIMAP2_URL (Join-Path $root "minimap2") } else { Write-Warning "MINIMAP2_URL not set: the build will have no minimap2" }
