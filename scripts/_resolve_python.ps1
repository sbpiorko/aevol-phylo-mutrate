function Resolve-ProjectPython {
    if ($env:PYTHON -and (Test-Path $env:PYTHON)) {
        return $env:PYTHON
    }

    foreach ($name in @("python", "python3", "py")) {
        $cmd = Get-Command $name -ErrorAction SilentlyContinue
        if ($cmd) {
            $oldErrorActionPreference = $ErrorActionPreference
            $ErrorActionPreference = "Continue"
            try {
                & $cmd.Source --version *> $null
                $exitCode = $LASTEXITCODE
            }
            catch {
                $exitCode = 1
            }
            finally {
                $ErrorActionPreference = $oldErrorActionPreference
            }
            if ($exitCode -eq 0) {
                return $cmd.Source
            }
        }
    }

    $bundled = Join-Path $env:USERPROFILE ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
    if (Test-Path $bundled) {
        return $bundled
    }

    throw "Could not find Python. Install Python, set `$env:PYTHON, or run inside Codex with the bundled runtime."
}
