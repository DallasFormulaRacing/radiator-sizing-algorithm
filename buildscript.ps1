$env:PATH = "$env:PATH;C:\msys64\mingw64\bin"

Write-Host "`n Building..."
make

Write-Host "`n Running simulation..."
Start-Process -FilePath ".\radiator_sim.exe" -Wait
