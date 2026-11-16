Option Explicit
Dim sh, fso, here, rscript, rc, logf
Set sh = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
here = fso.GetParentFolderName(WScript.ScriptFullName)
rscript = here & "\R\bin\Rscript.exe"
logf = sh.ExpandEnvironmentStrings("%LOCALAPPDATA%") & "\Sangerly\last_error.txt"
If Not fso.FileExists(rscript) Then
  MsgBox "R was not found in " & here & "\R. Reinstall Sangerly.", 16, "Sangerly"
  WScript.Quit 1
End If
sh.CurrentDirectory = here
sh.Environment("Process")("R_LIBS") = ""
sh.Environment("Process")("R_LIBS_USER") = ""
sh.Environment("Process")("SANGERLY_DESKTOP") = "1"
' window style 0 = hidden; True = wait until the app has stopped
rc = sh.Run("""" & rscript & """ --vanilla launch.R", 0, True)
If rc <> 0 Then
  If fso.FileExists(logf) Then
    MsgBox "Sangerly stopped unexpectedly (code " & rc & ")." & vbCrLf & "Details: " & logf, 16, "Sangerly"
  Else
    MsgBox "Sangerly stopped before it could write a log (code " & rc & ")." & vbCrLf & "Run Sangerly_debug.bat in the same folder to see the error.", 16, "Sangerly"
  End If
End If
