; Inno Setup 6.3+ script (https://jrsoftware.org/isinfo.php). Open this file in Inno Setup and press Compile.
; It is copied to dist\ by build_portable.R; the app files are in dist\Sangerly\.
#define AppName "Sangerly"
#define AppVer "1.0"

[Setup]
AppId={{B2F5C3E0-6B1A-4E8E-9D7A-5A1C0E5E1A01}
AppName={#AppName}
AppVersion={#AppVer}
DefaultDirName={localappdata}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=Output
OutputBaseFilename=Sangerly_Setup
Compression=lzma2/max
SolidCompression=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"

[Files]
Source: "Sangerly\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{sys}\wscript.exe"; Parameters: """{app}\Sangerly.vbs"""; WorkingDir: "{app}"; IconFilename: "{app}\Sangerly.ico"
Name: "{autodesktop}\{#AppName}"; Filename: "{sys}\wscript.exe"; Parameters: """{app}\Sangerly.vbs"""; WorkingDir: "{app}"; IconFilename: "{app}\Sangerly.ico"; Tasks: desktopicon

[Run]
Filename: "{sys}\wscript.exe"; Parameters: """{app}\Sangerly.vbs"""; Description: "Start Sangerly"; Flags: postinstall nowait skipifsilent
