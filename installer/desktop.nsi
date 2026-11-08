; Self-contained Windows desktop edition. Build with python build_setup.py.
Unicode True
!include "MUI2.nsh"
!include "x64.nsh"
!include "LogicLib.nsh"
!include "FileFunc.nsh"
!include "WinVer.nsh"
!include "StrFunc.nsh"
${Using:StrFunc} StrRep
${Using:StrFunc} StrStr
Var NoIntegration

!ifndef PAYLOAD
  !error "Build with build_setup.py to supply PAYLOAD and OUTPUT."
!endif
!ifndef EXE_SHA256
  !error "Build with build_setup.py to supply the verified executable SHA256."
!endif
Name "In-Silico Docking Studio"
OutFile "${OUTPUT}"
InstallDir "$LOCALAPPDATA\Programs\InSilicoTeachingSuite"
InstallDirRegKey HKCU "Software\InSilicoTeachingSuite" "InstallDir"
RequestExecutionLevel user
SetCompressor /SOLID lzma
SetCompressorDictSize 32
BrandingText "Medical education · Learn, explore and interpret"
VIProductVersion "1.3.0.0"
VIAddVersionKey "ProductName" "In-Silico Docking Studio"
VIAddVersionKey "FileDescription" "In-Silico Docking Studio Setup"
VIAddVersionKey "ProductVersion" "1.3.0"
VIAddVersionKey "FileVersion" "1.3.0"
VIAddVersionKey "LegalCopyright" "Dr Manu Kumar Shetty, 2026"

!define MUI_ICON "assets\icon.ico"
!define MUI_UNICON "assets\icon.ico"
!define MUI_ABORTWARNING
!define MUI_WELCOMEPAGE_TITLE "Your step-by-step docking workspace"
!define MUI_WELCOMEPAGE_TEXT "Run studies with multiple protein targets and drug molecules. Follow simple instructions, prepare each binding site and review docking results.$\r$\n$\r$\nSetup includes the scientific tools, logo, references and offline 3D viewer. You do not need Python or Anaconda.$\r$\n$\r$\nFor Windows 10/11 (64-bit). Microsoft Edge WebView2 Runtime is needed for the desktop window."
!define MUI_FINISHPAGE_TITLE "Ready to learn"
!define MUI_FINISHPAGE_TEXT "Open the docking studio from your desktop or Start menu. Create a study, add targets and drugs, then prepare and run docking.$\r$\n$\r$\nYour studies are saved separately in your user profile and are kept when the app is uninstalled."
!define MUI_FINISHPAGE_RUN "$INSTDIR\InSilicoTeachingSuite.exe"
!define MUI_FINISHPAGE_RUN_TEXT "Open In-Silico Docking Studio"
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\LICENSE"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

Function .onInit
  ${IfNot} ${RunningX64}
    MessageBox MB_ICONSTOP "This version needs 64-bit Windows."
    Abort
  ${EndIf}
  ${IfNot} ${AtLeastWin10}
    MessageBox MB_ICONSTOP "This version needs Windows 10 or Windows 11."
    Abort
  ${EndIf}
  SetShellVarContext current
  StrCpy $NoIntegration "0"
  ${GetParameters} $0
  ClearErrors
  ${GetOptions} $0 "/NOINTEGRATION" $1
  ${IfNot} ${Errors}
    StrCpy $NoIntegration "1"
  ${EndIf}
FunctionEnd

Section "Teaching Suite" MainSection
  ; Check the final chosen folder before extracting anything. FileOpen(a)
  ; opens read/write without truncating or writing; Windows refuses this
  ; while the executable is mapped by a running desktop or server process.
  IfFileExists "$INSTDIR\InSilicoTeachingSuite.exe" 0 upgrade_unlocked
  ClearErrors
  FileOpen $0 "$INSTDIR\InSilicoTeachingSuite.exe" a
  IfErrors upgrade_locked
  FileClose $0
  Goto upgrade_unlocked
  upgrade_locked:
    DetailPrint "Close the running In-Silico app before updating. No files were copied."
    IfSilent upgrade_locked_exit
    MessageBox MB_OK|MB_ICONEXCLAMATION "Close the running In-Silico app, then run setup again. No application files have been changed."
  upgrade_locked_exit:
    SetErrorLevel 10
    Quit
  upgrade_unlocked:
  SetOutPath "$INSTDIR"
  ClearErrors
  File /r "${PAYLOAD}\*.*"
  IfErrors extraction_failed
  ; Verify the entire copied executable before updating any installed
  ; version, shortcut or uninstaller metadata. Invoke certutil directly.
  nsExec::ExecToStack /TIMEOUT=30000 '"$SYSDIR\certutil.exe" -hashfile "$INSTDIR\InSilicoTeachingSuite.exe" SHA256'
  Pop $0
  Pop $1
  StrCmp $0 "0" 0 verification_failed
  ; Normalise CRLF, then match a complete hash line, not a path/heading.
  ${StrRep} $2 "$1" "$\r" ""
  ${StrStr} $3 "$2" "$\n${EXE_SHA256}$\n"
  StrCmp $3 "" verification_failed
  Goto payload_verified
  extraction_failed:
    DetailPrint "Setup could not replace every application file. Version metadata was not updated."
    IfSilent extraction_failed_exit
    MessageBox MB_OK|MB_ICONSTOP "Setup could not replace every application file. Close the app and run setup again. Installed version information has not been updated."
  extraction_failed_exit:
    SetErrorLevel 11
    Quit
  verification_failed:
    DetailPrint "The copied application failed SHA256 verification. Version metadata was not updated."
    IfSilent verification_failed_exit
    MessageBox MB_OK|MB_ICONSTOP "The copied application could not be verified. Run setup again before opening the app. Installed version information has not been updated."
  verification_failed_exit:
    SetErrorLevel 12
    Quit
  payload_verified:
  ClearErrors
  WriteUninstaller "$INSTDIR\Uninstall.exe"
  WriteINIStr "$INSTDIR\installation.ini" "Setup" "NoIntegration" "$NoIntegration"
  ${If} $NoIntegration == "1"
    Goto installed
  ${EndIf}
  CreateDirectory "$SMPROGRAMS\In-Silico Teaching Suite"
  CreateShortcut "$SMPROGRAMS\In-Silico Teaching Suite\In-Silico Teaching Suite.lnk" "$INSTDIR\InSilicoTeachingSuite.exe"
  CreateShortcut "$DESKTOP\In-Silico Teaching Suite.lnk" "$INSTDIR\InSilicoTeachingSuite.exe"
  CreateShortcut "$SMPROGRAMS\In-Silico Teaching Suite\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  WriteRegStr HKCU "Software\InSilicoTeachingSuite" "InstallDir" "$INSTDIR"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite" "DisplayName" "In-Silico Docking Studio"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite" "DisplayVersion" "1.3.0"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite" "Publisher" "Dr Manu Kumar Shetty and Dr Pratibha Nadig"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite" "UninstallString" '$\"$INSTDIR\Uninstall.exe$\"'
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite" "DisplayIcon" "$INSTDIR\InSilicoTeachingSuite.exe"
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite" "NoModify" 1
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite" "NoRepair" 1
  installed:
SectionEnd

Section "Uninstall"
  SetShellVarContext current
  ReadINIStr $NoIntegration "$INSTDIR\installation.ini" "Setup" "NoIntegration"
  ; Generated explicit file list: never recursively delete a user directory.
  !include "${UNINSTALL_FILES}"
  Delete "$INSTDIR\Uninstall.exe"
  Delete "$INSTDIR\installation.ini"
  RMDir "$INSTDIR"
  ${If} $NoIntegration == "1"
    Goto uninstalled
  ${EndIf}
  Delete "$DESKTOP\In-Silico Teaching Suite.lnk"
  Delete "$SMPROGRAMS\In-Silico Teaching Suite\In-Silico Teaching Suite.lnk"
  Delete "$SMPROGRAMS\In-Silico Teaching Suite\Uninstall.lnk"
  RMDir "$SMPROGRAMS\In-Silico Teaching Suite"
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\InSilicoTeachingSuite"
  DeleteRegKey HKCU "Software\InSilicoTeachingSuite"
  uninstalled:
SectionEnd
