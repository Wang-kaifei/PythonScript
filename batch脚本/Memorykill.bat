SET procName=pAnno.exe
SET RAMLimit=52340
:loop
FOR /F "tokens=*" %%F IN ('tasklist^|findstr %procName%') DO SET foundString=%%F
FOR /F "tokens=5" %%F IN ("%foundString%") DO SET RAMConsumption=%%F
echo %RAMConsumption%
echo RAMConsumption=%RAMConsumption:,=%
set RAMConsumption=%RAMConsumption:,=%
echo %RAMConsumption%
IF %RAMConsumption% GEQ %RAMLimit% TASKKILL /IM %procName%
timeout /T 3 /NOBREAK
goto loop
pause