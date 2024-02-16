"C:\Program Files\Cppcheck\cppcheck.exe" --std=c++17 --xml --suppress=cstyleCast --enable=all -I src src 2> cppcheck_result.xml
"C:\Program Files\Cppcheck\cppcheckgui.exe" cppcheck_result.xml