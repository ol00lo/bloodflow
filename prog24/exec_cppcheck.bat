"C:\Program Files\Cppcheck\cppcheck.exe" --std=c++17 --xml --suppress=cstyleCast --suppress=missingIncludeSystem --enable=all -I src src 2> cppcheck_result.xml
rem python fix_xml_cpp_check.py cppcheck_result.xml
"C:\Program Files\Cppcheck\cppcheckgui.exe" cppcheck_result.xml