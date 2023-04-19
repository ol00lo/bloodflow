file ./bloodflow_test

set args [upwind]

set breakpoint pending on
b transport_test.cpp:93
run

