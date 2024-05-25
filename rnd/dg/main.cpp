#define CATCH_CONFIG_RUNNER
#include "common_test.hpp"

TEST_CASE("bflow ping", "[ping]"){
    CHECK(1 == 1);
}


int main(int argc, char* argv[]){
	int result = Catch::Session().run(argc, argv);
	std::cout << "DONE" << std::endl;
	return result;
}
