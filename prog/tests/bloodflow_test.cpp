#define CATCH_CONFIG_RUNNER
#include "testutils.hpp"

TEST_CASE("Ping", "[ping]"){
	REQUIRE(true);
}

int main(int argc, char* argv[]){
	int result = Catch::Session().run(argc, argv);
	std::cout << "DONE" << std::endl;

	return result;
}

