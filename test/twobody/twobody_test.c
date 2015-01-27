#include <twobody/twobody.h>
#include "../numtest.h"

uint64_t numtest_num_cases_default = 1 << 23;

extern numtest_callback
    conic_test,
    dummy_test;

const struct numtest_case numtest_cases[] = {
    { "conic", conic_test, 3, 0 },
    { 0, 0, 0, 0 }
    };

int main(int argc, char *argv[]) {
    twobody_version();
    return numtest_main(argc, argv);
}

