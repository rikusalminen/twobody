#include <twobody/twobody.h>
#include "../numtest.h"

extern numtest_callback
    dummy_test;

const struct numtest_case numtest_cases[] = {
    { 0, 0, 0, 0 }
    };

int main(int argc, char *argv[]) {
    twobody_version();
    return numtest_main(argc, argv);
}

