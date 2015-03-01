#include <twobody/twobody.h>
#include "../numtest.h"

uint64_t numtest_num_cases_default = 1 << 23;

extern numtest_callback
    conic_test,
    anomaly_test,
    true_anomaly_test,
    eccentric_anomaly_test,
    orientation_test,
    orbit_from_state_test,
    orbit_from_elements_test,
    orbit_radial_test,
    stumpff_test,
    universal_test,
    dummy_test;

const struct numtest_case numtest_cases[] = {
    { "conic", conic_test, 3, 0 },
    { "anomaly", anomaly_test, 2, 0 },
    { "true_anomaly", true_anomaly_test, 4, 0 },
    { "eccentric_anomaly", eccentric_anomaly_test, 4, 0 },
    { "orientation", orientation_test, 3, 0 },
    { "orbit_from_state", orbit_from_state_test, 7, 0 },
    { "orbit_from_elements", orbit_from_elements_test, 6, 0 },
    { "orbit_radial", orbit_radial_test, 5, 0 },
    { "stumpff", stumpff_test, 2, 0 },
    { "universal", universal_test, 5, 0 },
    { 0, 0, 0, 0 }
    };

int main(int argc, char *argv[]) {
    twobody_version();
    return numtest_main(argc, argv);
}

