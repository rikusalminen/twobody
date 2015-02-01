#ifndef NUMTEST_H
#define NUMTEST_H

#include <stdint.h>

#define ZEROF(x) ((x)*(x) < 1.0e-15)
#define EQF(a, b) ((ZEROF(a) && ZEROF(b)) || ZEROF(((a)-(b))*((a)-(b))/((a)*(a) + ((b)*(b)))))
#define LTF(a, b) ((a) < (b) || EQF((a), (b)))

#define ASSERT(cond, msg, ...) \
    do { \
        numtest_assert( \
            (cond), test_ctx, __FILE__, __LINE__, __FUNCTION__, \
            (msg), ##__VA_ARGS__); \
    } while(0)
#define ASSERT_EQF(a, b, msg, ...) ASSERT(EQF((a), (b)), msg, ##__VA_ARGS__)
#define ASSERT_LTF(a, b, msg, ...) ASSERT(LTF((a), (b)), msg, ##__VA_ARGS__)
#define ASSERT_RANGEF(x, min, max, msg, ...) ASSERT(LTF((min), (x)) && LTF((x), (max)), msg, ##__VA_ARGS__)

struct numtest_ctx;

typedef void (numtest_callback)(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx);

struct numtest_case {
    const char *name;
    numtest_callback *func;
    int num_params;
    void *extra_args;
};

extern uint64_t numtest_num_cases_default;
extern const struct numtest_case numtest_cases[];

void numtest_assert(
    int cond,
    struct numtest_ctx *ctx,
    const char *file, int line, const char *function,
    const char *msg, ...);

int numtest_main(int argc, char *argv[]);

#endif
