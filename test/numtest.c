#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "numtest.h"

struct numtest_args {
    uint64_t first, last;
    int random, random_seed;
    int silent, verbose;

    char * const *tests;
    int num_tests;

    int list_tests;
};

struct numtest_ctx {
    const char *test_case_name;
    uint64_t seed;

    int asserts_passed;
    int asserts_failed;

    uint64_t cases_passed;
    uint64_t cases_failed;

    int tests_run;

    const struct numtest_args *args;
};

void numtest_assert_failed(
    struct numtest_ctx *ctx,
    const char *file, int line, const char *function,
    const char *msg,
    va_list va)
{
    uint64_t max_failures = 100;
    if(ctx->args->silent || (!ctx->args->verbose && ctx->cases_failed >= max_failures))
        return;     // only print first few failures to avoid spamming logs

    FILE *out = stdout;
    fprintf(out, "%s(%lu): ASSERT FAILED (%s:%d %s)  \n\t",
            ctx->test_case_name,
            ctx->seed,
            file,
            line,
            function);
    vfprintf(out, msg, va);
    fprintf(out, "\n");
}

void numtest_assert(
    int cond,
    struct numtest_ctx *ctx,
    const char *file, int line, const char *function,
    const char *msg, ...) {

    if(cond) {
        ctx->asserts_passed += 1;
        return;
    }

    ctx->asserts_failed += 1;

    va_list va;
    va_start(va, msg);
    numtest_assert_failed(ctx, file, line, function, msg, va);
    va_end(va);
}

static double test_pattern_1d(uint64_t seed) {
    if(seed < 2)
        return (double)seed;

    seed -= 1;

    int level = 64 - __builtin_clzl(seed);
    uint64_t numer = 1 + (seed - (1ull << (level - 1))) * 2;
    uint64_t denom = 1ull << level;

    return numer / (double)denom;
}

static void test_pattern(uint64_t seed, int dim, double *params) {
    if(dim == 1) {
        params[0] = test_pattern_1d(seed);
        return;
    }

    uint64_t seeds[dim];
    for(int i = 0; i < dim; ++i)
        seeds[i] = 0;

    // interpret seed as Morton code
    for(int i = 0; i < 64 && (1ull << i) <= seed; ++i)
        if(seed & (1ull << i))
            seeds[i % dim] |= 1ull << (i / dim);

    for(int i = 0; i < dim; ++i)
        params[i] = test_pattern_1d(seeds[i]);
}

static bool numtest_run_tests(const struct numtest_args *args) {
    struct numtest_ctx ctx = { 0, 0, 0, 0, 0, 0, 0, 0 };
    ctx.args = args;

    time_t time_begin = time(0), time_output = time_begin;
    uint64_t total_pass = 0, total_fail = 0;

    for(const struct numtest_case *test_case = numtest_cases + 0;
        test_case->name != 0;
        ++test_case) {
        bool skip = args->num_tests != 0;
        for(int i = 0; i < args->num_tests; ++i)
            if(strcmp(test_case->name, args->tests[i]) == 0)
                skip = false;
        if(skip)
            continue;

        ctx.test_case_name = test_case->name;
        ctx.cases_passed = ctx.cases_failed = 0;

        uint64_t num = 1 << 23;
        uint64_t first = args->first;
        uint64_t last = args->last != 0 ? args->last : num;

        time_t time_test_begin = time(0);

        for(uint64_t xxx = first; xxx <= last; ++xxx) {
            uint64_t seed = xxx; // TODO: random and random seed!

            double params[test_case->num_params];
            test_pattern(seed, test_case->num_params, params);

            // TODO: optionally add noise

            ctx.seed = seed;
            ctx.asserts_passed = ctx.asserts_failed = 0;
            test_case->func(params, test_case->num_params, test_case->extra_args, &ctx);

            if(ctx.asserts_failed == 0)
                ctx.cases_passed += 1;
            else
                ctx.cases_failed += 1;

            if(!args->silent && ctx.cases_passed + ctx.cases_failed % 10000) {
                time_t time_now = time(0);

                if(time_now - time_output >= 60) { // write a message once a minute
                    time_t total = time_now - time_begin,
                           seconds = total % 60,
                           minutes = (total / 60) % 60,
                           hours = total / 3600;
                    fprintf(stderr,
                        "%02luh%02lum%02lus %s (%02lu%%, %lu pass, %lu fail)\n",
                        hours, minutes, seconds,
                        ctx.test_case_name, 100*(xxx-first)/(last-first),
                        ctx.cases_passed, ctx.cases_failed);
                    time_output = time_now;
                }
            }
        }

        total_pass += ctx.cases_passed;
        total_fail += ctx.cases_failed;
        ctx.tests_run += 1;

        time_t time_test_end = time(0);
        if(!args->silent) {
            time_t total = time_test_end - time_test_begin,
                   seconds = total % 60,
                   minutes = (total / 60) % 60,
                   hours = total / 3600;
            fprintf(stderr,
                "%s %s  %lu%% (%lu pass, %lu fail, %02luh%02lum%02lus)\n",
                ctx.cases_failed == 0 ? "PASS" : "FAIL", ctx.test_case_name,
                100 * ctx.cases_passed / (ctx.cases_passed + ctx.cases_failed),
                ctx.cases_passed, ctx.cases_failed,
                hours, minutes, seconds);
            time_output = time_test_end;
        }
    }

    time_t time_end = time(0);
    if(!args->silent && total_fail + total_pass) {
        time_t total = time_end - time_begin,
               seconds = total % 60,
               minutes = (total / 60) % 60,
               hours = total / 3600;

        fprintf(stdout,
            "TESTS %s  %lu%%  "
            "(%d tests, %lu cases pass, %lu cases fail, %02luh%02lum%02lus)\n",
            total_fail == 0 ? "PASS" : "FAIL",
            100 * total_pass / (total_pass + total_fail),
            ctx.tests_run, total_pass, total_fail,
            hours, minutes, seconds);
    } else {
        fprintf(stderr, "NO TESTS RUN\n");
        return 1;
    }

    return total_fail == 0;
}

static void numtest_list_tests() {
    for(const struct numtest_case *test_case = numtest_cases + 0;
        test_case->name != 0;
        ++test_case)
        fprintf(stdout, "%s\n", test_case->name);
}

#include <getopt.h>
#include <stdlib.h>
#include <time.h>

static void usage() {
    printf("*** USAGE ***\n");
    exit(EXIT_FAILURE);
}

static struct numtest_args parse_args(int argc, char * const argv[]) {
    const struct option long_options[] = {
        {"first", required_argument, 0, 0 },
        {"last", required_argument, 0, 0 },
        {"random", optional_argument, 0, 0 },
        {"list", no_argument, 0, 0 },
        {"silent", no_argument, 0, 0 },
        {"verbose", no_argument, 0, 0 },
        { 0, 0, 0, 0}
    };

    const char *short_options = "f:l:r:Lsv";

    struct numtest_args args = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    while(1) {
        int option_index = 0;
        int c = getopt_long(argc, argv, short_options, long_options, &option_index);

        if(c == -1)
            break;

        if((c == 0 && strcmp(long_options[option_index].name, "first") == 0) ||
            c == 'f') {
            if(!optarg || sscanf(optarg, "%lu", &args.first) != 1)
                usage();
        } else if((c == 0 && strcmp(long_options[option_index].name, "last") == 0) ||
            c == 'l') {
            if(!optarg || sscanf(optarg, "%lu", &args.last) != 1)
                usage();
        } else if((c == 0 && strcmp(long_options[option_index].name, "random") == 0) ||
            c == 'r') {
            args.random = 1;
            if(optarg && sscanf(optarg, "%u", &args.random_seed) != 1)
                usage();
        } else if((c == 0 && strcmp(long_options[option_index].name, "list") == 0) ||
            c == 'L') {
            args.list_tests = 1;
        } else if((c == 0 && strcmp(long_options[option_index].name, "silent") == 0) ||
            c == 's') {
            args.silent = 1;
        } else if((c == 0 && strcmp(long_options[option_index].name, "verbose") == 0) ||
            c == 'v') {
            args.verbose = 1;
        } else {
            usage();
        }
    }

    args.tests = argv + optind;
    args.num_tests = argc - optind;

    if(args.random && args.random_seed == 0)
        args.random_seed = time(NULL);

    return args;
}

int numtest_main(int argc, char *argv[]) {
    struct numtest_args args = parse_args(argc, argv);

    if(args.list_tests) {
        numtest_list_tests();
        return EXIT_SUCCESS;
    }

    return numtest_run_tests(&args) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#if 0
void dummy_test(double *params, int num_params, void *extra_args, struct numtest_ctx* test_ctx) {
    (void)extra_args;

    for(int i = 0; i < num_params; ++i)
        ASSERT_RANGEF(params[i], 0.0, 1.0, "Parameter %d not within range: %lf\n", i, params[i]);

    /*
    for(int i = 0; i < num_params; ++i)
        printf("%2.4f\t", params[i]);
    printf("\n");
    */
}

const struct numtest_case numtest_cases[] = {
    { "dummy_test", dummy_test, 3, 0 },
    { "dummy_test2", dummy_test, 2, 0 },
    { 0, 0, 0, 0 }
};

int main(int argc, char *argv[]) { return numtest_main(argc, argv); }

#endif
