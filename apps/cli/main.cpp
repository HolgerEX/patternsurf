#include <iostream>
#include <string_view>

#ifndef PATTERNSURF_VERSION
#define PATTERNSURF_VERSION "unknown"
#endif

namespace {

void print_usage(std::ostream& output)
{
    output << "Usage: pattern_cli [--version | --help]\n";
}

}  // namespace

int main(int argc, char* argv[])
{
    if (argc == 2 && std::string_view(argv[1]) == "--version") {
        std::cout << "PatternSurface " << PATTERNSURF_VERSION << '\n';
        return 0;
    }

    if (argc == 1 || (argc == 2 && std::string_view(argv[1]) == "--help")) {
        print_usage(std::cout);
        return 0;
    }

    print_usage(std::cerr);
    return 2;
}
