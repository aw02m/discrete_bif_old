#include "sys_common.hpp"

int main(int argc, char* argv[])
{
    if (argc != 2){
        std::cerr << "Put a path to the input json." << std::endl;
        std::exit(1);
    }

    std::ifstream ifs(argv[1]);
    if (ifs.fail()){
        std::cerr << "File does NOT exist." << std::endl;
        std::exit(1);
    }

    nlohmann::json json;
    ifs >> json;

    dynamical_system *ds = new dynamical_system(json);

    newton(ds);

    return 0;
}